/*
 *  CardioMechanics.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "filesystem.h"

#include "CardioMechanics.h"
#include "CBSolver.h"
#include "CBSolverNewmarkBeta.h"
#include "CBSolverEquilibrium.h"


void CardioMechanics::Init1() {
    logFilename_ = parameters_->Get<std::string>("General.LogFile", "");
    DCCtrl::SetLogFile(logFilename_);
    
    if(parameters_->Get<bool>("General.Verbose", false))
        DCCtrl::VerboseOn();
    if(parameters_->Get<bool>("General.Debug", false))
        DCCtrl::DebugOn();
    
    DCCtrl::print << "\n\n\n----------------------------- Loading Setting ... ------------------------------\n\n";
    
    if(logFilename_!="")
        DCCtrl::print << "\tLog file: " << logFilename_ << "\n\n";
    
    DCCtrl::print << "\tSetting: " << parameters_->Get<std::string>("Setting","") << "\n";
    DCCtrl::print << "\tParameter file: " << parameterFile_ << "\n";
    DCCtrl::print << "\n\tLoading model...";
    
    InitModel();
    
    DCCtrl::print << "\t\tAmount of nodes: " << model_->GetNodes().size() << " \n";
    DCCtrl::print << "\t\tAmount of elements: " << model_->GetElements().size() << " \n\n";
    
    InitModelExporter();
}

void CardioMechanics::Init2()
{
    InitSolver();
    
    if(parameters_->Get<bool>("Solver.DomainDecomposition",false))
    {
        DomainDecomposition();
        InitSolver(); // Solver must be reinitialized with new node indexing;
    }
    
    solver_->Init(parameters_, model_);
}


void CardioMechanics::DeInit()
{
    if(modelExporter_)
        modelExporter_->DeInit(); // Very important since the exporter thread might still be active;
    if(model_)
        delete model_;
    if(solver_)
        delete solver_;
    if(parameters_)
        delete parameters_;
    
}

void CardioMechanics::Run()
{
    DCCtrl::print << "\n\n       Solver: " << solver_->GetType() << "\n\n";
    DCCtrl::print << "       Writing results to: " << parameters_->Get<std::string>("Export.Prefix","") << ".[timestep].vtu\n\n";
    
    solver_->Run();
    
}

void CardioMechanics::ReadParameterFile(std::string parameterFile)
{
    parameterFile_ = parameterFile;
    
    if(parameters_)
        delete parameters_;
    parameters_ = new ParameterMap;
    parameters_->ReadXml(parameterFile);
}

void CardioMechanics::SetParameter(std::string key, std::string value)
{
    if(parameters_ != 0)
        parameters_->Set(key, value);
}

void CardioMechanics::PrintParameters()
{
    DCCtrl::verbose<< "\n\n-----  Requested Parameters   ----\n";
    std::vector<std::string> s = parameters_->GetRequestedParameters();
    for(auto i:s)
        DCCtrl::cverbose << i <<"\n";
    DCCtrl::verbose << std::flush;
    
    DCCtrl::verbose << "-----  Ignored Parameters   ----\n";
    s = parameters_->GetIgnoredParameters();
    for(auto i: s)
        DCCtrl::cverbose << i <<"\n";
    DCCtrl::verbose << std::flush;
    DCCtrl::verbose<< "----- ----- ----- ----- ----- ----\n";
}

void CardioMechanics::InitModelExporter()
{
    if(parameters_ == 0 || model_ == 0)
        throw std::runtime_error("CardioMechanics::InitModelExporter(): CardioMechanics::ReadParameterFile(std::string parameterFile) and  CardioMechanics::InitModel() have to be run first !!!!");
    else
    {
        std::string modelExporterFormat = parameters_->Get<std::string>("Export.Format", "VTK");
        if(modelExporterFormat == "VTK")
        {
            modelExporter_ = new CBModelExporterVTK(model_, parameters_);
        }
        else
            throw std::runtime_error("CardioMechanics::InitModelExporter(): Output format " + modelExporterFormat + " unkown !!!");
    }
}

void CardioMechanics::InitModel()
{
    CBModelLoader* modelLoader;
    std::string                        format = parameters_->Get<std::string>("Mesh.Format");
    
    if(format == "Tetgen")
        modelLoader = new CBModelLoaderTetgen(parameters_);
    else
        throw std::runtime_error("Mesh Format: " + format + " is unkown !");
    model_ = modelLoader->Load();
    
    std::string sorting = parameters_->Get<std::string>("Mesh.Sorting", "None");
    
    model_->InitMapping();
    
    if(sorting != "None")
    {
        if(sorting == "PCA")
            model_->ApplySpatialSortPCA();
        else
            throw std::runtime_error("Sorting Method: " + sorting + " is unkown !");
    }
    delete modelLoader;
}

void CardioMechanics::InitSolver()
{
    if(solver_)
        delete solver_;
    solver_ = 0;
    
    std::string solverType = parameters_->Get<std::string>("Solver.Type");
    
    if(solverType == "Static")
        solver_ = new CBSolverEquilibrium();
    else if(solverType == "NewmarkBeta")
        solver_ = new CBSolverNewmarkBeta();
    
    if(!solver_)
        throw std::runtime_error("CardioMechanics::InitSolver(): Solver type: " + solverType + " is unkown ");
}

void CardioMechanics::DomainDecomposition()
{
    std::vector<std::pair<int,int> >mapping;
    
    solver_->GetDomainDecomposition(parameters_,model_,mapping);
    model_->ApplyDomainDecomposition(mapping);
    
    MPI_Barrier(Petsc::Comm());
}


