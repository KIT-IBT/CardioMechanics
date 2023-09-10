/*
 *  CBModelExporter.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#include "CBModelExporter.h"
#include "CBModel.h"
#include <unistd.h>



CBModelExporter::CBModelExporter(CBModel* model, ParameterMap* parameters)
{
    model_ = model;
    model_->SetExporter(this);
    parameters_ = parameters;
    Init();
    writeToFileThread_ = 0;
}


CBModelExporter::~CBModelExporter()
{
    if(isInitialized_)
    {
        std::cout << "CBModelExporter::~CBModelExporter() was not deinitialzed correctly. This may cause wrong results !!!!\n";
        exit(1);
    }
    
}

void CBModelExporter::DeInit()
{
    if(!model_)
        throw std::runtime_error("void CBModelExporter::DeInit(): the model object is already destroyed ... very bad");
    if(writeToFileThread_)
    {
        { std::lock_guard<std::mutex>(model_->GetLock()); }
        writeToFileThread_->join();
        delete writeToFileThread_;
    }
    writeToFileThread_ = 0;
    isInitialized_ = false;
}

void CBModelExporter::WriteToFileThreadFunction()
{
    std::lock_guard<std::mutex>(model_->GetLock());
    ExportModel();
    ExportGlobalData();
}

void CBModelExporter::ExportGlobalData()
{
    for(auto it : model_->GetGlobalData())
    {
        std::ofstream file;
        
        if(counter_ == 1)
            file.open((exportGlobalDataFileDir_+"/" + it.first + ".dat"));
        else
            file.open((exportGlobalDataFileDir_+"/" + it.first + ".dat"), std::ios::app);
        
        file << model_->GetCurrentTime() << " " << it.second << "\n";
        file.close();
    }
}


bool CBModelExporter::GetExportOption(std::string str, bool defaultValue)
{
    if(!parameters_)
        throw std::runtime_error("bool CBModelExporter::GetExportOption(std::string str,bool defaultValue = false): No parameters object given");
    else
        return(parameters_->Get<bool>("Export.Options." + str, defaultValue));
}


void CBModelExporter::Init()
{
    counter_                  = 0;
    if(parameters_)
    {
        exportTimeStep_ = parameters_->Get<TFloat>("Export.TimeStep", 1e-3);
        std::string str = parameters_->Get<std::string>("Export.Prefix", "");
        if(str == "")
            throw std::runtime_error("CBModelExporter::Init(): No export directory and prefix given !!!");
        size_t pos = str.find_last_of("/\\");
        if(pos == std::string::npos)
        {
            exportFileDir_        = std::string("./");
            exportFilenamePrefix_ = str;
        }
        else
        {
            exportFileDir_        = str.substr(0, pos);
            exportFilenamePrefix_ = str.substr(pos+1);
        }
    }
    else
    {
        throw std::runtime_error("CBModelExporter::Init(): parameters_ not set");
    }
    isInitialized_ = true;
    exportGlobalDataFileDir_ = exportFileDir_ + "/" + exportFilenamePrefix_ + "_data";
    if(!frizzle::filesystem::CreateDirectory(exportGlobalDataFileDir_))
        throw std::runtime_error("void CBModelExporter::InitPvdFile(): Path: " + exportGlobalDataFileDir_ + " exists but is not a directory");
}


void CBModelExporter::WriteToFile()
{
    if(exportFilenamePrefix_ != "")
    {
        if(writeToFileThread_ != 0)
        {
            writeToFileThread_->join();
            delete writeToFileThread_;
            writeToFileThread_ = 0;
        }
        writeToFileThread_ = new std::thread(std::bind(&CBModelExporter::WriteToFileThreadFunction, this));
    }
}

