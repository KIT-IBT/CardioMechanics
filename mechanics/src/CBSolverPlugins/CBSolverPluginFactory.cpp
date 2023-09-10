/*
 *  CBSolverPluginFactory.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBSolverPlugin.h"
#include "CBSolverPluginFactory.h"

#include "CBContactHandling.h"
#include "CBApplyPressureFromFunction.h"
#include "CBApplyPressureFromFunctionNodeExport.h"
#include "CBCirculation.h"
#include "CBApplyPressure.h"
#include "CBReferenceRecovery.h"
#include "CBLoadUnloadedState.h"
#include "CBRobinBoundary.h"
#include "CBRobinBoundaryGeneral.h"
#include "CBacCELLerate.h"

class CBSolverPlugin;

void CBSolverPluginFactory::LoadAllPlugins(std::vector<CBSolverPlugin *> &plugins, ParameterMap *parameters) {
    
    if(parameters->Get<bool>("Solver.Plugins.acCELLerate", false)) {
        CBSolverPlugin *solverPlugin = new CBacCELLerate();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.LoadUnloadedState", false)) {
        CBSolverPlugin *solverPlugin = new CBLoadUnloadedState();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.ReferenceRecovery", false)) {
        CBSolverPlugin *solverPlugin = new CBReferenceRecovery();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.Circulation", false)) {
        CBSolverPlugin *solverPlugin = new CBCirculation();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.ContactHandling", false)) {
        CBSolverPlugin *solverPlugin = new CBContactHandling();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.RobinBoundary", false)) {
        CBSolverPlugin *solverPlugin = new CBRobinBoundary();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.RobinBoundaryGeneral", false)) {
        CBSolverPlugin *solverPlugin = new CBRobinBoundaryGeneral();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.ApplyPressureFromFunction", false)) {
        CBSolverPlugin *solverPlugin = new CBApplyPressureFromFunction();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.ApplyPressureFromFunctionNodeExport", false)) {
        CBSolverPlugin *solverPlugin = new CBApplyPressureFromFunctionNodeExport();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    if (parameters->Get<bool>("Solver.Plugins.ApplyPressure", false)) {
        CBSolverPlugin *solverPlugin = new CBApplyPressure();
        solverPlugin->SetParameters(parameters);
        plugins_.push_back(solverPlugin);
    }
    
    plugins = plugins_;
} // CBSolverPluginFactory::LoadAllPlugins
