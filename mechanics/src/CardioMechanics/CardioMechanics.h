/*
 *  CardioMechanics.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 19.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CARDIO_MECHANICS
#define CARDIO_MECHANICS

#include <string>
#include <sstream>
#include <vector>

#include "ParameterMap.h"

#include "CBModel.h"
#include "CBModelLoader.h"
#include "CBModelLoaderTetgen.h"
#include "CBSolver.h"
#include "CBModelExporterVTK.h"


using namespace math_pack;

class CardioMechanics
{
public:
    CardioMechanics(){}
    ~CardioMechanics()
    {
        DeInit();
    }
    
    void Init1();
    void Init2();
    void DeInit();
    void Run();
    
    void ReadParameterFile(std::string parameterFile);
    void SetParameter(std::string key, std::string value);
    void PrintParameters();
    
protected:
private:
    void InitModel();
    void InitSolver();
    void InitModelExporter();
    void InitStdoutHandler();
    void StdoutHandler();
    void DomainDecomposition();
    
    CBModel*         model_         = 0;
    CBModelExporter* modelExporter_ = 0;
    ParameterMap*    parameters_    = 0;
    CBSolver*        solver_        = 0;
    std::string      logFilename_;
    std::string      parameterFile_;
};
#endif
