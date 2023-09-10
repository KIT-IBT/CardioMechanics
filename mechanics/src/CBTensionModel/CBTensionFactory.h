/*
 *  CBTensionFactory.h
 *  CardioMechanics
 *
 *  Created by Lukas Baron on Thu Apr 27 2017.
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#pragma once

#include "CBTensionModel.h"
#include "CBElementSolid.h"
class CBElementSolid;

// ===== declaration =====

class CBTensionFactory {
protected:
    using ProducerFunction = std::function<CBTensionModel *(CBElementSolid *)>;
    std::map<const std::string, ProducerFunction> producers_;
    
    void Register(std::string tensionName, ProducerFunction function) {
        producers_[tensionName] = function;
    }
    
    CBFileManager *fileManager_ = 0;
    ParameterMap *parameters_ = 0;
    CBTiming *timing_ = 0;
    bool isInitialized_ = false;
    
public:
    CBTensionFactory();
    ~CBTensionFactory() {}
    
    void Init(ParameterMap *parameters,  CBTiming *timing,  CBFileManager *fileManager);
    CBTensionModel *New(CBElementSolid *ele);
};
