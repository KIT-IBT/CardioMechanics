/*
 * File: CBTensionFactory.h
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
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
