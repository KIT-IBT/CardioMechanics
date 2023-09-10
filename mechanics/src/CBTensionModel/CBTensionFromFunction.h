/*
 *  CBDTensionFromFunction.h
 *  CardioMechanics
 *
 *  Created by Armin Mueller on 09. August 2017.
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#pragma once

#include "CBTensionModel.h"
#include "ParameterMap.h"


class CBAnalyticFunction;

///Tension model were the force is calculated by different functions (Linear, Sinus, DoubleHill, SinusDriver)
class CBTensionFromFunction : public CBTensionModel {
    
public:
    CBTensionFromFunction(ParameterMap* parameters, CBElement* e);
    ~CBTensionFromFunction();
    virtual double CalcActiveTension(const math_pack::Matrix3<double>& deformation, const double time) override;
    
    
protected:
    
private:
    CBAnalyticFunction* getImpl_ = 0;
    std::string GetWithKeyOrWithFallback(ParameterMap* parameters, std::string key, std::string keyFallback, std::string tail, std::string defaultValue);
};
