/*
 * File: CBTensionFromFunction.h
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
