/*
 * File: CBDataFromFunction.cpp
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


#include "CBDataFromFunction.h"

#include "CBAnalyticFunction.h"


CBDataFromFunction::~CBDataFromFunction() {
    delete getImpl_;
}


CBDataFromFunction::CBDataFromFunction(ParameterMap* parameters, std::string parameterKey)
{
    std::string type = parameters->Get<std::string>(parameterKey + ".Type", "Linear");
    offset_ = parameters->Get<TFloat>(parameterKey + ".Offset", 0.);
    
    // note: this could be better solved using a factory class for analytic
    // functions, see e.g. CBElementFactory
    if(type == "Linear")
        getImpl_ = new Linear(parameters, parameterKey);
    else if(type == "Sinus")
        getImpl_ = new Sinus(parameters, parameterKey);
    else if(type == "DoubleHill")
        getImpl_ = new DoubleHill(parameters, parameterKey);
    else if(type == "SinusDriver")
        getImpl_ = new SinusDriver(parameters, parameterKey);
    else
        throw std::runtime_error("CBAnalyticFunction: Type " + type + " is unkown");
}


TFloat CBDataFromFunction::Get(TFloat time, TInt index)
{
    time = time-offset_;
    return((*getImpl_)(time));
};
