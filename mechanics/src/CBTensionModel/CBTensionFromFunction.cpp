/*
 *  CBTensionFromFunction.h
 *  CardioMechanics
 *
 *  Created by Armin Mueller on 09. August 2017
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#include "CBTensionFromFunction.h"
#include "CBAnalyticFunction.h"


CBTensionFromFunction::~CBTensionFromFunction() {
    delete getImpl_;
}

/// Method to get a value from xml with a fallback and a default value. Priorities: key+tail > keyFallback+tail > defaultValue
std::string CBTensionFromFunction::GetWithKeyOrWithFallback(ParameterMap* parameters, std::string key, std::string keyFallback, std::string tail, std::string defaultValue) {
    
    if(!(parameters->IsAvailable(key + tail)) && (parameters->IsAvailable(keyFallback + tail))) {
        //key is not available, but the fallback
        return parameters->Get<std::string>(keyFallback);
    } else {
        //key is available, so use it (or the default value)
        return parameters->Get<std::string>(key+tail, defaultValue);
    }
}

CBTensionFromFunction::CBTensionFromFunction(ParameterMap* parameters, CBElement* e) {
    //Get Material Index and use it (or the Mat_Default)
    int matIndex = e->GetMaterialIndex();
    std::string key = "Materials.Mat_" + std::to_string(matIndex);
    std::string keyFallback = "Materials.Mat_Default";
    //if type is not given, this will throw an exception
    std::string type = GetWithKeyOrWithFallback(parameters, key, keyFallback, ".FromFunction.Type", "None");
    // max tension, variable is defined in CBTensionModel
    Tmax_ = parameters->Get<TFloat>(key+".TensionMax");
    
    if(type == "Linear")
        getImpl_ = new Linear(parameters, key+".FromFunction");
    else if(type == "Sinus")
        getImpl_ = new Sinus(parameters, key+".FromFunction");
    else if(type == "DoubleHill")
        getImpl_ = new DoubleHill(parameters, key+".FromFunction");
    else if(type == "SinusDriver")
        getImpl_ = new SinusDriver(parameters, key+".FromFunction");
    else
        throw std::runtime_error("CBTensionFromFunction: Type " + type + " is unkown");
}


double CBTensionFromFunction::CalcActiveTension(const math_pack::Matrix3<double>& deformation, const double time) {
    return(Tmax_ * (*getImpl_)(time - GetActivationTime()));
}
