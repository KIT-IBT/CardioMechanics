/*
 *  CBConstitutiveModel.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#include "CBConstitutiveModel.h"

double relaxedExp(double a, double criticalValue)
{
    double r = exp(a);
    
    if(r < criticalValue)
        return r;
    else
        return criticalValue + criticalValue * ( a - log(criticalValue));
}

double relaxedExpDerivative(double a, double criticalValue)
{
    double r = exp(a);
    
    if(r < criticalValue)
        return a*r;
    else
        return criticalValue;
}

void CBConstitutiveModel::Init(ParameterMap* parameters, TInt materialIndex)
{
    std::stringstream mi;
    mi << materialIndex;
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".IgnoreCorruptElements") == false && parameters->IsAvailable("Materials.Mat_Default.IgnoreCorruptElements") == true)
        ignoreCorruptElements_ =  parameters->Get<bool>("Materials.Mat_Default.IgnoreCorruptElements",false);
    else
        ignoreCorruptElements_ =  parameters->Get<bool>("Materials.Mat_" + mi.str() + ".IgnoreCorruptElements",false);
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".CriticalVolumeChange") == false && parameters->IsAvailable("Materials.Mat_Default.CriticalVolumeChange") == true)
        criticalVolumeChange_ =  parameters->Get<TFloat>("Materials.Mat_Default.CriticalVolumeChange",0.1);
    else
        criticalVolumeChange_ =  parameters->Get<TFloat>("Materials.Mat_" + mi.str() + ".CriticalVolumeChange",0.1);
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".TensionMax") == false && parameters->IsAvailable("Materials.Mat_Default.TensionMax") == true)
        Tmax_ =  parameters->Get<double>("Materials.Mat_Default.TensionMax");
    else
        Tmax_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".TensionMax", 0.);
    
}

