/*
 *  CBCircIsovolumetric.cpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 29.01.20.
 *  Copyright 2020 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCircIsovolumetric.h"

CBCircIsovolumetric::CBCircIsovolumetric(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "DummyStateVar";
    
    //----- Results -----
    
    results_[0] = &pressure_; resultsNames_[0] = "Pressure";
    results_[1] = &volume_; resultsNames_[1] = "Volume";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &pressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.Surface");
    
    cavityPreloadingPressures_[0] = 0.0;
    cavityPreloadingPressures_[1] = 0.0;
    
    WriteHeaderToFile();
}


void CBCircIsovolumetric::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    initialCavityVolume_ = cavityVolumes[0];
    volume_ = initialCavityVolume_;
}


void CBCircIsovolumetric::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    *(estCavityPressures_[0]) = cavityPressures[0];
    
    cavityVolumes[0] = initialCavityVolume_;
}


void CBCircIsovolumetric::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
}


void CBCircIsovolumetric::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
}
