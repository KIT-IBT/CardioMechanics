/*
 *  CBCircWindkessel3AsymptoticPressure.cpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 26.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCircWindkessel3AsymptoticPressure.h"

CBCircWindkessel3AsymptoticPressure::CBCircWindkessel3AsymptoticPressure(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "VentrVolume";
    stateVarsNames_[1] = "ArtVolume";
    
    //----- Results -----
    
    results_[0] = &ventrPressure_;  resultsNames_[0] = "VentrPressure";
    results_[1] = &artPressure_;    resultsNames_[1] = "ArtPressure";
    results_[2] = &artFlow_;        resultsNames_[2] = "ArtFlow";
    results_[3] = &perFlow_;        resultsNames_[3] = "PerFlow";
    results_[4] = &venFlow_;        resultsNames_[4] = "VenFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &ventrPressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);
    
    artResist_          = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtResist");
    artCompli_          = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtCompli");
    perResist_          = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PerResist");
    asymptoticPressure_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.AsymptoticPressure");
    venResist_          = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenResist");
    venPressure_        = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenPressure");
    
    stateVars_[1]   = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.ArtVolume");
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.VentrSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.VentrPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.VentrPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, -1);
    
    WriteHeaderToFile();
}


void CBCircWindkessel3AsymptoticPressure::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
}


void CBCircWindkessel3AsymptoticPressure::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
}


void CBCircWindkessel3AsymptoticPressure::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat ventrPressure = cavityPressures[0];
    
    TFloat artVolume = stateVars[1];
    
    TFloat artCompliPressure = artVolume/artCompli_;
    
    if(ventrPressure >= artCompliPressure)
    {
        artFlow_ = (ventrPressure-artCompliPressure)/artResist_;
        artPressure_ = ventrPressure;
    }
    else
    {
        artFlow_ = 0.0;
        artPressure_ = artCompliPressure;
    }
    
    perFlow_ = (artCompliPressure-asymptoticPressure_)/perResist_;
    
    if(venPressure_ >= ventrPressure)
        venFlow_ = (venPressure_-ventrPressure)/venResist_;
    else
        venFlow_ = 0.0;
}


void CBCircWindkessel3AsymptoticPressure::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] = venFlow_ - artFlow_;
    stateVarsDot[1] = artFlow_ - perFlow_;
}
