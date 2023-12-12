/*
 * File: CBCircOneVentricle.cpp
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


#include "CBCircOneVentricle.h"

CBCircOneVentricle::CBCircOneVentricle(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "VentrVolume";
    stateVarsNames_[1] = "ArtVolume";
    stateVarsNames_[2] = "VenVolume";
    
    //----- Results -----
    
    results_[0] = &ventrPressure_;  resultsNames_[0] = "VentrPressure";
    results_[1] = &artPressure_;    resultsNames_[1] = "ArtPressure";
    results_[2] = &venPressure_;    resultsNames_[2] = "VenPressure";
    results_[3] = &artFlow_;        resultsNames_[3] = "ArtFlow";
    results_[4] = &perFlow_;        resultsNames_[4] = "PerFlow";
    results_[5] = &venFlow_;        resultsNames_[5] = "VenFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &ventrPressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);
    
    artValveResist_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtValveResist");
    artResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtResist");
    artCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtCompli");
    artVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtVolumeUnstr");
    perResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PerResist");
    venResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenResist");
    venCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenCompli");
    venVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenVolumeUnstr");
    
    totalVolume_  = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.ArtVolume");
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.VentrSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.VentrPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.VentrPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, -1);
    
    WriteHeaderToFile();
}


void CBCircOneVentricle::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[1];
}


void CBCircOneVentricle::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
}


void CBCircOneVentricle::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat ventrPressure = cavityPressures[0];
    
    TFloat artVolume = stateVars[1];
    TFloat venVolume = stateVars[2];
    
    TFloat artCompliPressure = (artVolume-artVolumeUnstr_)/artCompli_;
    venPressure_ = (venVolume-venVolumeUnstr_)/venCompli_;
    
    if(ventrPressure >= artCompliPressure)
        artFlow_ = (ventrPressure-artCompliPressure)/(artValveResist_+artResist_);
    else
        artFlow_ = 0.0;
    artPressure_ = artCompliPressure + artFlow_ * artResist_;
    
    perFlow_ = (artCompliPressure-venPressure_)/perResist_;
    
    if(venPressure_ >= ventrPressure)
        venFlow_ = (venPressure_-ventrPressure)/venResist_;
    else
        venFlow_ = 0.0;
}


void CBCircOneVentricle::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] = venFlow_ - artFlow_;
    stateVarsDot[1] = artFlow_ - perFlow_;
    stateVarsDot[2] = perFlow_ - venFlow_;
}
