/*
 *  CBCircTwoVentricles.cpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 29.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCircTwoVentricles.h"

CBCircTwoVentricles::CBCircTwoVentricles(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "LvVolume";
    stateVarsNames_[1] = "SysArtVolume";
    stateVarsNames_[2] = "SysVenVolume";
    stateVarsNames_[3] = "RvVolume";
    stateVarsNames_[4] = "PulArtVolume";
    stateVarsNames_[5] = "PulVenVolume";
    stateVarsNames_[6] = "StrokeVolumeDiff";
    
    //----- Results -----
    
    results_[ 0] = &lvPressure_;        resultsNames_[ 0] = "LvPressure";
    results_[ 1] = &sysArtPressure_;    resultsNames_[ 1] = "SysArtPressure";
    results_[ 2] = &sysVenPressure_;    resultsNames_[ 2] = "SysVenPressure";
    results_[ 3] = &rvPressure_;        resultsNames_[ 3] = "RvPressure";
    results_[ 4] = &pulArtPressure_;    resultsNames_[ 4] = "PulArtPressure";
    results_[ 5] = &pulVenPressure_;    resultsNames_[ 5] = "PulVenPressure";
    results_[ 6] = &sysArtFlow_;        resultsNames_[ 6] = "SysArtFlow";
    results_[ 7] = &sysPerFlow_;        resultsNames_[ 7] = "SysPerFlow";
    results_[ 8] = &sysVenFlow_;        resultsNames_[ 8] = "SysVenFlow";
    results_[ 9] = &pulArtFlow_;        resultsNames_[ 9] = "PulArtFlow";
    results_[10] = &pulPerFlow_;        resultsNames_[10] = "PulPerFlow";
    results_[11] = &pulVenFlow_;        resultsNames_[11] = "PulVenFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &lvPressure_;
    estCavityPressures_[1] = &rvPressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);
    
    sysArtValveResist_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveResist");
    sysArtResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtResist");
    sysArtCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtCompli");
    sysArtVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtVolumeUnstr");
    sysPerResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysPerResist");
    sysVenResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenResist");
    sysVenCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenCompli");
    sysVenVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenVolumeUnstr");
    pulArtValveResist_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveResist");
    pulArtResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtResist");
    pulArtCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtCompli");
    pulArtVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtVolumeUnstr");
    pulPerResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulPerResist");
    pulVenResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenResist");
    pulVenCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenCompli");
    pulVenVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenVolumeUnstr");
    
    totalVolume_  = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtVolume");
    stateVars_[4] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtVolume");
    stateVars_[5] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulVenVolume");
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.LvSurface");
    cavitySurfaceIndices_[1] = parameters->Get<TInt>(parameterPrefix + ".Cavities.RvSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.LvPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.RvPreloading", 0.0);
    cavityPreloadingPressures_[2] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LvPressure", 0.0);
    cavityPreloadingPressures_[3] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RvPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, 6);
    
    WriteHeaderToFile();
}


void CBCircTwoVentricles::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[3] = cavityVolumes[1];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[1]-stateVars_[3]-stateVars_[4]-stateVars_[5];
}


void CBCircTwoVentricles::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
    cavityVolumes[1] = estStateVars_[3];
}


void CBCircTwoVentricles::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat lvPressure = cavityPressures[0];
    TFloat rvPressure = cavityPressures[1];
    
    TFloat sysArtVolume = stateVars[1];
    TFloat sysVenVolume = stateVars[2];
    TFloat pulArtVolume = stateVars[4];
    TFloat pulVenVolume = stateVars[5];
    
    TFloat sysArtCompliPressure = (sysArtVolume-sysArtVolumeUnstr_)/sysArtCompli_;
    sysVenPressure_ = (sysVenVolume-sysVenVolumeUnstr_)/sysVenCompli_;
    TFloat pulArtCompliPressure = (pulArtVolume-pulArtVolumeUnstr_)/pulArtCompli_;
    pulVenPressure_ = (pulVenVolume-pulVenVolumeUnstr_)/pulVenCompli_;
    
    if(lvPressure >= sysArtCompliPressure)
        sysArtFlow_ = (lvPressure-sysArtCompliPressure)/(sysArtValveResist_+sysArtResist_);
    else
        sysArtFlow_ = 0.0;
    sysArtPressure_ = sysArtCompliPressure + sysArtFlow_ * sysArtResist_;
    
    sysPerFlow_ = (sysArtCompliPressure-sysVenPressure_)/sysPerResist_;
    
    if(sysVenPressure_ >= rvPressure)
        sysVenFlow_ = (sysVenPressure_-rvPressure)/sysVenResist_;
    else
        sysVenFlow_ = 0.0;
    
    if(rvPressure >= pulArtCompliPressure)
        pulArtFlow_ = (rvPressure-pulArtCompliPressure)/(pulArtValveResist_+pulArtResist_);
    else
        pulArtFlow_ = 0.0;
    pulArtPressure_ = pulArtCompliPressure + pulArtFlow_ * pulArtResist_;
    
    pulPerFlow_ = (pulArtCompliPressure-pulVenPressure_)/pulPerResist_;
    
    if(pulVenPressure_ >= lvPressure)
        pulVenFlow_ = (pulVenPressure_-lvPressure)/pulVenResist_;
    else
        pulVenFlow_ = 0.0;
}


void CBCircTwoVentricles::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] = pulVenFlow_ - sysArtFlow_;
    stateVarsDot[1] = sysArtFlow_ - sysPerFlow_;
    stateVarsDot[2] = sysPerFlow_ - sysVenFlow_;
    stateVarsDot[3] = sysVenFlow_ - pulArtFlow_;
    stateVarsDot[4] = pulArtFlow_ - pulPerFlow_;
    stateVarsDot[5] = pulPerFlow_ - pulVenFlow_;
    stateVarsDot[6] = sysArtFlow_ - pulArtFlow_;
}
