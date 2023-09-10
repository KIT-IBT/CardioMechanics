/*
 *  CBCircWholeHeartValves.cpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 10.03.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCircWholeHeartValves.h"

CBCircWholeHeartValves::CBCircWholeHeartValves(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[ 0] = "LvVolume";
    stateVarsNames_[ 1] = "SysArtVolume";
    stateVarsNames_[ 2] = "SysVenVolume";
    stateVarsNames_[ 3] = "RaVolume";
    stateVarsNames_[ 4] = "RvVolume";
    stateVarsNames_[ 5] = "PulArtVolume";
    stateVarsNames_[ 6] = "PulVenVolume";
    stateVarsNames_[ 7] = "LaVolume";
    stateVarsNames_[ 8] = "SysArtFlow";
    stateVarsNames_[ 9] = "RavFlow";
    stateVarsNames_[10] = "PulArtFlow";
    stateVarsNames_[11] = "LavFlow";
    stateVarsNames_[12] = "StrokeVolumeDiff";
    
    //----- Results -----
    
    results_[ 0] = &lvPressure_;        resultsNames_[ 0] = "LvPressure";
    results_[ 1] = &sysArtPressure_;    resultsNames_[ 1] = "SysArtPressure";
    results_[ 2] = &sysVenPressure_;    resultsNames_[ 2] = "SysVenPressure";
    results_[ 3] = &raPressure_;        resultsNames_[ 3] = "RaPressure";
    results_[ 4] = &rvPressure_;        resultsNames_[ 4] = "RvPressure";
    results_[ 5] = &pulArtPressure_;    resultsNames_[ 5] = "PulArtPressure";
    results_[ 6] = &pulVenPressure_;    resultsNames_[ 6] = "PulVenPressure";
    results_[ 7] = &laPressure_;        resultsNames_[ 7] = "LaPressure";
    results_[ 8] = &sysPerFlow_;        resultsNames_[ 8] = "SysPerFlow";
    results_[ 9] = &sysVenFlow_;        resultsNames_[ 9] = "SysVenFlow";
    results_[10] = &pulPerFlow_;        resultsNames_[10] = "PulPerFlow";
    results_[11] = &pulVenFlow_;        resultsNames_[11] = "PulVenFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &lvPressure_;
    estCavityPressures_[1] = &raPressure_;
    estCavityPressures_[2] = &rvPressure_;
    estCavityPressures_[3] = &laPressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);    
    
    bloodDensity_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.BloodDensity");
    
    sysArtValveMax_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveMax");
    sysArtValveMin_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveMin");
    sysArtValveAreaRef_     = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveAreaRef");
    sysArtResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtResist");
    sysArtCompli_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtCompli");
    sysArtVolumeUnstr_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtVolumeUnstr");
    sysPerResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysPerResist");
    sysVenResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenResist");
    sysVenCompli_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenCompli");
    sysVenVolumeUnstr_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysVenVolumeUnstr");
    ravValveMax_            = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveMax");
    ravValveMin_            = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveMin");
    ravValveAreaRef_        = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveAreaRef");
    pulArtValveMax_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveMax");
    pulArtValveMin_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveMin");
    pulArtValveAreaRef_     = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveAreaRef");
    pulArtResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtResist");
    pulArtCompli_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtCompli");
    pulArtVolumeUnstr_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtVolumeUnstr");
    pulPerResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulPerResist");
    pulVenResist_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenResist");
    pulVenCompli_           = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenCompli");
    pulVenVolumeUnstr_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenVolumeUnstr");
    lavValveMax_            = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveMax");
    lavValveMin_            = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveMin");
    lavValveAreaRef_        = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveAreaRef");
    
    totalVolume_   = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[ 1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtVolume");
    stateVars_[ 5] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtVolume");
    stateVars_[ 6] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulVenVolume");
    stateVars_[ 8] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtFlow", 0.0);
    stateVars_[ 9] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RavFlow", 0.0);
    stateVars_[10] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtFlow", 0.0);
    stateVars_[11] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LavFlow", 0.0);
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.LvSurface");
    cavitySurfaceIndices_[1] = parameters->Get<TInt>(parameterPrefix + ".Cavities.RaSurface");
    cavitySurfaceIndices_[2] = parameters->Get<TInt>(parameterPrefix + ".Cavities.RvSurface");
    cavitySurfaceIndices_[3] = parameters->Get<TInt>(parameterPrefix + ".Cavities.LaSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.LvPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.RaPreloading", 0.0);
    cavityPreloadingPressures_[2] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.RvPreloading", 0.0);
    cavityPreloadingPressures_[3] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.LaPreloading", 0.0);
    cavityPreloadingPressures_[4] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LvPressure", 0.0);
    cavityPreloadingPressures_[5] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RaPressure", 0.0);
    cavityPreloadingPressures_[6] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RvPressure", 0.0);
    cavityPreloadingPressures_[7] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LaPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, 12);
    
    WriteHeaderToFile();
}


void CBCircWholeHeartValves::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[3] = cavityVolumes[1];
    stateVars_[4] = cavityVolumes[2];
    stateVars_[7] = cavityVolumes[3];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[1]-stateVars_[3]-stateVars_[4]-stateVars_[5]-stateVars_[6]-stateVars_[7];
}


void CBCircWholeHeartValves::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
    cavityVolumes[1] = estStateVars_[3];
    cavityVolumes[2] = estStateVars_[4];
    cavityVolumes[3] = estStateVars_[7];
}


void CBCircWholeHeartValves::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat raPressure = cavityPressures[1];
    TFloat laPressure = cavityPressures[3];
    
    TFloat sysArtVolume = stateVars[ 1];
    TFloat sysVenVolume = stateVars[ 2];
    TFloat pulArtVolume = stateVars[ 5];
    TFloat pulVenVolume = stateVars[ 6];
    TFloat sysArtFlow   = stateVars[ 8];
    TFloat pulArtFlow   = stateVars[10];
    
    TFloat sysArtCompliPressure = (sysArtVolume-sysArtVolumeUnstr_)/sysArtCompli_;
    sysArtPressure_ = sysArtCompliPressure + sysArtFlow * sysArtResist_;
    sysVenPressure_ = (sysVenVolume-sysVenVolumeUnstr_)/sysVenCompli_;
    
    TFloat pulArtCompliPressure = (pulArtVolume-pulArtVolumeUnstr_)/pulArtCompli_;
    pulArtPressure_ = pulArtCompliPressure + pulArtFlow * pulArtResist_;
    pulVenPressure_ = (pulVenVolume-pulVenVolumeUnstr_)/pulVenCompli_;
    
    sysPerFlow_ = (sysArtCompliPressure-sysVenPressure_)/sysPerResist_;
    sysVenFlow_ = (sysVenPressure_-raPressure)/sysVenResist_;
    
    pulPerFlow_ = (pulArtCompliPressure-pulVenPressure_)/pulPerResist_;
    pulVenFlow_ = (pulVenPressure_-laPressure)/pulVenResist_;
}


void CBCircWholeHeartValves::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    TFloat lvPressure = cavityPressures[0];
    TFloat raPressure = cavityPressures[1];
    TFloat rvPressure = cavityPressures[2];
    TFloat laPressure = cavityPressures[3];
    
    Algebraics(cavityPressures, stateVars);
    
    TFloat sysArtVolume     = stateVars[ 1];
    TFloat pulArtVolume     = stateVars[ 5];
    TFloat sysArtFlow       = stateVars[ 8];
    TFloat ravFlow          = stateVars[ 9];
    TFloat pulArtFlow       = stateVars[10];
    TFloat lavFlow          = stateVars[11];
    
    TFloat sysArtCompliPressure = (sysArtVolume-sysArtVolumeUnstr_)/sysArtCompli_;
    TFloat pulArtCompliPressure = (pulArtVolume-pulArtVolumeUnstr_)/pulArtCompli_;
    
    //---------- SysArtValve ----------
    
    TFloat sysArtValveAreaEff;
    if(lvPressure > sysArtCompliPressure || sysArtFlow > 0)
    {
        TFloat sysArtValveAreaMax = sysArtValveMax_ * sysArtValveAreaRef_;
        sysArtValveAreaEff = sysArtValveAreaMax*sysArtValveAreaRef_ / (sysArtValveAreaRef_-sysArtValveAreaMax);
    }
    else
    {
        TFloat sysArtValveAreaMin = sysArtValveMin_ * sysArtValveAreaRef_;
        sysArtValveAreaEff = sysArtValveAreaMin*sysArtValveAreaRef_/(sysArtValveAreaRef_-sysArtValveAreaMin);
    }
    TFloat sysArtFlowDot = 0;
    if(sysArtValveAreaEff > 5e-3)
    {
        TFloat sysArtValveInert = 6.28 * bloodDensity_ / sqrt(sysArtValveAreaEff);
        TFloat sysArtValveBorda = bloodDensity_ / (2 * sysArtValveAreaEff*sysArtValveAreaEff);
        sysArtFlowDot = (lvPressure - sysArtCompliPressure - sysArtResist_*sysArtFlow - sysArtValveBorda*sysArtFlow*std::abs(sysArtFlow)) / sysArtValveInert;
    }
    
    //---------- RavValve ----------
    
    TFloat ravValveAreaEff;
    if(raPressure > rvPressure || ravFlow > 0)
    {
        TFloat ravValveAreaMax = ravValveMax_ * ravValveAreaRef_;
        ravValveAreaEff = ravValveAreaMax*ravValveAreaRef_ / (ravValveAreaRef_-ravValveAreaMax);
    }
    else
    {
        TFloat ravValveAreaMin = ravValveMin_ * ravValveAreaRef_;
        ravValveAreaEff = ravValveAreaMin*ravValveAreaRef_/(ravValveAreaRef_-ravValveAreaMin);
    }
    TFloat ravFlowDot = 0;
    if(ravValveAreaEff > 5e-3)
    {
        TFloat ravValveInert = 6.28 * bloodDensity_ / sqrt(ravValveAreaEff);
        TFloat ravValveBorda = bloodDensity_ / (2 * ravValveAreaEff*ravValveAreaEff);
        ravFlowDot = (raPressure - rvPressure - ravValveBorda*ravFlow*std::abs(ravFlow)) / ravValveInert;
    }
    
    //---------- PulArtValve ----------
    
    TFloat pulArtValveAreaEff;
    if(rvPressure > pulArtCompliPressure || pulArtFlow > 0)
    {
        TFloat pulArtValveAreaMax = pulArtValveMax_ * pulArtValveAreaRef_;
        pulArtValveAreaEff = pulArtValveAreaMax*pulArtValveAreaRef_ / (pulArtValveAreaRef_-pulArtValveAreaMax);
    }
    else
    {
        TFloat pulArtValveAreaMin = pulArtValveMin_ * pulArtValveAreaRef_;
        pulArtValveAreaEff = pulArtValveAreaMin*pulArtValveAreaRef_/(pulArtValveAreaRef_-pulArtValveAreaMin);
    }
    TFloat pulArtFlowDot = 0;
    if(pulArtValveAreaEff > 5e-3)
    {
        TFloat pulArtValveInert = 6.28 * bloodDensity_ / sqrt(pulArtValveAreaEff);
        TFloat pulArtValveBorda = bloodDensity_ / (2 * pulArtValveAreaEff*pulArtValveAreaEff);
        pulArtFlowDot = (rvPressure - pulArtCompliPressure - pulArtResist_*pulArtFlow - pulArtValveBorda*pulArtFlow*std::abs(pulArtFlow)) / pulArtValveInert;
    }
    
    //---------- LavValve ----------
    
    TFloat lavValveAreaEff;
    if(laPressure > lvPressure || lavFlow > 0)
    {
        TFloat lavValveAreaMax = lavValveMax_ * lavValveAreaRef_;
        lavValveAreaEff = lavValveAreaMax*lavValveAreaRef_ / (lavValveAreaRef_-lavValveAreaMax);
    }
    else
    {
        TFloat lavValveAreaMin = lavValveMin_ * lavValveAreaRef_;
        lavValveAreaEff = lavValveAreaMin*lavValveAreaRef_/(lavValveAreaRef_-lavValveAreaMin);
    }
    TFloat lavFlowDot = 0;
    if(lavValveAreaEff > 5e-3)
    {
        TFloat lavValveInert = 6.28 * bloodDensity_ / sqrt(lavValveAreaEff);
        TFloat lavValveBorda = bloodDensity_ / (2 * lavValveAreaEff*lavValveAreaEff);
        lavFlowDot = (laPressure - lvPressure - lavValveBorda*lavFlow*std::abs(lavFlow)) / lavValveInert;
    }
    
    //----------
    
    stateVarsDot[ 0] =    lavFlow  - sysArtFlow;
    stateVarsDot[ 1] = sysArtFlow  - sysPerFlow_;
    stateVarsDot[ 2] = sysPerFlow_ - sysVenFlow_;
    stateVarsDot[ 3] = sysVenFlow_ -    ravFlow;
    stateVarsDot[ 4] =    ravFlow  - pulArtFlow;
    stateVarsDot[ 5] = pulArtFlow  - pulPerFlow_;
    stateVarsDot[ 6] = pulPerFlow_ - pulVenFlow_;
    stateVarsDot[ 7] = pulVenFlow_ -    lavFlow;
    stateVarsDot[ 8] = sysArtFlowDot;
    stateVarsDot[ 9] =    ravFlowDot;
    stateVarsDot[10] = pulArtFlowDot;
    stateVarsDot[11] =    lavFlowDot;
    stateVarsDot[12] = sysArtFlow - pulArtFlow;
}
