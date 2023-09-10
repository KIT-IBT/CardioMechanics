/*
 *  CBCircWholeHeartValvesDynamic.cpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 10.03.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCircWholeHeartValvesDynamic.h"

CBCircWholeHeartValvesDynamic::CBCircWholeHeartValvesDynamic(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
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
    stateVarsNames_[12] = "SysArtValveState";
    stateVarsNames_[13] = "RavValveState";
    stateVarsNames_[14] = "PulArtValveState";
    stateVarsNames_[15] = "LavValveState";
    stateVarsNames_[16] = "StrokeVolumeDiff";
    
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
    sysArtValveRateOpening_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveRateOpening");
    sysArtValveRateClosing_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveRateClosing");
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
    ravValveRateOpening_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveRateOpening");
    ravValveRateClosing_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveRateClosing");
    pulArtValveMax_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveMax");
    pulArtValveMin_         = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveMin");
    pulArtValveAreaRef_     = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveAreaRef");
    pulArtValveRateOpening_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveRateOpening");
    pulArtValveRateClosing_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveRateClosing");
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
    lavValveRateOpening_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveRateOpening");
    lavValveRateClosing_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveRateClosing");
    
    totalVolume_   = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[ 1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtVolume");
    stateVars_[ 5] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtVolume");
    stateVars_[ 6] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulVenVolume");
    stateVars_[ 8] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtFlow", 0.0);
    stateVars_[ 9] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RavFlow", 0.0);
    stateVars_[10] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtFlow", 0.0);
    stateVars_[11] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LavFlow", 0.0);
    stateVars_[12] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtValveState", 0.0);
    stateVars_[13] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.RavValveState", 0.0);
    stateVars_[14] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtValveState", 0.0);
    stateVars_[15] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LavValveState", 0.0);
    
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
    
    InitSteadyStateCheck(parameters, parameterPrefix, 16);
    
    WriteHeaderToFile();
}


void CBCircWholeHeartValvesDynamic::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[3] = cavityVolumes[1];
    stateVars_[4] = cavityVolumes[2];
    stateVars_[7] = cavityVolumes[3];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[1]-stateVars_[3]-stateVars_[4]-stateVars_[5]-stateVars_[6]-stateVars_[7];
}


void CBCircWholeHeartValvesDynamic::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
    cavityVolumes[1] = estStateVars_[3];
    cavityVolumes[2] = estStateVars_[4];
    cavityVolumes[3] = estStateVars_[7];
}


void CBCircWholeHeartValvesDynamic::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
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


void CBCircWholeHeartValvesDynamic::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
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
    TFloat sysArtValveState = stateVars[12];
    TFloat ravValveState    = stateVars[13];
    TFloat pulArtValveState = stateVars[14];
    TFloat lavValveState    = stateVars[15];
    
    TFloat sysArtCompliPressure = (sysArtVolume-sysArtVolumeUnstr_)/sysArtCompli_;
    TFloat pulArtCompliPressure = (pulArtVolume-pulArtVolumeUnstr_)/pulArtCompli_;
    
    //---------- SysArtValve ----------
    
    TFloat sysArtValveIndex = (sysArtValveMax_-sysArtValveMin_)*sysArtValveState + sysArtValveMin_;
    
    TFloat sysArtValveAreaEff = sysArtValveAreaRef_ * sysArtValveIndex / (1.0-sysArtValveIndex);
    
    TFloat sysArtFlowDot = 0;
    if(sysArtValveAreaEff > 5e-3) // make sure not to divide by zero
    {
        TFloat sysArtValveInert = 6.28 * bloodDensity_ / sqrt(sysArtValveAreaEff);
        TFloat sysArtValveBorda = bloodDensity_ / (2 * sysArtValveAreaEff*sysArtValveAreaEff);
        sysArtFlowDot = (lvPressure - sysArtCompliPressure - sysArtResist_*sysArtFlow - sysArtValveBorda*sysArtFlow*std::abs(sysArtFlow)) / sysArtValveInert;
    }
    
    TFloat sysArtValveStateDot;
    if(lvPressure > sysArtCompliPressure)
        sysArtValveStateDot = (1-sysArtValveState) * sysArtValveRateOpening_ * (lvPressure-sysArtCompliPressure);
    else
        sysArtValveStateDot = sysArtValveState * sysArtValveRateClosing_ * (lvPressure-sysArtCompliPressure);
    
    //---------- RavValve ----------
    
    TFloat ravValveIndex = (ravValveMax_-ravValveMin_)*ravValveState + ravValveMin_;
    
    TFloat ravValveAreaEff = ravValveAreaRef_ * ravValveIndex / (1.0-ravValveIndex);
    
    TFloat ravFlowDot = 0;
    if(ravValveAreaEff > 5e-3)
    {
        TFloat ravValveInert = 6.28 * bloodDensity_ / sqrt(ravValveAreaEff);
        TFloat ravValveBorda = bloodDensity_ / (2 * ravValveAreaEff*ravValveAreaEff);
        ravFlowDot = (raPressure - rvPressure - ravValveBorda*ravFlow*std::abs(ravFlow) - 0.000*ravFlow) / ravValveInert;
    }
    
    TFloat ravValveStateDot;
    if(raPressure > rvPressure)
        ravValveStateDot = (1-ravValveState) * ravValveRateOpening_ * (raPressure-rvPressure);
    else
        ravValveStateDot = ravValveState * ravValveRateClosing_ * (raPressure-rvPressure);
    
    //---------- PulArtValve ----------
    
    TFloat pulArtValveIndex = (pulArtValveMax_-pulArtValveMin_)*pulArtValveState + pulArtValveMin_;
    
    TFloat pulArtValveAreaEff = pulArtValveAreaRef_ * pulArtValveIndex / (1.0-pulArtValveIndex);
    
    TFloat pulArtFlowDot = 0;
    if(pulArtValveAreaEff > 5e-3)
    {
        TFloat pulArtValveInert = 6.28 * bloodDensity_ / sqrt(pulArtValveAreaEff);
        TFloat pulArtValveBorda = bloodDensity_ / (2 * pulArtValveAreaEff*pulArtValveAreaEff);
        pulArtFlowDot = (rvPressure - pulArtCompliPressure - pulArtResist_*pulArtFlow - pulArtValveBorda*pulArtFlow*std::abs(pulArtFlow)) / pulArtValveInert;
    }
    
    TFloat pulArtValveStateDot;
    if(rvPressure > pulArtCompliPressure)
        pulArtValveStateDot = (1-pulArtValveState) * pulArtValveRateOpening_ * (rvPressure-pulArtCompliPressure);
    else
        pulArtValveStateDot = pulArtValveState * pulArtValveRateClosing_ * (rvPressure-pulArtCompliPressure);
    
    //---------- LavValve ----------
    
    TFloat lavValveIndex = (lavValveMax_-lavValveMin_)*lavValveState + lavValveMin_;
    
    TFloat lavValveAreaEff = lavValveAreaRef_ * lavValveIndex / (1.0-lavValveIndex);
    
    TFloat lavFlowDot = 0;
    if(lavValveAreaEff > 5e-3)
    {
        TFloat lavValveInert = 6.28 * bloodDensity_ / sqrt(lavValveAreaEff);
        TFloat lavValveBorda = bloodDensity_ / (2 * lavValveAreaEff*lavValveAreaEff);
        lavFlowDot = (laPressure - lvPressure - lavValveBorda*lavFlow*std::abs(lavFlow)) / lavValveInert;
    }
    
    TFloat lavValveStateDot;
    if(laPressure > lvPressure)
        lavValveStateDot = (1-lavValveState) * lavValveRateOpening_ * (laPressure-lvPressure);
    else
        lavValveStateDot = lavValveState * lavValveRateClosing_ * (laPressure-lvPressure);
    
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
    stateVarsDot[12] = sysArtValveStateDot;
    stateVarsDot[13] =    ravValveStateDot;
    stateVarsDot[14] = pulArtValveStateDot;
    stateVarsDot[15] =    lavValveStateDot;
    stateVarsDot[16] = sysArtFlow - pulArtFlow;
}
