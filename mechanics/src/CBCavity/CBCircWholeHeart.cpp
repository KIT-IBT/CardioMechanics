/*
 * File: CBCircWholeHeart.cpp
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


#include "CBCircWholeHeart.h"

CBCircWholeHeart::CBCircWholeHeart(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "LvVolume";
    stateVarsNames_[1] = "SysArtVolume";
    stateVarsNames_[2] = "SysVenVolume";
    stateVarsNames_[3] = "RaVolume";
    stateVarsNames_[4] = "RvVolume";
    stateVarsNames_[5] = "PulArtVolume";
    stateVarsNames_[6] = "PulVenVolume";
    stateVarsNames_[7] = "LaVolume";
    stateVarsNames_[8] = "StrokeVolumeDiff";
    
    //----- Results -----
    
    results_[ 0] = &lvPressure_;        resultsNames_[ 0] = "LvPressure";
    results_[ 1] = &sysArtPressure_;    resultsNames_[ 1] = "SysArtPressure";
    results_[ 2] = &sysVenPressure_;    resultsNames_[ 2] = "SysVenPressure";
    results_[ 3] = &raPressure_;        resultsNames_[ 3] = "RaPressure";
    results_[ 4] = &rvPressure_;        resultsNames_[ 4] = "RvPressure";
    results_[ 5] = &pulArtPressure_;    resultsNames_[ 5] = "PulArtPressure";
    results_[ 6] = &pulVenPressure_;    resultsNames_[ 6] = "PulVenPressure";
    results_[ 7] = &laPressure_;        resultsNames_[ 7] = "LaPressure";
    results_[ 8] = &sysArtFlow_;        resultsNames_[ 8] = "SysArtFlow";
    results_[ 9] = &sysPerFlow_;        resultsNames_[ 9] = "SysPerFlow";
    results_[10] = &sysVenFlow_;        resultsNames_[10] = "SysVenFlow";
    results_[11] = &ravFlow_;           resultsNames_[11] = "RavFlow";
    results_[12] = &pulArtFlow_;        resultsNames_[12] = "PulArtFlow";
    results_[13] = &pulPerFlow_;        resultsNames_[13] = "PulPerFlow";
    results_[14] = &pulVenFlow_;        resultsNames_[14] = "PulVenFlow";
    results_[15] = &lavFlow_;           resultsNames_[15] = "LavFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &lvPressure_;
    estCavityPressures_[1] = &raPressure_;
    estCavityPressures_[2] = &rvPressure_;
    estCavityPressures_[3] = &laPressure_;
    
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
    ravValveResist_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.RavValveResist");
    pulArtValveResist_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtValveResist");
    pulArtResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtResist");
    pulArtCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtCompli");
    pulArtVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulArtVolumeUnstr");
    pulPerResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulPerResist");
    pulVenResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenResist");
    pulVenCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenCompli");
    pulVenVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenVolumeUnstr");
    lavValveResist_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveResist");
    
    totalVolume_  = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtVolume");
    stateVars_[5] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulArtVolume");
    stateVars_[6] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.PulVenVolume");
    stateVars_[8] = 0;
    
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
    
    InitSteadyStateCheck(parameters, parameterPrefix, 8);
    
    WriteHeaderToFile();
}


void CBCircWholeHeart::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[3] = cavityVolumes[1];
    stateVars_[4] = cavityVolumes[2];
    stateVars_[7] = cavityVolumes[3];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[1]-stateVars_[3]-stateVars_[4]-stateVars_[5]-stateVars_[6]-stateVars_[7];
}


void CBCircWholeHeart::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
    cavityVolumes[1] = estStateVars_[3];
    cavityVolumes[2] = estStateVars_[4];
    cavityVolumes[3] = estStateVars_[7];
}


void CBCircWholeHeart::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat lvPressure = cavityPressures[0];
    TFloat raPressure = cavityPressures[1];
    TFloat rvPressure = cavityPressures[2];
    TFloat laPressure = cavityPressures[3];
    
    TFloat sysArtVolume = stateVars[1];
    TFloat sysVenVolume = stateVars[2];
    TFloat pulArtVolume = stateVars[5];
    TFloat pulVenVolume = stateVars[6];
    
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
    sysVenFlow_ = (sysVenPressure_-raPressure)/sysVenResist_;
    
    if(raPressure >= rvPressure)
        ravFlow_ = (raPressure-rvPressure)/ravValveResist_;
    else
        ravFlow_ = 0.0;
    
    if(rvPressure >= pulArtCompliPressure)
        pulArtFlow_ = (rvPressure-pulArtCompliPressure)/(pulArtValveResist_+pulArtResist_);
    else
        pulArtFlow_ = 0.0;
    pulArtPressure_ = pulArtCompliPressure + pulArtFlow_ * pulArtResist_;
    
    pulPerFlow_ = (pulArtCompliPressure-pulVenPressure_)/pulPerResist_;
    pulVenFlow_ = (pulVenPressure_-laPressure)/pulVenResist_;
    
    if(laPressure >= lvPressure)
        lavFlow_ = (laPressure-lvPressure)/lavValveResist_;
    else
        lavFlow_ = 0.0;
}


void CBCircWholeHeart::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] =    lavFlow_ - sysArtFlow_;
    stateVarsDot[1] = sysArtFlow_ - sysPerFlow_;
    stateVarsDot[2] = sysPerFlow_ - sysVenFlow_;
    stateVarsDot[3] = sysVenFlow_ -    ravFlow_;
    stateVarsDot[4] =    ravFlow_ - pulArtFlow_;
    stateVarsDot[5] = pulArtFlow_ - pulPerFlow_;
    stateVarsDot[6] = pulPerFlow_ - pulVenFlow_;
    stateVarsDot[7] = pulVenFlow_ -    lavFlow_;
    stateVarsDot[8] = sysArtFlow_ - pulArtFlow_;
}
