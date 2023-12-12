/*
 * File: CBCircLeftHeart.cpp
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


#include "CBCircLeftHeart.h"

CBCircLeftHeart::CBCircLeftHeart(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
{
    //----- State Variables -----
    
    stateVarsNames_[0] = "LvVolume";
    stateVarsNames_[1] = "SysArtVolume";
    stateVarsNames_[2] = "PulVenVolume";
    stateVarsNames_[3] = "LaVolume";
    
    //----- Results -----
    
    results_[0] = &lvPressure_;        resultsNames_[0] = "LvPressure";
    results_[1] = &sysArtPressure_;    resultsNames_[1] = "SysArtPressure";
    results_[2] = &pulVenPressure_;    resultsNames_[2] = "PulVenPressure";
    results_[3] = &laPressure_;        resultsNames_[3] = "LaPressure";
    results_[4] = &sysArtFlow_;        resultsNames_[4] = "SysArtFlow";
    results_[5] = &sysPerFlow_;        resultsNames_[5] = "SysPerFlow";
    results_[6] = &pulVenFlow_;        resultsNames_[6] = "PulVenFlow";
    results_[7] = &lavFlow_;           resultsNames_[7] = "LavFlow";
    
    //----- Cavity Pressures and Indices -----
    
    estCavityPressures_[0] = &lvPressure_;
    estCavityPressures_[1] = &laPressure_;
    
    //----- Parameters -----
    
    auto pos = parameterPrefix.find_last_of(".")+1;
    name_ = parameterPrefix.substr(pos, parameterPrefix.size()-pos);
    
    sysArtValveResist_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtValveResist");
    sysArtResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtResist");
    sysArtCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtCompli");
    sysArtVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysArtVolumeUnstr");
    sysPerResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.SysPerResist");
    pulVenResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenResist");
    pulVenCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenCompli");
    pulVenVolumeUnstr_ = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PulVenVolumeUnstr");
    lavValveResist_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.LavValveResist");
    
    totalVolume_  = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.TotalVolume");
    stateVars_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.SysArtVolume");
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.LvSurface");
    cavitySurfaceIndices_[1] = parameters->Get<TInt>(parameterPrefix + ".Cavities.LaSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.LvPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.LaPreloading", 0.0);
    cavityPreloadingPressures_[2] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LvPressure", 0.0);
    cavityPreloadingPressures_[3] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.LaPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, -1);
    
    WriteHeaderToFile();
}


void CBCircLeftHeart::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
    stateVars_[3] = cavityVolumes[1];
    stateVars_[2] = totalVolume_-stateVars_[0]-stateVars_[3]-stateVars_[1];
}


void CBCircLeftHeart::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
    cavityVolumes[1] = estStateVars_[3];
}


void CBCircLeftHeart::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
{
    TFloat lvPressure = cavityPressures[0];
    TFloat laPressure = cavityPressures[1];
    
    TFloat sysArtVolume = stateVars[1];
    TFloat pulVenVolume = stateVars[2];
    
    TFloat sysArtCompliPressure = (sysArtVolume-sysArtVolumeUnstr_)/sysArtCompli_;
    pulVenPressure_ = (pulVenVolume-pulVenVolumeUnstr_)/pulVenCompli_;
    
    if(lvPressure >= sysArtCompliPressure)
        sysArtFlow_ = (lvPressure-sysArtCompliPressure)/(sysArtValveResist_+sysArtResist_);
    else
        sysArtFlow_ = 0.0;
    sysArtPressure_ = sysArtCompliPressure + sysArtFlow_ * sysArtResist_;
    
    sysPerFlow_ = (sysArtCompliPressure-pulVenPressure_)/sysPerResist_;
    
    pulVenFlow_ = (pulVenPressure_-laPressure)/pulVenResist_;
    
    if(laPressure >= lvPressure)
        lavFlow_ = (laPressure-lvPressure)/lavValveResist_;
    else
        lavFlow_ = 0.0;
}


void CBCircLeftHeart::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] =    lavFlow_ - sysArtFlow_;
    stateVarsDot[1] = sysArtFlow_ - sysPerFlow_;
    stateVarsDot[2] = sysPerFlow_ - pulVenFlow_;
    stateVarsDot[3] = pulVenFlow_ -    lavFlow_;
}
