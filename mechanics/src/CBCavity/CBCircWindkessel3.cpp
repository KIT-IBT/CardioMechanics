/*
 * File: CBCircWindkessel3.cpp
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


#include "CBCircWindkessel3.h"

CBCircWindkessel3::CBCircWindkessel3(ParameterMap* parameters, std::string parameterPrefix) : CBCircModel(parameters, parameterPrefix)
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
    
    artResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtResist");
    artCompli_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.ArtCompli");
    perResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.PerResist");
    venResist_      = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenResist");
    venPressure_    = parameters->Get<TFloat>(parameterPrefix + ".CircParameters.VenPressure");
    
    initialArtVolume_ = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.ArtVolume");
    resetArtVolume_   = parameters->Get<bool>(parameterPrefix + ".InitialConditions.ResetArtVolume", false);
    stateVars_[1]   = initialArtVolume_;
    
    
    cavitySurfaceIndices_[0] = parameters->Get<TInt>(parameterPrefix + ".Cavities.VentrSurface");
    
    cavityPreloadingPressures_[0] = parameters->Get<TFloat>(parameterPrefix + ".Cavities.VentrPreloading", 0.0);
    cavityPreloadingPressures_[1] = parameters->Get<TFloat>(parameterPrefix + ".InitialConditions.VentrPressure", 0.0);
    
    InitSteadyStateCheck(parameters, parameterPrefix, -1);
    
    WriteHeaderToFile();
}


void CBCircWindkessel3::SetInitialCavityVolumes(TFloat cavityVolumes[])
{
    stateVars_[0] = cavityVolumes[0];
}


void CBCircWindkessel3::GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[])
{
    Integrate(cavityPressures, timeStep);
    
    cavityVolumes[0] = estStateVars_[0];
}


void CBCircWindkessel3::Algebraics(TFloat cavityPressures[], TFloat stateVars[])
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
        
        // reset arterial pressure to be equal to cavity pressure, avoids opening of the valve by hard resetting
        // use case: atrial afterload with large perResist
        // this prevents increasing arterial volume over time when simulating multiple heart beats,
        // but should not be active when used as ventricular afterload (there the capacitor unloads over perResist)
        if (resetArtVolume_)
        {
            estStateVars_[1] = std::max(initialArtVolume_, artCompli_*ventrPressure);
            artPressure_ = estStateVars_[1]/artCompli_;
        }
    }
    
    perFlow_ = artCompliPressure/perResist_;
    
    if(venPressure_ >= ventrPressure)
        venFlow_ = (venPressure_-ventrPressure)/venResist_;
    else
        venFlow_ = 0.0;
}


void CBCircWindkessel3::Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[])
{
    Algebraics(cavityPressures, stateVars);
    
    stateVarsDot[0] = venFlow_ - artFlow_;
    stateVarsDot[1] = artFlow_ - perFlow_;
}
