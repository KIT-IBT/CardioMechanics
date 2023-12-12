/*
 * File: CBCircWindkessel3AsymptoticPressure.h
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


#ifndef CB_CIRC_WINDKESSEL_3_ASYMPTOTIC_PRESSURE_H
#define CB_CIRC_WINDKESSEL_3_ASYMPTOTIC_PRESSURE_H

#include "CBCircModel.h"

class CBCircWindkessel3AsymptoticPressure : public CBCircModel<1,2,5>
{
public:
    CBCircWindkessel3AsymptoticPressure(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircWindkessel3AsymptoticPressure() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat artResist_, artCompli_;
    TFloat perResist_, asymptoticPressure_;
    TFloat venResist_, venPressure_;
    
    //----- Results -----
    
    TFloat ventrPressure_, artPressure_;
    TFloat artFlow_, perFlow_, venFlow_;
};

#endif
