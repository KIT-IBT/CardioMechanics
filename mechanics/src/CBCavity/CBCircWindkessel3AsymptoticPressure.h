/*
 *  CBCircWindkessel3AsymptoticPressure.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 26.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
