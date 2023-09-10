/*
 *  CBCircWindkessel3.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 26.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRC_WINDKESSEL_3_H
#define CB_CIRC_WINDKESSEL_3_H

#include "CBCircModel.h"

class CBCircWindkessel3 : public CBCircModel<1,2,5>
{
public:
    CBCircWindkessel3(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircWindkessel3() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat artResist_, artCompli_;
    TFloat perResist_;
    TFloat venResist_, venPressure_;
    
    TFloat initialArtVolume_;
    bool resetArtVolume_;
    
    //----- Results -----
    
    TFloat ventrPressure_, artPressure_;
    TFloat artFlow_, perFlow_, venFlow_;
};

#endif
