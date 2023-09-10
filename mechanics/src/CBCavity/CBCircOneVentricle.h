/*
 *  CBCircOneVentricle.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 26.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRC_ONE_VENTRICLE_H
#define CB_CIRC_ONE_VENTRICLE_H

#include "CBCircModel.h"

class CBCircOneVentricle : public CBCircModel<1,3,6>
{
public:
    CBCircOneVentricle(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircOneVentricle() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat totalVolume_;
    
    TFloat artValveResist_;
    TFloat artResist_, artCompli_, artVolumeUnstr_;
    TFloat perResist_;
    TFloat venResist_, venCompli_, venVolumeUnstr_;
    
    //----- Results -----
    
    TFloat ventrPressure_, artPressure_, venPressure_;
    TFloat artFlow_, perFlow_, venFlow_;
};

#endif
