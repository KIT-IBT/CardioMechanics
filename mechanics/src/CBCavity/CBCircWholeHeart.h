/*
 *  CBCircWholeHeart.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 16.01.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRC_WHOLE_HEART_H
#define CB_CIRC_WHOLE_HEART_H

#include "CBCircModel.h"

class CBCircWholeHeart : public CBCircModel<4,9,16>
{
public:
    CBCircWholeHeart(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircWholeHeart() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat totalVolume_;
    
    TFloat sysArtValveResist_;
    TFloat sysArtResist_, sysArtCompli_, sysArtVolumeUnstr_;
    TFloat sysPerResist_;
    TFloat sysVenResist_, sysVenCompli_, sysVenVolumeUnstr_;
    TFloat ravValveResist_;
    TFloat pulArtValveResist_;
    TFloat pulArtResist_, pulArtCompli_, pulArtVolumeUnstr_;
    TFloat pulPerResist_;
    TFloat pulVenResist_, pulVenCompli_, pulVenVolumeUnstr_;
    TFloat lavValveResist_;
    
    //----- Results -----
    
    TFloat lvPressure_, sysArtPressure_, sysVenPressure_, raPressure_;
    TFloat rvPressure_, pulArtPressure_, pulVenPressure_, laPressure_;
    TFloat sysArtFlow_, sysPerFlow_, sysVenFlow_, ravFlow_;
    TFloat pulArtFlow_, pulPerFlow_, pulVenFlow_, lavFlow_;
};

#endif
