/*
 *  CBCircWholeHeartValves.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 10.03.16.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRC_WHOLE_HEART_VALVES_H
#define CB_CIRC_WHOLE_HEART_VALVES_H

#include "CBCircModel.h"

class CBCircWholeHeartValves : public CBCircModel<4,13,12>
{
public:
    CBCircWholeHeartValves(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircWholeHeartValves() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat totalVolume_;
    
    TFloat bloodDensity_;
    
    TFloat sysArtValveMax_, sysArtValveMin_, sysArtValveAreaRef_;
    TFloat sysArtResist_, sysArtCompli_, sysArtVolumeUnstr_;
    TFloat sysPerResist_;
    TFloat sysVenResist_, sysVenCompli_, sysVenVolumeUnstr_;
    TFloat ravValveMax_, ravValveMin_, ravValveAreaRef_;
    TFloat pulArtValveMax_, pulArtValveMin_, pulArtValveAreaRef_;
    TFloat pulArtResist_, pulArtCompli_, pulArtVolumeUnstr_;
    TFloat pulPerResist_;
    TFloat pulVenResist_, pulVenCompli_, pulVenVolumeUnstr_;
    TFloat lavValveMax_, lavValveMin_, lavValveAreaRef_;
    
    //----- Results -----
    
    TFloat lvPressure_, sysArtPressure_, sysVenPressure_, raPressure_;
    TFloat rvPressure_, pulArtPressure_, pulVenPressure_, laPressure_;
    TFloat sysPerFlow_, sysVenFlow_, pulPerFlow_, pulVenFlow_;
};

#endif
