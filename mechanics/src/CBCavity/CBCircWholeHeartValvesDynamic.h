/*
 * File: CBCircWholeHeartValvesDynamic.h
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


#ifndef CB_CIRC_WHOLE_HEART_VALVES_DYNAMIC_H
#define CB_CIRC_WHOLE_HEART_VALVES_DYNAMIC_H

#include "CBCircModel.h"

class CBCircWholeHeartValvesDynamic : public CBCircModel<4,17,12>
{
public:
    CBCircWholeHeartValvesDynamic(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircWholeHeartValvesDynamic() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    //----- Parameters -----
    
    TFloat totalVolume_;
    
    TFloat bloodDensity_;
    
    TFloat sysArtValveMax_, sysArtValveMin_, sysArtValveAreaRef_;
    TFloat sysArtValveRateOpening_, sysArtValveRateClosing_;
    TFloat sysArtResist_, sysArtCompli_, sysArtVolumeUnstr_;
    TFloat sysPerResist_;
    TFloat sysVenResist_, sysVenCompli_, sysVenVolumeUnstr_;
    TFloat ravValveMax_, ravValveMin_, ravValveAreaRef_;
    TFloat ravValveRateOpening_, ravValveRateClosing_;
    TFloat pulArtValveMax_, pulArtValveMin_, pulArtValveAreaRef_;
    TFloat pulArtValveRateOpening_, pulArtValveRateClosing_;
    TFloat pulArtResist_, pulArtCompli_, pulArtVolumeUnstr_;
    TFloat pulPerResist_;
    TFloat pulVenResist_, pulVenCompli_, pulVenVolumeUnstr_;
    TFloat lavValveMax_, lavValveMin_, lavValveAreaRef_;
    TFloat lavValveRateOpening_, lavValveRateClosing_;
    
    //----- Results -----
    
    TFloat lvPressure_, sysArtPressure_, sysVenPressure_, raPressure_;
    TFloat rvPressure_, pulArtPressure_, pulVenPressure_, laPressure_;
    TFloat sysPerFlow_, sysVenFlow_, pulPerFlow_, pulVenFlow_;
};

#endif
