/*
 * File: CBCircLeftHeart.h
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


#ifndef CB_CIRC_LEFT_HEART_H
#define CB_CIRC_LEFT_HEART_H

#include "CBCircModel.h"

class CBCircLeftHeart : public CBCircModel<2,4,8>
{
public:
    CBCircLeftHeart(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircLeftHeart() { };
    
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
    TFloat pulVenResist_, pulVenCompli_, pulVenVolumeUnstr_;
    TFloat lavValveResist_;
    
    //----- Results -----
    
    TFloat lvPressure_, sysArtPressure_;
    TFloat pulVenPressure_, laPressure_;
    TFloat sysArtFlow_, sysPerFlow_;
    TFloat pulVenFlow_, lavFlow_;
};

#endif
