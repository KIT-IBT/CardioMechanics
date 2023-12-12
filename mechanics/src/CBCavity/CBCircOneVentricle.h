/*
 * File: CBCircOneVentricle.h
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
