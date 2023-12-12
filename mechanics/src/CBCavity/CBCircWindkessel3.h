/*
 * File: CBCircWindkessel3.h
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
