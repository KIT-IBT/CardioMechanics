/*
 * File: CBCircIsovolumetric.h
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


#ifndef CB_CIRC_ISOVOLUMETRIC_H
#define CB_CIRC_ISOVOLUMETRIC_H

#include "CBCircModel.h"

class CBCircIsovolumetric : public CBCircModel<1,1,2>
{
public:
    CBCircIsovolumetric(ParameterMap* parameters, std::string parameterPrefix);
    ~CBCircIsovolumetric() { };
    
    void SetInitialCavityVolumes(TFloat cavityVolumes[]);
    void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]);
    
private:
    void Algebraics(TFloat cavityPressures[], TFloat stateVars[]);
    void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]);
    
    TFloat initialCavityVolume_;
    
    //----- Results -----
    
    TFloat pressure_, volume_;
};

#endif
