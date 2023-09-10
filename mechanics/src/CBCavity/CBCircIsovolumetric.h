/*
 *  CBCircIsovolumetric.h
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 29.01.20.
 *  Copyright 2020 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
