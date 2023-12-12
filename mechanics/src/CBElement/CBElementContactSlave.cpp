/*
 * File: CBElementContactSlave.cpp
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



#include <cmath>

#include "CBElementContactSlave.h"


CBElementContactSlave::CBElementContactSlave(CBElementContactSlave& other) :CBElementSurfaceT3(other) {
}


CBElement* CBElementContactSlave::New() {
    return(new CBElementContactSlave);
}

CBElement* CBElementContactSlave::Clone() {
    return(new CBElementContactSlave(*this));
}


void CBElementContactSlave::GetNodesCoords(TFloat* nodesCoords)
{
    TInt   nodesCoordsIndices[9];
    
    for(unsigned int i = 0; i < 3; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * Ancestor::nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * Ancestor::nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * Ancestor::nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
}

