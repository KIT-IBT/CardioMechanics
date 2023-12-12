/*
 * File: CBElementContactMaster.cpp
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


#include "CBElementContactMaster.h"
#include "CBElementContactSlave.h"

CBElementContactMaster::CBElementContactMaster(CBElementContactMaster& other) : CBElementSurfaceT3(other) {
    distanceVectorToSlave_ = other.distanceVectorToSlave_;
    distanceToSlave_ = other.distanceToSlave_;
    slaveAtGaussPoint_ = other.slaveAtGaussPoint_;
    slaveAtVertex_ = other.slaveAtVertex_;
}

CBElement* CBElementContactMaster::New() {
    return(new CBElementContactMaster);
}

CBElement* CBElementContactMaster::Clone() {
    return(new CBElementContactMaster(*this));
}

void CBElementContactMaster::GetNodesCoords(TFloat* nodesCoords) {
    TInt nodesCoordsIndices[9];
    
    for(unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * Ancestor::nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * Ancestor::nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * Ancestor::nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
}
