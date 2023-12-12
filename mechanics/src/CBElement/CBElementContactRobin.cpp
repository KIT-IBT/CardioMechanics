/*
 * File: CBElementContactRobin.cpp
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


#include "CBElementContactRobin.h"

CBElementContactRobin::CBElementContactRobin(CBElementContactRobin &other) : CBElementSurfaceT3(other) {
    referenceNormal_ = other.referenceNormal_;
}

CBElement *CBElementContactRobin::New() {
    return new CBElementContactRobin;
}

CBElement *CBElementContactRobin::Clone() {
    return new CBElementContactRobin(*this);
}

void CBElementContactRobin::GetNodesCoords(TFloat *nodesCoords) {
    TInt nodesCoordsIndices[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * Ancestor::nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * Ancestor::nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * Ancestor::nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
}
