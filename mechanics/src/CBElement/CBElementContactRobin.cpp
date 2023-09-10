/* -------------------------------------------------------
 
 CBElementContactRobin.cpp
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach (06.09.2021)
 Last modified: Tobias Gerach (06.09.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

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
