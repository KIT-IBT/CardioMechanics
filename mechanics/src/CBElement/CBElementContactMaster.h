/*
 *  CBContactHandlingInterfaceElement.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 31.03.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CONTACT_MASTER_ELEMENT_H
#define CB_CONTACT_MASTER_ELEMENT_H

#include "CBElementSurfaceT3.h"
#include "CBElementAdapter.h"
#include <array>

class CBElementContactSlave;

using namespace math_pack;

class CBElementContactMaster : public CBElementSurfaceT3 {
public:
    CBElementContactMaster(CBElementContactMaster &other);
    CBElementContactMaster() {}
    
    virtual ~CBElementContactMaster() {}
    
    CBElement *Clone() override;
    static CBElement *New();
    
    std::string GetType() override {return std::string("CONTACT_MASTER");}
    
    void SetSlaveAtGaussPoint(TInt i, TInt p) {slaveAtGaussPoint_[i] = p;}
    
    TInt GetSlaveAtGaussPoint(TInt i) {return slaveAtGaussPoint_[i];}
    
    void SetSlaveAtVertex(TInt i, TInt p) {slaveAtVertex_[i] = p;}
    
    TInt GetSlaveAtVertex(TInt i) {return slaveAtVertex_[i];}
    
    void GetNodesCoords(TFloat *nodesCoords);
    
    void SetDistanceToSlave(TFloat d) {distanceToSlave_ = d;}
    
    void SetDistanceVectorToSlave(int i, Vector3<TFloat> v) {distanceVectorToSlave_[i] = v;}
    
    Vector3<TFloat> GetDistanceVectorToSlave(int i) {return distanceVectorToSlave_[i];}
    
    TFloat GetDistanceToSlave() {return distanceToSlave_;}
    
    TFloat GetPrevDistanceToSlave() {return prevDistanceToSlave_;}
    
    void UpdatePrevDistanceToSlave() {prevDistanceToSlave_ = distanceToSlave_;}
    
    void SetReferenceNormal(Vector3<TFloat> referenceNormal) {
        referenceNormal_ = referenceNormal;
    }
    
    Vector3<TFloat> GetReferenceNormal() {return referenceNormal_; }
    
protected:
private:
    CBElementContactMaster(const CBElementContactMaster &);
    
    void operator=(const CBElementContactMaster &);
    
    typedef CBElementSurfaceT3   Ancestor;
    typedef CBElement            Base;
    std::array<Vector3<TFloat>, 3> distanceVectorToSlave_ = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    TFloat distanceToSlave_ = 0;
    TFloat prevDistanceToSlave_ = 0;
    std::array<TInt, 3>   slaveAtGaussPoint_ = {{-1, -1, -1}};
    std::array<TInt, 3>   slaveAtVertex_ = {{-1, -1, -1}};
    
    Vector3<TFloat> referenceNormal_ = {0, 0, 0};
};
#endif // ifndef CB_CONTACT_MASTER_ELEMENT_H
