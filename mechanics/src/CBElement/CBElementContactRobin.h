/* -------------------------------------------------------
 
 CBElementContactRobin.h
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach (06.09.2021)
 Last modified: Tobias Gerach (08.11.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#ifndef CB_ELEMENT_CONTACT_ROBIN_H
#define CB_ELEMENT_CONTACT_ROBIN_H

#include "CBElementSurfaceT3.h"
#include "CBElementAdapter.h"
#include <array>

class CBElementContactRobin : public CBElementSurfaceT3 {
public:
    CBElementContactRobin(CBElementContactRobin &other);
    CBElementContactRobin() {}
    
    virtual ~CBElementContactRobin() {}
    
    CBElement *Clone() override;
    static CBElement *New();
    
    std::string GetType() override {return std::string("CONTACT_ROBIN");}
    
    void GetNodesCoords(TFloat *nodesCoords);
    
    void SetReferenceNormal(Vector3<TFloat> referenceNormal) {
        referenceNormal_ = referenceNormal;
    }
    
    Vector3<TFloat> GetReferenceNormal() {return referenceNormal_; }
    
protected:
private:
    CBElementContactRobin(const CBElementContactRobin &);
    
    void operator=(const CBElementContactRobin &);
    
    typedef CBElementSurfaceT3   Ancestor;
    typedef CBElement            Base;
    
    Vector3<TFloat> referenceNormal_ = {0, 0, 0};
};
#endif  // ifndef CB_ELEMENT_CONTACT_ROBIN_H
