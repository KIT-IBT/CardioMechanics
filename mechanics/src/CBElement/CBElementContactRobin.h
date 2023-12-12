/*
 * File: CBElementContactRobin.h
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
