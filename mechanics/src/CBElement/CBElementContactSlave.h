/*
 *  CBContactHandlingSlaveElement.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 20.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#ifndef CB_CONTACT_SLAVE_ELEMENT_H
#define CB_CONTACT_SLAVE_ELEMENT_H

#include "CBElementSurfaceT3.h"
#include "CBElementAdapter.h"

using namespace math_pack;

class CBElementContactMaster;


class CBElementContactSlave : public CBElementSurfaceT3
{
public:
    CBElementContactSlave( CBElementContactSlave& other);
    CBElementContactSlave(){}
    virtual ~CBElementContactSlave(){}
    CBElement* Clone() override;
    
    static CBElement* New();
    std::string GetType() override {return(std::string("CONTACT_SLAVE")); }
    
    void GetNodesCoords(TFloat* nodesCoords);
    
protected:
private:
    friend class CBElementContactMaster;
    
    void operator=(const CBElementContactSlave&);
    
    typedef CBElementSurfaceT3   Ancestor;
    typedef CBElement            Base;
};

#endif
