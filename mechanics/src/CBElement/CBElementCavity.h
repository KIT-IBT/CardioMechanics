/*
 *  CBElementCavity.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 31.03.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_CAVITY_H
#define CB_ELEMENT_CAVITY_H

#include "CBElementSurface.h"

#include "Triangle.h"


/// CBElementCavity adds support to calculate volumes and apply pressure
class CBElementCavity : public CBElementSurface
{
public:
    CBElementCavity() : CBElementSurface() { }
    CBElementCavity(CBElementCavity& other) : CBElementSurface(other) { }
    virtual ~CBElementCavity(){}
    
    TFloat CalcContributionToVolume();
    virtual TFloat CalcContributionToVolume(const TFloat* referenceCoords) = 0;
    virtual void CalcContributionToVolumeJacobian(TFloat* volumeJacobianEntries, TInt* volumeJacobianEntriesIndices) = 0;
    
    virtual void ApplyPressure(TFloat pressure) = 0;
    virtual void CalcPressureJacobian(TFloat pressure) = 0;
    
    virtual math_pack::Triangle<TFloat> GetTriangle();
    
    // private:
    // CBElementCavity(const CBElementCavity&);
    //   void operator=(const CBElementCavity&);
};

#endif

