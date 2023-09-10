/*
 *  CBFormulation.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#ifndef CB_FORMULATION_H
#define CB_FORMULATION_H

#include "Matrix3.h"


#include "CBStatus.h"


using namespace math_pack;

class CBSolver;


class CBFormulation
{
public:
    CBFormulation(CBSolver* solver){solver_ = solver; }
    virtual ~CBFormulation(){}
    virtual CBStatus  CalcNodalForces() = 0;
    virtual CBStatus  CalcNodalForcesJacobian() = 0;
protected:
    CBSolver* solver_;    // No ownership
private:
    CBFormulation();
};

#endif

