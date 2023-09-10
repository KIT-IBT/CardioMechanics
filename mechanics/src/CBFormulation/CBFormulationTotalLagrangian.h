/*
 *  CBFormulationTotalLagrangian.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#ifndef CB_FORMULATION_TOTAL_LAGRANGIAN_H
#define CB_FORMULATION_TOTAL_LAGRANGIAN_H

#include "CBFormulation.h"



class CBFormulationTotalLagrangian : public CBFormulation
{
public:
    CBFormulationTotalLagrangian(CBSolver* solver) : CBFormulation(solver){}
    virtual ~CBFormulationTotalLagrangian(){}
    CBStatus CalcNodalForces();
    CBStatus CalcNodalForcesJacobian();
protected:
private:
    typedef CBFormulation   Base;
};

#endif
