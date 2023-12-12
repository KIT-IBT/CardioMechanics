/*
 * File: CBFormulationTotalLagrangian.h
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
