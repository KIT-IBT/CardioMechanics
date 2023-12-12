/*
 * File: CBFormulation.h
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

