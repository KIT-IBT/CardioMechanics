/*
 *  CBSolverEquilibrium
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_SOLVER_PETSC_EQUILIBRIUM
#define CB_SOLVER_PETSC_EQUILIBRIUM

#include "CBSolver.h"
#include <cstdio>


class CBSolverEquilibrium : public CBSolver
{
public:
    
    CBSolverEquilibrium() : CBSolver(){}
    virtual void Init(ParameterMap* parameters, CBModel* model);
    void DeInit();
    virtual std::string GetType(){return("Static Solver"); }
    friend PetscErrorCode CBSolverEquilibriumSNESHelperFunctionForces(SNES snes, Vec x, Vec f, void* solver);
    friend PetscErrorCode CBSolverEquilibriumSNESHelperFunctionForcesJacobian(SNES snes, Vec x, Mat jacobian, Mat preconditionerMatrix, void* solver);
    
protected:
    void InitVectors();
    void InitMatrices();
    void InitPETScSolver();
    void UpdateInitialGuess(TFloat time);
    CBStatus SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false);
    
    CBStatus CalcNodalForces(Vec displacement, Vec forces);
    CBStatus CalcNodalForcesJacobian(Vec displacement, Mat jacobian);
    typedef CBSolver   Base;
    
    SNES snes_;
    int snesStep_ = 0;
    KSP  ksp_;
    PC   pc_;
    Vec  residuum_;
    Vec  displacement_;
    Vec  initialGuess_;
    Vec  tmpVector_;
    
    // time tracking variables, for initial guess computation
    TFloat timeLast_;
    TFloat stepLast_;
private:
};

#endif
