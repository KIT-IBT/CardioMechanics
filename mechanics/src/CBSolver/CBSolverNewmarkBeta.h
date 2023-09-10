/*
 *  CBSolverNewmarkBeta
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_SOLVER_PETSC_NEWMARK_BETA
#define CB_SOLVER_PETSC_NEWMARK_BETA

#include "CBSolver.h"
#include "CBStatus.h"
#include <cstdio>


class CBSolverNewmarkBeta : public CBSolver {
public:
    CBSolverNewmarkBeta() : CBSolver(), isInitDampingParametersDone_(false), isInitMassMatrixDone_(false) {}
    
    ~CBSolverNewmarkBeta() {}
    
    TFloat GetKineticEnergy() override {return kineticEnergy_;}
    
    TFloat GetDampingEnergyDissipation() override {return dampingEnergyDissipation_;}
    
    void Init(ParameterMap *_parameter, CBModel *_model) override;
    void DeInit() override;
    
    std::string GetType() override {return "Newmark Beta Solver";  }
    
    void SetZeroVelocityAndAcceleration() override;
    void SetZeroDisplacement() override;
    friend PetscErrorCode CBSolverNewmarkBetaSNESHelperFunctionForces(SNES snes, Vec x, Vec f, void *solver);
    friend PetscErrorCode CBSolverNewmarkBetaSNESHelperFunctionForcesJacobian(SNES snes, Vec x, Mat jacobian,
                                                                              Mat preconditionerMatrix, void *_solver);
    void SetVelocity(std::vector<Vector3<TFloat>> vel) override;
    void SetAcceleration(std::vector<Vector3<TFloat>> acc) override;
    void ExportSNESMatrix(TFloat time) override;
    
protected:
    CBStatus SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false) override;
    CBStatus CalcNodalForcesJacobian(Vec displacement, Mat jacobian) override;
    CBStatus CalcNodalForcesJacobianAndDamping(Vec displacement, Mat jacobian);
    CBStatus CalcDampingMatrix();
    
    bool        useConsistentMassMatrix_;
    Vec         velocity_;
    Vec         acceleration_;
    Vec         displacement_;
    Vec         absDisplacement_;
    Vec         residuum_;
    Vec         tmpDisplacement_;
    Vec         tmpVelocity_;
    Vec         tmpVector_;
    Vec         initialGuess_;
    PetscScalar kineticEnergy_;
    PetscScalar dampingEnergyDissipation_ = 0;
    PetscScalar beta_;
    PetscScalar gamma_;
    PetscScalar globalRayleighAlpha_;
    PetscScalar globalRayleighBeta_;
    
private:
    CBStatus CalcNodalForces(Vec displacement, Vec forces);
    void InitVectors() override;
    void InitMatrices() override;
    void InitMassMatrix();
    void InitMassMatrixLumped();
    void InitMassMatrixConsistent();
    void InitParameters() override;
    void InitNewmarkBetaParameter();
    void InitDampingMatrix();
    void InitDampingParameter();
    void InitPETScSolver() override;
    void Export(PetscScalar time) override;
    
    typedef CBSolver   Base;
    
    SNES        snes_;
    int         snesStep_ = 0;
    KSP         ksp_;
    PC          pc_;
    
    Mat         massMatrix_;
    Mat         dampingMatrix_;
    
    PetscScalar prevTime_ = INFINITY;
    
    bool        isInitDampingParametersDone_;
    bool        isInitMassMatrixDone_;
    bool        updateJacobian_ = true;
    std::string solverType_ = "mumps";
    
    PetscInt    snesIts_ = 0;
}; // class CBSolverNewmarkBeta

#endif // ifndef CB_SOLVER_PETSC_NEWMARK_BETA
