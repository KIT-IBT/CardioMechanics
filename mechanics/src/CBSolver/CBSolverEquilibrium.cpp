/**@file CBSolverEquilibrium.cpp
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.kit.edu/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

/*
 *  CBSolverEquilibrium.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 18.08.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#include "CBSolverEquilibrium.h"
#include <stdlib.h>

PetscErrorCode CBSolverEquilibriumSNESHelperFunctionForces(SNES snes, Vec u, Vec f, void *_solver) {
    CBSolverEquilibrium *solver = reinterpret_cast<CBSolverEquilibrium *>(_solver);
    CBStatus              rc     = solver->CalcNodalForces(u, f);
    
    if (rc != CBStatus::SUCCESS)
        SNESSetFunctionDomainError(snes);
    return 0;
}

PetscErrorCode CBSolverEquilibriumSNESHelperFunctionForcesJacobian(SNES snes, Vec u, Mat jacobian,
                                                                   Mat preconditionerMatrix, void *_solver) {
    CBSolverEquilibrium *solver = reinterpret_cast<CBSolverEquilibrium *>(_solver);
    CBStatus              rc     = solver->CalcNodalForcesJacobian(u, jacobian);
    
    if (rc != CBStatus::SUCCESS)
        SNESSetFunctionDomainError(snes);
    return 0;
}

void CBSolverEquilibrium::Init(ParameterMap *parameters, CBModel *model) {
    Base::Init(parameters, model);
    Base::status_ = CBStatus::WAITING;
    stepLast_ = 0;
    timeLast_ = 0;
}

CBStatus CBSolverEquilibrium::CalcNodalForces(Vec displacement, Vec forces) {
    snesStep_++;
    DCCtrl::debug<< "SNES Step:" << snesStep_ << " - Residuum: ";
    
    
    Vec localDisplacedNodesSeq;
    Vec displacedNodes;
    VecAssemblyBegin(displacement);
    VecAssemblyEnd(displacement);
    
    VecDuplicate(Base::nodes_, &displacedNodes);
    
    VecAssemblyBegin(displacedNodes);
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    VecAssemblyEnd(displacedNodes);
    
    VecCopy(Base::nodes_, displacedNodes);
    VecAXPY(displacedNodes, 1, displacement);
    
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(displacedNodes, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(displacedNodes, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    VecGhostGetLocalForm(displacedNodes, &localDisplacedNodesSeq);
    
    
    VecZeroEntries(Base::nodalForces_);
    
    Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    Base::adapter_->LinkNodalForces(Base::nodalForces_);
    CBStatus rc = Base::formulation_->CalcNodalForces();
    
    if (rc != CBStatus::SUCCESS)
        DCCtrl::cverbose << CBStatusToStr(rc);
    
    for (auto p : plugins_)
        p->ApplyToNodalForces();
    
    VecAssemblyBegin(Base::nodalForces_);
    VecAssemblyEnd(Base::nodalForces_);
    
    VecCopy(Base::nodalForces_, forces);
    VecDestroy(&displacedNodes);
    VecDestroy(&localDisplacedNodesSeq);
    
    PetscScalar r;
    VecNorm(forces, NORM_2, &r);
    DCCtrl::debug<< r <<"\n";
    
    int localSuccess;
    int globalSuccess;
    
    if (rc != CBStatus::SUCCESS)
        localSuccess = false;
    else
        localSuccess = true;
    MPI_Allreduce(&localSuccess, &globalSuccess, 1, MPI_INT, MPI_LAND, Petsc::Comm());
    if (!globalSuccess)
        return CBStatus::FAILED;
    else
        return CBStatus::SUCCESS;
} // CBSolverEquilibrium::CalcNodalForces

CBStatus CBSolverEquilibrium::CalcNodalForcesJacobian(Vec displacement, Mat jacobian) {
    double t = MPI_Wtime();
    DCCtrl::debug<< "\n\nCalculating System Matrix:\n";
    
    Base::adapter_->SetFiniteDifferencesEpsilon(CalcFiniteDifferencesEpsilon(displacement));
    Vec localDisplacedNodesSeq;
    Vec displacedNodes;
    
    VecDuplicate(Base::nodes_, &displacedNodes);
    VecAssemblyBegin(displacedNodes);
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    VecAssemblyEnd(displacedNodes);
    
    VecCopy(Base::nodes_, displacedNodes);
    VecAXPY(displacedNodes, 1, displacement);
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(displacedNodes, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(displacedNodes, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    VecGhostGetLocalForm(displacedNodes, &localDisplacedNodesSeq);
    
    MatZeroEntries(jacobian);
    MatSetLocalToGlobalMapping(jacobian, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    Base::adapter_->LinkNodalForcesJacobian(jacobian);
    
    double t1 = MPI_Wtime();
    DCCtrl::debug << "\t Elements ...";
    
    CBStatus rc = Base::formulation_->CalcNodalForcesJacobian();
    
    if (rc != CBStatus::SUCCESS)
        DCCtrl::cverbose << CBStatusToStr(rc);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]\n";
    
    for (std::vector<CBSolverPlugin *>::iterator it = Base::plugins_.begin(); it != Base::plugins_.end(); it++) {
        t1 = MPI_Wtime();
        DCCtrl::debug << "\t Plugin: " <<  (*it)->GetName() << " ...";
        
        (*it)->ApplyToNodalForcesJacobian();
        
        DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]\n";
    }
    
    VecDestroy(&displacedNodes);
    VecDestroy(&localDisplacedNodesSeq);
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "\t Assembling ...";
    
    MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
    
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]\n";
    
    MatAXPY(jacobian, 1, boundaryConditionsNodalForcesJacobianDiagonalComponents_, DIFFERENT_NONZERO_PATTERN);
    
    DCCtrl::debug << "Total time [" << MPI_Wtime() - t << " s]\n\n";
    
    int localSuccess;
    int globalSuccess;
    
    if (rc != CBStatus::SUCCESS)
        localSuccess = false;
    else
        localSuccess = true;
    MPI_Allreduce(&localSuccess, &globalSuccess, 1, MPI_INT, MPI_LAND, Petsc::Comm());
    if (!globalSuccess)
        return CBStatus::FAILED;
    else
        return CBStatus::SUCCESS;
} // CBSolverEquilibrium::CalcNodalForcesJacobian

void CBSolverEquilibrium::InitPETScSolver() {
    SNESCreate(PETSC_COMM_WORLD, &snes_);
    SNESSetTolerances(snes_, precision_, precision_, precision_, maxSnesIts_, maxFunEval_);
    SNESSetFunction(snes_, residuum_, CBSolverEquilibriumSNESHelperFunctionForces, (void *)this);
    SNESSetJacobian(snes_, Base::nodalForcesJacobian_, Base::nodalForcesJacobian_,
                    CBSolverEquilibriumSNESHelperFunctionForcesJacobian, (void *)this);
    SNESGetKSP(snes_, &ksp_);
    KSPGetPC(ksp_, &pc_);
    
    if (Base::parameters_->Get<bool>("Solver.LU", true)) {
        PCSetType(pc_, PCLU);
        KSPSetType(ksp_, "preonly");
    }
    
    PCFactorSetMatSolverPackage(pc_, "mumps");
    SNESSetFromOptions(snes_);
    KSPSetFromOptions(ksp_);
    PCSetFromOptions(pc_);
}

void CBSolverEquilibrium::InitVectors() {
    Base::InitNodalForces();
    
    VecDuplicate(Base::nodalForces_, &residuum_);
    VecDuplicate(Base::nodes_, &displacement_);
    VecDuplicate(Base::nodes_, &initialGuess_);
    VecDuplicate(Base::nodes_, &tmpVector_);
    VecZeroEntries(residuum_);
    VecZeroEntries(displacement_);
    VecZeroEntries(initialGuess_);
    VecZeroEntries(tmpVector_);
}

void CBSolverEquilibrium::InitMatrices() {
    Base::InitNodalForcesJacobian();
}

void CBSolverEquilibrium::DeInit() {
    VecDestroy(&residuum_);
    VecDestroy(&displacement_);
    VecDestroy(&initialGuess_);
    VecDestroy(&tmpVector_);
    SNESDestroy(&snes_);
    KSPDestroy(&ksp_);
    PCDestroy(&pc_);
    Base::DeInit();
}

void CBSolverEquilibrium::UpdateInitialGuess(TFloat time) {
    // intelligently detect STEPBACK, REPEAT, SUCCESS from changing time
    // variables and adapt the initial guess accordingly
    TFloat step = timing_.GetTimeStep();
    
    if (time > timeLast_) {
        // SUCCESS, save displacement from last computation in a normalized form
        // initialGuess_ = displacement_/stepLast_
        VecCopy(displacement_, initialGuess_);
        VecScale(initialGuess_, 1./stepLast_);
        
        // displacement_ = step * initialGuess_
        VecScale(displacement_, 1.*step/stepLast_);
    } else if (time == timeLast_) {
        // REPEAT, solution might be near the final one:
        // displacement_ = 0
        VecZeroEntries(displacement_);
    } else if (time < timeLast_) {
        // STEPBACK, reconstruct displacement_ from initialGuess_ variable
        // displacement_ = step * initialGuess_
        VecCopy(initialGuess_, displacement_);
        VecScale(displacement_, step);
    } else {
        throw std::runtime_error("CBSolverStatic::SolverStep() Could not detect the reason for doing this step.");
    }
    
    // save time and step, used during the next run
    timeLast_ = time;
    stepLast_ = step;
} // CBSolverEquilibrium::UpdateInitialGuess

CBStatus CBSolverEquilibrium::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation) {
    // the solution is the displacement to previous time step node coords
    
    SNESConvergedReason snesReason;
    KSPConvergedReason  kspReason;
    PetscInt            snesIts    = 0;
    PetscInt            kspIts     = 0;
    
    UpdateInitialGuess(time);
    
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    
    UpdateGhostNodesAndLinkToAdapter();
    
    for (std::vector<CBSolverPlugin *>::iterator it = Base::plugins_.begin(); it != Base::plugins_.end(); it++)
        (*it)->Apply(time);
    
    SNESDestroy(&snes_);
    InitPETScSolver();
    
    snesStep_ = 0;
    SNESSolve(snes_, PETSC_NULL, displacement_);
    
    SNESGetIterationNumber(snes_, &snesIts);
    KSPGetTotalIterations(ksp_, &kspIts);
    SNESGetConvergedReason(snes_, &snesReason);
    KSPGetConvergedReason(ksp_, &kspReason);
    
    DCCtrl::debug << "\n\t--- --- Solver step finished --- ---\n";
    DCCtrl::debug << "\tSNES Reason: "     << Petsc::SNESReasonToString(snesReason) << "\n";
    DCCtrl::debug << "\tSNES Iterations: " << snesIts    << "\n";
    DCCtrl::debug << "\tKSP  Reason: "     << Petsc::KSPReasonToString(kspReason)  << "\n";
    DCCtrl::debug << "\tKSP  Iterations: " << kspIts     << "\n\n";
    Base::convergedReason_ = snesReason;
    Base::snesIterations_ += snesIts;
    Base::kspIterations_ += kspIts;
    
    if ((snesReason <= 0) || (kspReason < 0)) {
        VecSet(displacement_, 0);
        return CBStatus::FAILED;
    }
    
    // ---------------- Give the plugins the chance to analyse the results and to share their honest opinions --------------
    
    bool evaluate = false;
    
    for (auto p : plugins_)
        if (p->WantsToAnalyzeResults())
            evaluate = true;
    
    if (evaluate) {
        Vec localDisplacedNodesSeq = 0;
        VecCopy(Base::nodes_, tmpVector_);
        VecAXPY(tmpVector_, 1, displacement_);
        
        if (DCCtrl::IsParallel()) {
            VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
            Base::adapter_->LinkNodes(localDisplacedNodesSeq);
        } else {
            Base::adapter_->LinkNodes(tmpVector_);
        }
        
        for (auto p : plugins_)
            p->AnalyzeResults();
    }
    
    // --------------------------------------------------
    
    CBStatus pluginsFeedback = CBStatus::DACCORD;
    
    for (auto p : plugins_) {
        CBStatus f = p->GetStatus();
        
        if ((f == CBStatus::REPEAT) && (pluginsFeedback != CBStatus::FAILED))
            pluginsFeedback = CBStatus::REPEAT;
        
        if (f == CBStatus::FAILED)
            pluginsFeedback = CBStatus::FAILED;
        
        if (f == CBStatus::PREPARING_SIMULATION)
            pluginsFeedback = CBStatus::PREPARING_SIMULATION;
    }
    
    switch (pluginsFeedback) {
        case CBStatus::FAILED:
        case CBStatus::REPEAT:
            return pluginsFeedback;
            
        default:
            
            // Apply displacements to global nodes vector
            VecAXPY(Base::nodes_, 1, displacement_);
            
            return CBStatus::SUCCESS;
    }
} // CBSolverEquilibrium::SolverStep
