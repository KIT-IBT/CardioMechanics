/*
 *  CBSolverNewmarkBeta.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 29.09.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBSolverNewmarkBeta.h"

PetscErrorCode CBSolverNewmarkBetaSNESHelperFunctionForces(SNES snes, Vec u, Vec f, void *_solver) {
    CBSolverNewmarkBeta *solver = reinterpret_cast<CBSolverNewmarkBeta *>(_solver);
    CBStatus             rc     = solver->CalcNodalForces(u, f);
    
    if (rc != CBStatus::SUCCESS)
        SNESSetFunctionDomainError(snes);
    
    return 0;
}

PetscErrorCode CBSolverNewmarkBetaSNESHelperFunctionForcesJacobian(SNES snes, Vec u, Mat jacobian,
                                                                   Mat preconditionerMatrix, void *_solver) {
    CBSolverNewmarkBeta *solver = reinterpret_cast<CBSolverNewmarkBeta *>(_solver);
    CBStatus             rc     = solver->CalcNodalForcesJacobian(u, jacobian);
    
    //    CBStatus             rc     = solver->CalcNodalForcesJacobianAndDamping(u, jacobian);
    
    if (rc != CBStatus::SUCCESS)
        SNESSetFunctionDomainError(snes);
    
    return 0;
}

void CBSolverNewmarkBeta::Init(ParameterMap *_parameter, CBModel *_model) {
    Base::Init(_parameter, _model);
    Base::status_ = CBStatus::WAITING;
}

void CBSolverNewmarkBeta::DeInit() {
    VecDestroy(&residuum_);
    VecDestroy(&displacement_);
    VecDestroy(&absDisplacement_);
    VecDestroy(&tmpDisplacement_);
    VecDestroy(&tmpVelocity_);
    
    VecDestroy(&velocity_);
    VecDestroy(&acceleration_);
    VecDestroy(&tmpVector_);
    SNESDestroy(&snes_);
    KSPDestroy(&ksp_);
    PCDestroy(&pc_);
    Base::DeInit();
}

void CBSolverNewmarkBeta::InitVectors() {
    Base::InitNodalForces();
    
    VecDuplicate(Base::nodalForces_, &residuum_);
    VecDuplicate(Base::nodes_, &displacement_);
    VecDuplicate(Base::nodes_, &absDisplacement_);
    VecDuplicate(Base::nodes_, &initialGuess_);
    VecDuplicate(Base::nodes_, &tmpDisplacement_);
    VecDuplicate(Base::nodes_, &tmpVelocity_);
    VecDuplicate(Base::nodes_, &velocity_);
    VecDuplicate(Base::nodes_, &acceleration_);
    VecDuplicate(Base::nodes_, &tmpVector_);
    VecZeroEntries(residuum_);
    VecZeroEntries(displacement_);
    VecZeroEntries(absDisplacement_);
    VecZeroEntries(initialGuess_);
    VecZeroEntries(tmpDisplacement_);
    VecZeroEntries(tmpVelocity_);
    VecZeroEntries(velocity_);
    VecZeroEntries(acceleration_);
    VecZeroEntries(tmpVector_);
}

void CBSolverNewmarkBeta::SetVelocity(std::vector<Vector3<TFloat>> vel) {
    PetscInt from, to;
    
    VecGetOwnershipRange(velocity_, &from, &to);
    from /= 3;
    to   /= 3;
    
    for (int i = 0; i < (to-from); i++) {
        int cur = i+from;
        TFloat v[3];
        v[0] = vel[cur](0);
        v[1] = vel[cur](1);
        v[2] = vel[cur](2);
        int index[3];
        index[0] = 3*i+0;
        index[1] = 3*i+1;
        index[2] = 3*i+2;
        
        VecSetValues(velocity_, 3, index, v, INSERT_VALUES);
    }
}

void CBSolverNewmarkBeta::SetAcceleration(std::vector<Vector3<TFloat>> acc) {
    PetscInt from, to;
    
    VecGetOwnershipRange(acceleration_, &from, &to);
    from /= 3;
    to   /= 3;
    
    for (int i = 0; i < (to-from); i++) {
        int cur = i+from;
        TFloat v[3];
        v[0] = acc[cur](0);
        v[1] = acc[cur](1);
        v[2] = acc[cur](2);
        int index[3];
        index[0] = 3*i+0;
        index[1] = 3*i+1;
        index[2] = 3*i+2;
        
        VecSetValues(acceleration_, 3, index, v, INSERT_VALUES);
    }
}

void CBSolverNewmarkBeta::InitMatrices() {
    Base::InitNodalForcesJacobian();
    
    InitMassMatrix();
    InitDampingMatrix();
}

void CBSolverNewmarkBeta::InitMassMatrix() {
    if (useConsistentMassMatrix_) {
        InitMassMatrixConsistent();
    } else {
        InitMassMatrixLumped();
    }
    isInitMassMatrixDone_ = true;
    
    //    DCCtrl::print << "MassMatrix >>>\n";
    //    PetscViewer viewer;
    //    MatView(massMatrix_, PETSC_VIEWER_STDOUT_SELF);
    //    PetscViewerDestroy(&viewer);
    //    DCCtrl::print << "<<< MassMatrix\n";
}

void CBSolverNewmarkBeta::InitPETScSolver() {
    // Set up the linear model of the nonlinear equations
    // r + Ad = 0
    // m1 = 1 / (beta * dt^2)
    // c1 = gamma / (beta * dt)
    SNESCreate(PETSC_COMM_WORLD, &snes_);
    
    // function to calc residual r(d_n+1) = 0 = M * [m1 * (d_n+1 - d~_n+1)] + C * [c1 * (d_n+1 - d~_n+1) + v~_n+1] + f_int
    // - f_ext
    SNESSetFunction(snes_, residuum_, CBSolverNewmarkBetaSNESHelperFunctionForces, (void *)this);
    
    // function to calc jacobian A = m1 * M + c1 * C + K_int - K_ext
    SNESSetJacobian(snes_, Base::nodalForcesJacobian_, Base::nodalForcesJacobian_,
                    CBSolverNewmarkBetaSNESHelperFunctionForcesJacobian, (void *)this);
    
    SNESSetTolerances(snes_, precision_, precision_, precision_, maxSnesIts_, maxFunEval_);
    
    SNESGetKSP(snes_, &ksp_);
    KSPGetPC(ksp_, &pc_);
    
    // Set how often the Jacobian and Preconditioner are rebuilt:
    // -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2
    // means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again
    SNESSetLagJacobian(snes_, -2);
    SNESSetLagPreconditioner(snes_, -2);
    
    if (Base::parameters_->Get<bool>("Solver.LU", true)) {
        PCSetType(pc_, PCLU);
        KSPSetType(ksp_, "preonly");
    }
    
    if (solverType_ == "mumps") {
        PCFactorSetMatSolverType(pc_, "mumps");
    } else if (solverType_ == "superlu") {
        if (DCCtrl::IsParallel())
            PCFactorSetMatSolverType(pc_, "superlu_dist");
        else
            PCFactorSetMatSolverType(pc_, "superlu");
    } else {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitPETScSolver(): unkown Solver.Type " + solverType_ +
                                 ". Choose either mumps or superlu.");
    }
    
    SNESSetFromOptions(snes_);
    KSPSetFromOptions(ksp_);
    PCSetFromOptions(pc_);
}  // CBSolverNewmarkBeta::InitPETScSolver

void CBSolverNewmarkBeta::InitParameters() {
    Base::InitParameters();
    
    InitDampingParameter();
    InitNewmarkBetaParameter();
}

void CBSolverNewmarkBeta::InitNewmarkBetaParameter() {
    // gamma > 0.5 enables numerical dissipation, however it degrades the accuracy of the solver
    beta_                    = Base::parameters_->Get<double>("Solver.NewmarkBeta.Beta", 0.25);
    gamma_                   = Base::parameters_->Get<double>("Solver.NewmarkBeta.Gamma", 0.5);
    useConsistentMassMatrix_ = Base::parameters_->Get<bool>("Solver.NewmarkBeta.ConsistentMassMatrix", true);
    solverType_              = Base::parameters_->Get<std::string>("Solver.NewmarkBeta.Type", "mumps");
}

void CBSolverNewmarkBeta::InitDampingMatrix() {
    if (!isInitDampingParametersDone_) {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitDampingMatrix(): void CBSolverNewmarkBeta::InitDampingParameter() has to be run first");
    }
    
    if (!Base::isInitNodalForcesJacobianDone_) {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitDampingMatrix(): void CBSolver::InitNodalForcesJacobian() has to be run first");
    }
    
    if (!Base::isInitFormulationDone_) {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitDampingMatrix(): void CBSolver::InitFormulation() has to be run first");
    }
    
    if (!isInitMassMatrixDone_) {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitDampingMatrix(): void CBSolver::InitMassMatrix() has to be run first");
    }
    
    if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0)) {
        // Calculate nodes forces jacobian and initialize the damping matrix
        
        if (DCCtrl::IsParallel()) {
            MatCreateAIJ(
                         Petsc::Comm(), 3 * numLocalNodes_, 3 * numLocalNodes_, PETSC_DETERMINE, PETSC_DETERMINE, 0,
                         model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, 0,
                         model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, &dampingMatrix_);
            MatSetLocalToGlobalMapping(dampingMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
        } else {
            MatCreateSeqAIJ(Petsc::Comm(), 3 * Base::numNodes_, 3 * Base::numNodes_, 0,
                            model_->GetNodeNeighborsForNnz().data(), &dampingMatrix_);
            MatSetLocalToGlobalMapping(dampingMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
        }
        
        UpdateGhostNodesAndLinkToAdapter();
        
        if (globalRayleighBeta_ != 0) {
            Base::adapter_->LinkNodalForcesJacobian(dampingMatrix_);
            Base::formulation_->CalcNodalForcesJacobian();
            
            MatAssemblyBegin(dampingMatrix_, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dampingMatrix_, MAT_FINAL_ASSEMBLY);
            
            if (globalRayleighAlpha_ != 0) {
                MatAXPY(dampingMatrix_, globalRayleighAlpha_/globalRayleighBeta_, massMatrix_, DIFFERENT_NONZERO_PATTERN);
                MatScale(dampingMatrix_, globalRayleighBeta_);
            } else {
                MatScale(dampingMatrix_, globalRayleighBeta_);
            }
        } else {
            MatCopy(massMatrix_, dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
            MatScale(dampingMatrix_, globalRayleighAlpha_);
        }
        
        MatSetLocalToGlobalMapping(dampingMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
        MatAssemblyBegin(dampingMatrix_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(dampingMatrix_, MAT_FINAL_ASSEMBLY);
    }
}  // CBSolverNewmarkBeta::InitDampingMatrix

void CBSolverNewmarkBeta::InitDampingParameter() {
    std::string dampingType = Base::parameters_->Get<std::string>("Material.Global.Damping.Type", "Rayleigh");
    
    if (dampingType == "Rayleigh") {
        globalRayleighAlpha_ = Base::parameters_->Get<double>("Materials.Global.Damping.Rayleigh.Alpha", 0);
        globalRayleighBeta_  = Base::parameters_->Get<double>("Materials.Global.Damping.Rayleigh.Beta", 0);
    } else {
        throw std::runtime_error(
                                 "CBSolverNewmarkBeta::InitGlobalDampingParameter(): Global damping type " + dampingType + " is unkown !!");
    }
    
    isInitDampingParametersDone_ = true;
}

/// Write matrix, vector and initial guess of the ksp solver to files. These
/// can be useful for finding good solver parameters using some external tools.
void CBSolverNewmarkBeta::ExportSNESMatrix(TFloat time) {
    std::string exportDir = model_->GetExporter()->GetExportDirPrefix();
    std::string SNESMatrixDir = exportDir + "_SNES";
    
    if (!frizzle::filesystem::CreateDirectory(SNESMatrixDir) && (time == 0)) {
        throw std::runtime_error(
                                 "void CBSolverNewmarkBeta::ExportSNESMatrix(): Path " + SNESMatrixDir + " exists but is not a directory");
    }
    
    std::string AmatFilename            = SNESMatrixDir + "/Amat." + std::to_string(time) + ".bin";
    std::string rhsvecFilename          = SNESMatrixDir + "/rhs." + std::to_string(time) + ".bin";
    std::string initialGuessFilename    = SNESMatrixDir + "/initialGuess." + std::to_string(time) + ".bin";
    
    DCCtrl::debug << "AmatFilename: " << AmatFilename << std::endl;
    DCCtrl::debug << "rhsvecFilename: " << rhsvecFilename << std::endl;
    DCCtrl::debug << "initialGuessFilename: " << initialGuessFilename << std::endl;
    
    
    SNESGetKSP(snes_, &ksp_);
    Mat Amat, Pmat;
    Vec rhsvec;
    
    KSPGetOperators(ksp_, &Amat, &Pmat);
    
    PetscViewer viewMatrix;
    PetscViewer viewVector;
    
    DCCtrl::debug << "Saving matrix to file " << AmatFilename << "...\n";
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, AmatFilename.c_str(), FILE_MODE_WRITE, &viewMatrix);
    MatView(Amat, viewMatrix);
    DCCtrl::debug << "Saving matrix to file " << AmatFilename << "...OK\n";
    
    KSPGetRhs(ksp_, &rhsvec);
    
    DCCtrl::debug << "Saving vector rhs to file " << rhsvecFilename << "...\n";
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, rhsvecFilename.c_str(), FILE_MODE_WRITE, &viewVector);
    VecView(rhsvec, viewVector);
    DCCtrl::debug << "Saving vector rhs to file " << rhsvecFilename << "...OK\n";
    
    DCCtrl::debug << "Saving vector initialGuess to file " << initialGuessFilename << "...\n";
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, initialGuessFilename.c_str(), FILE_MODE_WRITE, &viewVector);
    VecView(initialGuess_, viewVector);
    DCCtrl::debug << "Saving vector initialGuess to file " << initialGuessFilename << "...OK\n";
    
    PetscViewerDestroy(&viewMatrix);
    PetscViewerDestroy(&viewVector);
}  // CBSolverNewmarkBeta::ExportSNESMatrix

CBStatus CBSolverNewmarkBeta::CalcNodalForces(Vec displacement, Vec forces) {
    snesStep_++;
    DCCtrl::debug << "SNES Step:" << snesStep_ << " - F-Norm: ";
    
    Vec localDisplacedNodesSeq = 0;
    Vec velocity               = 0;
    
    // Set node coordinates to d~_n+1
    VecCopy(Base::nodes_, tmpVector_);
    VecAXPY(tmpVector_, 1, displacement);
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
        Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    } else {
        Base::adapter_->LinkNodes(tmpVector_);
    }
    
    // add internal nodal forces f_int (passive and active stress contribution of the tissue)
    VecZeroEntries(Base::nodalForces_);
    Base::adapter_->LinkNodalForces(Base::nodalForces_);
    CBStatus rc = Base::formulation_->CalcNodalForces();
    
    if (rc != CBStatus::SUCCESS)
        DCCtrl::cverbose << CBStatusToStr(rc);
    
    // add external nodal forces -f_ext (contribution of plugins)
    for (auto p : plugins_)
        p->ApplyToNodalForces();
    
    VecAssemblyBegin(Base::nodalForces_);
    VecAssemblyEnd(Base::nodalForces_);
    
    // add mass M * [m1 * (d_n+1 - d~_n+1)] with m1 = 1 / (beta * dt^2)
    VecCopy(displacement, tmpVector_);
    VecAXPY(tmpVector_, -1, tmpDisplacement_);
    PetscScalar factor = 1 / (beta_*timing_.GetTimeStep()*timing_.GetTimeStep());
    VecScale(tmpVector_, factor);
    MatMult(massMatrix_, tmpVector_, forces);
    
    // add damping forces C * [c1 * (d_n+1 - d~_n+1) + v~_n+1] with c1 = gamma / (beta * dt)
    if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0)) {
        VecDuplicate(tmpVelocity_, &velocity);
        VecCopy(tmpVelocity_, velocity);
        VecAXPY(velocity, gamma_ * timing_.GetTimeStep(), tmpVector_);
        MatMult(dampingMatrix_, velocity, tmpVector_);
        VecAXPY(forces, 1, tmpVector_);
    }
    
    VecAXPY(forces, 1, Base::nodalForces_);
    
    VecDestroy(&velocity);
    
    if (localDisplacedNodesSeq)
        VecDestroy(&localDisplacedNodesSeq);
    
    PetscScalar r, s;
    
    // Calc ||F|| and ||x|| of SNES Step to print
    VecNorm(forces, NORM_2, &r);
    VecNorm(displacement_, NORM_2, &s);
    DCCtrl::debug << std::setprecision(12) << std::fixed << r << " - S-Norm: " << s << "\n";
    
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
}  // CBSolverNewmarkBeta::CalcNodalForces

CBStatus CBSolverNewmarkBeta::CalcNodalForcesJacobianAndDamping(Vec displacement, Mat jacobian) {
    DCCtrl::print << "CALC JACOBIAN AND DAMPING MAT" << std::endl;
    
    double t = MPI_Wtime();
    DCCtrl::debug<< "\nCalculating System Matrix:\n";
    
    if (displacement != 0)
        VecCopy(displacement, tmpVector_);
    VecZeroEntries(tmpVector_);
    
    Base::adapter_->SetFiniteDifferencesEpsilon(epsilon_);
    
    Vec localDisplacedNodesSeq = 0;
    
    VecCopy(Base::nodes_, tmpVector_);
    if (displacement != 0)
        VecAXPY(tmpVector_, 1, displacement);
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
        Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    } else {
        Base::adapter_->LinkNodes(tmpVector_);
    }
    
    MatZeroEntries(dampingMatrix_);
    MatSetLocalToGlobalMapping(dampingMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    Base::adapter_->LinkNodalForcesJacobian(dampingMatrix_);
    
    double t1 = MPI_Wtime();
    DCCtrl::debug << "Elements ...";
    CBStatus rc = Base::formulation_->CalcNodalForcesJacobian();
    if (rc != CBStatus::SUCCESS)
        DCCtrl::cverbose << CBStatusToStr(rc);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "Assembling ...";
    MatAssemblyBegin(dampingMatrix_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dampingMatrix_, MAT_FINAL_ASSEMBLY);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/tmp/dampingMatrix.m", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    std::cout << "\nwriting damping matrix to dampingMatrix.m ..." << std::flush;
    MatView(dampingMatrix_, viewer);
    std::cout << " done!\n";
    PetscViewerDestroy(&viewer);
    
    MatZeroEntries(jacobian);
    MatSetLocalToGlobalMapping(jacobian, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    Base::adapter_->LinkNodalForcesJacobian(jacobian);
    
    for (std::vector<CBSolverPlugin *>::iterator it = Base::plugins_.begin(); it != Base::plugins_.end(); it++) {
        t1 = MPI_Wtime();
        DCCtrl::debug << "Plugin: " <<  (*it)->GetName() << " ...";
        (*it)->ApplyToNodalForcesJacobian();
        DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    }
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "Assembling ...";
    MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/tmp/jacobian1.m", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    std::cout << "\nwriting jacobian to jacobian1.m ..." << std::flush;
    MatView(jacobian, viewer);
    std::cout << " done!\n";
    PetscViewerDestroy(&viewer);
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "Assembling ...";
    MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    MatAXPY(jacobian, 1.0, dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
    
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/tmp/jacobian2.m", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    std::cout << "\nwriting jacobian to jacobian2.m ..." << std::flush;
    MatView(jacobian, viewer);
    std::cout << " done!\n";
    PetscViewerDestroy(&viewer);
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "Damping matrix ... ";
    
    // -----
    if (globalRayleighAlpha_ != 0) {
        // C = rayleighAlpha * M + rayleighBeta * K
        MatAXPY(dampingMatrix_, globalRayleighAlpha_/globalRayleighBeta_, massMatrix_, DIFFERENT_NONZERO_PATTERN);
        MatScale(dampingMatrix_, globalRayleighBeta_);
    } else {
        MatScale(dampingMatrix_, globalRayleighBeta_);
    }
    
    // -----
    
    if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0)) {
        // A += gamma/(dt*beta) * C
        MatAXPY(jacobian, gamma_ / (beta_*timing_.GetTimeStep()), dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
    }
    
    // A += 1/(dt^2*beta) * M
    MatAXPY(jacobian, 1 / (beta_*timing_.GetTimeStep()*timing_.GetTimeStep()), massMatrix_, DIFFERENT_NONZERO_PATTERN);
    
    MatAXPY(jacobian, 1, boundaryConditionsNodalForcesJacobianDiagonalComponents_, DIFFERENT_NONZERO_PATTERN);
    
    MatSetLocalToGlobalMapping(jacobian, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    if (localDisplacedNodesSeq)
        VecDestroy(&localDisplacedNodesSeq);
    
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
}  // CBSolverNewmarkBeta::CalcNodalForcesJacobianAndDamping

CBStatus CBSolverNewmarkBeta::CalcNodalForcesJacobian(Vec displacement, Mat jacobian) {
    // DCCtrl::print << "CALC JACOBIAN" << std::endl;
    
    double t = MPI_Wtime();
    
    DCCtrl::debug<< "\n\nCalculating System Matrix:\n";
    
    if (displacement != 0)
        VecCopy(displacement, tmpVector_);
    VecZeroEntries(tmpVector_);
    
    Base::adapter_->SetFiniteDifferencesEpsilon(CalcFiniteDifferencesEpsilon(tmpVector_));
    
    Vec localDisplacedNodesSeq = 0;
    
    // Set nodes to d~_n+1
    VecCopy(Base::nodes_, tmpVector_);
    if (displacement != 0)
        VecAXPY(tmpVector_, 1, displacement);
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
        Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    } else {
        Base::adapter_->LinkNodes(tmpVector_);
    }
    
    MatZeroEntries(jacobian);
    
    MatSetLocalToGlobalMapping(jacobian, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    Base::adapter_->LinkNodalForcesJacobian(jacobian);
    
    // add the derivative of the internal nodal forces f_int w.r.t. current displacement
    double t1 = MPI_Wtime();
    DCCtrl::debug << "Elements ...";
    CBStatus rc = Base::formulation_->CalcNodalForcesJacobian();
    
    if (rc != CBStatus::SUCCESS)
        DCCtrl::cverbose << CBStatusToStr(rc);
    
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    // add the derivative of the external nodal forces -f_ext w.r.t. current displacement
    for (auto p : plugins_) {
        t1 = MPI_Wtime();
        DCCtrl::debug << "Plugin: " <<  p->GetName() << " ...";
        p->ApplyToNodalForcesJacobian();
        DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    }
    
    t1 = MPI_Wtime();
    DCCtrl::debug << "Assembling ...";
    MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    t1 = MPI_Wtime();
    
    DCCtrl::debug << "Damping matrix ... ";
    
    if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0)) {
        // A += gamma/(dt*beta) * C
        MatAXPY(jacobian, gamma_ / (beta_*timing_.GetTimeStep()), dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
    }
    
    // A += 1/(dt^2*beta) * M
    MatAXPY(jacobian, 1 / (beta_*timing_.GetTimeStep()*timing_.GetTimeStep()), massMatrix_, DIFFERENT_NONZERO_PATTERN);
    
    MatAXPY(jacobian, 1, boundaryConditionsNodalForcesJacobianDiagonalComponents_, DIFFERENT_NONZERO_PATTERN);
    MatSetLocalToGlobalMapping(jacobian, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    
    DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
    
    if (localDisplacedNodesSeq)
        VecDestroy(&localDisplacedNodesSeq);
    
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
}  // CBSolverNewmarkBeta::CalcNodalForcesJacobian

CBStatus CBSolverNewmarkBeta::CalcDampingMatrix() {
    if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0)) {
        CBStatus rc = CBStatus::FAILED;
        
        if (globalRayleighBeta_ != 0) {
            Base::adapter_->LinkNodalForcesJacobian(dampingMatrix_);
            rc = Base::formulation_->CalcNodalForcesJacobian();
            
            MatAssemblyBegin(dampingMatrix_, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dampingMatrix_, MAT_FINAL_ASSEMBLY);
            
            if (globalRayleighAlpha_ != 0) {
                // C = rayleighAlpha * M + rayleighBeta * K
                // PetscErrorCode MatAXPY(Mat Y,PetscScalar a,Mat X,MatStructure str) -> Y = aX + Y
                // dampingMatrix_  = dampingMatrix_ + (globalRayleighAlpha_/globalRayleighBeta_) * massMatrix_
                // MatAXPY(dampingMatrix_, globalRayleighAlpha_/globalRayleighBeta_, massMatrix_, DIFFERENT_NONZERO_PATTERN);
                // dampingMatrix_ *= globalRayleighBeta_ * dampingMatrix_
                MatScale(dampingMatrix_, globalRayleighBeta_);
                MatAXPY(dampingMatrix_, globalRayleighAlpha_, massMatrix_, DIFFERENT_NONZERO_PATTERN);
                
                // dampingMatrix_ *= globalRayleighBeta_ * (K + (globalRayleighAlpha_/globalRayleighBeta_) * massMatrix_)
                //                  - K*beta + alpha*M
                
                //       MatScale(dampingMatrix_, globalRayleighBeta_);
                //     MatAXPY(dampingMatrix_, globalRayleighAlpha_, massMatrix_, DIFFERENT_NONZERO_PATTERN);
                //              MatCopy(massMatrix_, dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
                //              MatScale(dampingMatrix_, globalRayleighAlpha_);
            } else {
                MatScale(dampingMatrix_, globalRayleighBeta_);
            }
            rc = CBStatus::SUCCESS;
        } else {
            MatCopy(massMatrix_, dampingMatrix_, DIFFERENT_NONZERO_PATTERN);
            MatScale(dampingMatrix_, globalRayleighAlpha_);
            rc = CBStatus::SUCCESS;
        }
        
        
        MatSetLocalToGlobalMapping(dampingMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
        MatAssemblyBegin(dampingMatrix_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(dampingMatrix_, MAT_FINAL_ASSEMBLY);
        
        return rc;
    }
    return CBStatus::NOTHING_DONE;
}  // CBSolverNewmarkBeta::CalcDampingMatrix

void CBSolverNewmarkBeta::InitMassMatrixLumped() {
    if (DCCtrl::IsParallel()) {
        MatCreateAIJ(
                     Petsc::Comm(), 3 * Base::numLocalNodes_, 3 * Base::numLocalNodes_, PETSC_DETERMINE, PETSC_DETERMINE, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, &massMatrix_);
        MatSetLocalToGlobalMapping(massMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    } else {
        MatCreateSeqAIJ(Petsc::Comm(), 3 * Base::numNodes_, 3 * Base::numNodes_, 0,
                        model_->GetNodeNeighborsForNnz().data(), &massMatrix_);
        MatSetLocalToGlobalMapping(massMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    }
    
    Base::adapter_->LinkMassMatrix(massMatrix_);
    
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    
    Vec localDisplacedNodesSeq = 0;
    
    if (DCCtrl::IsParallel()) {
        VecCopy(Base::nodes_, tmpVector_);
        
        if (DCCtrl::IsParallel()) {
            VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        }
        
        VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
        Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    } else {
        Base::adapter_->LinkNodes(Base::nodes_);
    }
    
    for (auto &it : solidElements_)
        it->CalcLumpedMassMatrix();
    
    MatAssemblyBegin(massMatrix_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(massMatrix_, MAT_FINAL_ASSEMBLY);
    
    if (localDisplacedNodesSeq)
        VecDestroy(&localDisplacedNodesSeq);
}  // CBSolverNewmarkBeta::InitMassMatrixLumped

void CBSolverNewmarkBeta::InitMassMatrixConsistent() {
    if (DCCtrl::IsParallel()) {
        MatCreateAIJ(
                     Petsc::Comm(), 3 * Base::numLocalNodes_, 3 * Base::numLocalNodes_, PETSC_DETERMINE, PETSC_DETERMINE, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, &massMatrix_);
        MatSetLocalToGlobalMapping(massMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    } else {
        MatCreateSeqAIJ(Petsc::Comm(), 3 * Base::numNodes_, 3 * Base::numNodes_, 0,
                        model_->GetNodeNeighborsForNnz().data(), &massMatrix_);
        MatSetLocalToGlobalMapping(massMatrix_, Base::nodesIndicesMapping_, Base::nodesIndicesMapping_);
    }
    
    Base::adapter_->LinkMassMatrix(massMatrix_);
    
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    Vec localDisplacedNodesSeq = 0;
    if (DCCtrl::IsParallel()) {
        VecCopy(Base::nodes_, tmpVector_);
        
        if (DCCtrl::IsParallel()) {
            VecGhostUpdateBegin(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(tmpVector_, INSERT_VALUES, SCATTER_FORWARD);
        }
        
        VecGhostGetLocalForm(tmpVector_, &localDisplacedNodesSeq);
        Base::adapter_->LinkNodes(localDisplacedNodesSeq);
    } else {
        Base::adapter_->LinkNodes(Base::nodes_);
    }
    
    for (auto &it : solidElements_) {
        it->CalcConsistentMassMatrix();
    }
    MatAssemblyBegin(massMatrix_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(massMatrix_, MAT_FINAL_ASSEMBLY);
    if (localDisplacedNodesSeq)
        VecDestroy(&localDisplacedNodesSeq);
}  // CBSolverNewmarkBeta::InitMassMatrixConsistent

void CBSolverNewmarkBeta::Export(PetscScalar time) {
    Base::Export(time);
    
    auto *exporter = model_->GetExporter();
    std::string str = "KineticEnergy";
    if (exporter->GetExportOption(str, false)) {
        model_->SetGlobalData(str, GetKineticEnergy());
    }
    
    str = "DampingEnergyDissipation";
    if (exporter->GetExportOption(str, false)) {
        model_->SetGlobalData(str, GetDampingEnergyDissipation());
    }
    
    str = "TotalEnergy";
    if (exporter->GetExportOption(str, false)) {
        model_->SetGlobalData(str, GetDampingEnergyDissipation()+GetDeformationEnergy()+GetKineticEnergy());
    }
    
    str = "Velocity";
    if (exporter->GetExportOption(str, true))
        Base::ExportNodesVectorData(str, velocity_);
    
    str = "Acceleration";
    if (exporter->GetExportOption(str, true))
        Base::ExportNodesVectorData(str, acceleration_);
    
    str = "Displacement";
    if (exporter->GetExportOption(str, true))
        Base::ExportNodesVectorData(str, displacement_);
    
    str = "AbsDisplacement";
    if (exporter->GetExportOption(str, true))
        Base::ExportNodesVectorData(str, absDisplacement_);
    
    str = "NodalForces";
    if (exporter->GetExportOption(str, false))
        Base::ExportNodesVectorData(str, Base::nodalForces_);
}  // CBSolverNewmarkBeta::Export

CBStatus CBSolverNewmarkBeta::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation) {
    SNESConvergedReason snesReason;
    KSPConvergedReason  kspReason;
    PetscInt            kspIts = 0;
    
    VecAssemblyBegin(Base::nodes_);
    VecAssemblyEnd(Base::nodes_);
    
    UpdateGhostNodesAndLinkToAdapter();
    
    for (auto p : plugins_)
        p->Apply(time);
    
    // Predictor phase (see TJR Hughes 1978: Implicit-explicit Finite elements in nonlinear transient analysis):
    // tmpDisplacement_ = d~_n+1 = dt * v_n + dt^2/2 * (1-2*beta) * a_n     (+ d_n, but displacements are relative to the
    // current configuration stored in Base::nodes_)
    VecAXPBYPCZ(tmpDisplacement_, timing_.GetTimeStep(),
                timing_.GetTimeStep() * timing_.GetTimeStep() * 0.5 * (1 - 2 * beta_), 0, velocity_, acceleration_);
    
    // tmpVelocity_     = v~_n+1 = v_n + dt * (1-gamma) * a_n
    VecAXPBYPCZ(tmpVelocity_, 1.0, (1 - gamma_) * timing_.GetTimeStep(), 0, velocity_, acceleration_);
    
    // update displacement_ with tmpDisplacement_ since we need to explicitly set an initial value for SNESSolve
    VecCopy(tmpDisplacement_, displacement_);
    
    // initialGuess_ = displacement_   -> save initialGuess to be exported by ExportSNESMatrix
    VecCopy(displacement_, initialGuess_);
    
    if ((std::abs(time-prevTime_) > std::numeric_limits<float>::epsilon()) || forceJacobianAndDampingRecalculation ||
        updateJacobian_) {
        // The damping matrix C has to be calculated in each time step, as C depends on the changing stiffness matrix Ks
        TFloat t1 = MPI_Wtime();
        DCCtrl::debug << "\n\nCalc Damping Matrix ...";
        CalcDampingMatrix();
        DCCtrl::debug << " done [" << MPI_Wtime() - t1 << " s]" << std::endl;
        
        prevTime_ = time;
    }
    
    bool newJacobian = false;
    if (updateJacobian_ || (time < prevTime_) || forceJacobianAndDampingRecalculation) {
        DCCtrl::debug << "Calculating new Jacobian ...";
        snesStep_ = 0;
        SNESDestroy(&snes_);
        InitPETScSolver();
        DCCtrl::debug << "done" << std::endl;
        updateJacobian_ = false;
        newJacobian     = true;
    } else {
        DCCtrl::debug << "Reusing old Jacobian." << std::endl;
    }
    
    // Solve following nonlinear system of equations r + A * dx = 0 with
    // r(d_n+1) = 0 = M * [m1 * (d_n+1 - d~_n+1)] + C * [c1 * (d_n+1 - d~_n+1) + v~_n+1] + f_int - f_ext
    // with Newton's method using the Jacobian:
    // A = m1 * M + c1 * C + K_int - K_ext
    // m1 = 1 / (beta * dt^2)
    // c1 = gamma / (beta * dt)
    // Note that f_tot = f_int + f_ext in our case, since all plugins that contribute to f_ext return negative values
    // In each iteration, the internal and external forces are linearized:
    // K * d_n+1 ~= f_int(d_n+1, t_n+1) + f_ext(d_n+1, t_n+1)
    SNESSolve(snes_, PETSC_NULL, displacement_);
    
    SNESGetIterationNumber(snes_, &snesIts_);
    KSPGetTotalIterations(ksp_, &kspIts);
    SNESGetConvergedReason(snes_, &snesReason);
    KSPGetConvergedReason(ksp_, &kspReason);
    
    DCCtrl::debug << "\n--- --- Solver step finished --- ---\n";
    DCCtrl::debug << "SNES Reason: "     << Petsc::SNESReasonToString(snesReason) << "\n";
    DCCtrl::debug << "SNES Iterations: " << snesIts_    << "\n";
    DCCtrl::debug << "KSP  Reason: "     << Petsc::KSPReasonToString(kspReason)  << "\n";
    DCCtrl::debug << "KSP  Iterations: " << kspIts     << "\n\n";
    Base::convergedReason_ = snesReason;
    Base::snesIterations_ += snesIts_;
    Base::kspIterations_  += kspIts;
    
    // Repeat & recalculate system matrix IF
    // 1) SNES diverges
    // 2) SNES converges with SNORM RELATIVE (if this is the case, the residual is often still significantly different
    // from 0, which indicates failure in the machanics step) -> enable by adding:  || (snesReason == 4)
    // 3) SNES took too many iterations to converge
    if ( (snesReason <= 0) && ((snesReason != -5) || !Base::ignoreMaxIt_) ) {
        DCCtrl::debug << "SNES Diverged: Forcing Jacobian calculation in the next time step" << std::endl;
        if (newJacobian) {
            updateJacobian_ = true;
            return CBStatus::FAILED;
        } else {
            updateJacobian_ = true;
            return CBStatus::REPEAT;
        }
    }
    
    if (snesIts_ > refSnesIts_) {
        DCCtrl::debug << "SNES Iterations: " << snesIts_ << " > " << refSnesIts_ <<
        " Forcing Jacobian calculation in the next time step" << std::endl;
        updateJacobian_ = true;
    }
    
    // ---------------- Give the plugins the chance to analyse the results and to share their honest opinions
    // --------------
    
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
            
            // Corrector phase:
            // Newton iterations (SNES) have calculated: displacement = d_n+1
            auto timestep = timing_.GetTimeStep();
            
            // acceleration_ = a_n+1 = (d_n+1 - d~_n+1) / (dt^2*beta)
            VecAXPBYPCZ(acceleration_, 1.0 / (beta_*timestep*timestep), -1.0 / (beta_*timestep*timestep), 0, displacement_,
                        tmpDisplacement_);
            
            // velocity_     = v_n+1 = v~_n+1 + dt * gamma * a_n+1
            VecCopy(tmpVelocity_, velocity_);
            VecAXPY(velocity_, gamma_ * timing_.GetTimeStep(), acceleration_);
            
            // Calculate kinetic energy and energy dissipated by damping
            Vec k;
            VecDuplicate(velocity_, &k);
            MatMult(massMatrix_, velocity_, k);
            VecDot(k, velocity_, &kineticEnergy_);
            kineticEnergy_ *= 0.5;
            if ((globalRayleighAlpha_ != 0) || (globalRayleighBeta_ != 0) ) {
                MatMult(dampingMatrix_, velocity_, k);
                PetscScalar ed;
                VecDot(k, displacement_, &ed);
                if (status_ != CBStatus::PREPARING_SIMULATION)
                    dampingEnergyDissipation_ += ed;
            }
            VecDestroy(&k);
            
            // Add displacement_ of relative configuration to absdisplacement_ of initial configuration
            if (time > 0.0) {
                VecAXPY(absDisplacement_, 1.0, displacement_);
            }
            
            // Apply displacements to global nodes vector
            if (evaluate)
                VecCopy(tmpVector_, Base::nodes_);
            else
                VecAXPY(Base::nodes_, 1, displacement_);
            
            return CBStatus::SUCCESS;
    }  // switch
}  // CBSolverNewmarkBeta::SolverStep

void CBSolverNewmarkBeta::SetZeroVelocityAndAcceleration() {
    VecSet(acceleration_, 0.0);
    VecSet(velocity_, 0.0);
}

void CBSolverNewmarkBeta::SetZeroDisplacement() {
    VecSet(absDisplacement_, 0.0);
}
