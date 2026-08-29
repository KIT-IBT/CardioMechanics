/*
 * File: CBSolverActiveStressEstimator.cpp
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


#include <algorithm>
#include <cstdlib>
#include <iomanip>

#include "CBSolverActiveStressEstimator.h"

void CBSolverActiveStressEstimator::Init(ParameterMap *parameters, CBModel *model) {
    CBSolverEquilibrium::Init(parameters, model);

    // The ContactHandling plugin provides the distance between the current (master) surface and the
    // prescribed target (slave) surface, which is the residual the estimator drives to zero.
    contact_ = nullptr;
    for (auto &p : plugins_)
        if (dynamic_cast<CBContactHandling *>(p) != nullptr)
            contact_ = dynamic_cast<CBContactHandling *>(p);
    if (contact_ == nullptr)
        throw std::runtime_error("CBSolverActiveStressEstimator::Init(): ContactHandling plugin is required");

    mat_ = parameters_->GetArray<TInt>("Solver.ActiveStressEstimator.Materials", {});

    std::string exportDir = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir();
    filename_       = parameters_->Get<std::string>("Solver.ActiveStressEstimator.ExportFile", exportDir + "/InvSlvNewtonIterNorm.dat");
    filenameLcurve_ = parameters_->Get<std::string>("Solver.ActiveStressEstimator.ExportLcurveDataFile", exportDir + "/LcurveData.dat");

    // Nodes and elements of interest: optionally restricted to a set of materials.
    std::set<TInt> nodesOfInterest;
    std::set<TInt> elementsOfInterest;
    for (auto &i : model_->GetSolidElements()) {
        if (mat_.size() != 0)
            if (std::find(mat_.begin(), mat_.end(), i->GetMaterialIndex()) == mat_.end())
                continue;
        for (unsigned int j = 0; j < i->GetNumberOfNodesIndices(); j++)
            nodesOfInterest.insert(i->GetNodeIndex(j));
        elementsOfInterest.insert(i->GetIndex());
    }

    numNodesOfInterestIndices_   = 3 * nodesOfInterest.size();
    numElementOfInterestIndices_ = elementsOfInterest.size();

    PetscInt *ni = new PetscInt[numNodesOfInterestIndices_];
    elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices_];

    nim_ = new PetscInt[numNodes_];
    for (int i = 0; i < numNodes_; i++)
        nim_[i] = -1;

    {
        auto it = nodesOfInterest.begin();
        for (size_t i = 0; i < nodesOfInterest.size(); i++) {
            nim_[*it] = i;
            ni[3*i + 0] = 3 * *it + 0;
            ni[3*i + 1] = 3 * *it + 1;
            ni[3*i + 2] = 3 * *it + 2;
            it++;
        }
    }

    TInt numTargets = parameters_->Get<TInt>("Solver.ActiveStressEstimator.NumberOfTargetNodes", -1);

    l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1", 0);
    l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2", 0);
    l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3", 0);
    l4_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l4", 0);
    l5_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l5", 0);
    l6_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l6", 0);
    l7_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l7", 0);

    srand(0); // deterministic subsampling of target nodes

    std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
    std::set<TInt> m2;
    for (auto i = m.begin(); i != m.end(); i++)
        if (nodesOfInterest.find(*i) != nodesOfInterest.end())
            m2.insert(*i);

    int r = numTargets > 0 ? m2.size() / numTargets : 0;
    if (r <= 1 || numTargets == -1) {
        mn_ = m2;
    } else {
        for (auto i : m2)
            if (rand() % r == 0)
                mn_.insert(i);
    }

    numMasterNodesIndices_ = 3 * mn_.size();
    masterNodesIndicesMapping_               = new PetscInt[numMasterNodesIndices_];
    masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

    {
        auto it = mn_.begin();
        for (size_t i = 0; i < mn_.size(); i++) {
            PetscInt j = nim_[*it];
            if (j == -1)
                throw std::runtime_error("CBSolverActiveStressEstimator::Init(): master node is not a node of interest");
            masterNodesIndicesMapping_[3*i + 0] = 3 * (*it) + 0;
            masterNodesIndicesMapping_[3*i + 1] = 3 * (*it) + 1;
            masterNodesIndicesMapping_[3*i + 2] = 3 * (*it) + 2;
            masterNodesIndicesNodesOfInterestMapping_[3*i + 0] = 3 * j + 0;
            masterNodesIndicesNodesOfInterestMapping_[3*i + 1] = 3 * j + 1;
            masterNodesIndicesNodesOfInterestMapping_[3*i + 2] = 3 * j + 2;
            it++;
        }
    }

    {
        auto it = elementsOfInterest.begin();
        for (size_t i = 0; i < elementsOfInterest.size(); i++) {
            elementsOfInterestMapping_[i] = *it;
            it++;
        }
    }

    ISCreateGeneral(MPI_COMM_SELF, numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
    ISCreateGeneral(MPI_COMM_SELF, numElementOfInterestIndices_, elementsOfInterestMapping_, PETSC_COPY_VALUES, &elementsOfInterestIndices_);
    ISCreateGeneral(MPI_COMM_SELF, numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
    VecDuplicate(nodes_, &dist_);

    GenerateElementLaplacian();

    // dfdtau gets one entry per (node DOF, adjacent solid element) pair, contributed by
    // CalcNodalForcesActiveStressJacobian. Count the adjacent elements per node so the matrix can be
    // preallocated exactly; the three DOFs of a node share the same count.
    dfdtauNnz_.assign(3 * numNodes_, 0);
    for (auto &e : GetSolidElementVector())
        for (unsigned int j = 0; j < e->GetNumberOfNodesIndices(); j++) {
            TInt n = e->GetNodeIndex(j);
            if ((n < 0) || (n >= numNodes_))
                throw std::runtime_error("CBSolverActiveStressEstimator::Init(): node index out of range");
            for (int k = 0; k < 3; k++)
                dfdtauNnz_[3*n + k]++;
        }

    DCPetsc::CreateSeqVector(numElementOfInterestIndices_, &ti_);
    VecDuplicate(ti_, &ti1_);
    VecDuplicate(ti_, &ti2_);
    VecDuplicate(nodes_, &tmpNodes_);

    delete[] ni;
}

void CBSolverActiveStressEstimator::UpdateMasterNodesOfInterest() {
    std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
    std::set<TInt> mn;
    for (auto i = m.begin(); i != m.end(); i++)
        if (mn_.find(*i) != mn_.end())
            mn.insert(*i);

    numMasterNodesIndices_ = 3 * mn.size();

    delete[] masterNodesIndicesMapping_;
    delete[] masterNodesIndicesNodesOfInterestMapping_;
    masterNodesIndicesMapping_                = new PetscInt[numMasterNodesIndices_];
    masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

    auto it = mn.begin();
    for (size_t i = 0; i < mn.size(); i++) {
        PetscInt j = nim_[*it];
        if (j == -1)
            throw std::runtime_error("CBSolverActiveStressEstimator::UpdateMasterNodesOfInterest(): master node is not a node of interest");
        masterNodesIndicesMapping_[3*i + 0] = 3 * (*it) + 0;
        masterNodesIndicesMapping_[3*i + 1] = 3 * (*it) + 1;
        masterNodesIndicesMapping_[3*i + 2] = 3 * (*it) + 2;
        masterNodesIndicesNodesOfInterestMapping_[3*i + 0] = 3 * j + 0;
        masterNodesIndicesNodesOfInterestMapping_[3*i + 1] = 3 * j + 1;
        masterNodesIndicesNodesOfInterestMapping_[3*i + 2] = 3 * j + 2;
        it++;
    }
    if (masterNodesIndices_)
        ISDestroy(&masterNodesIndices_);
    ISCreateGeneral(MPI_COMM_SELF, numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
}

CBStatus CBSolverActiveStressEstimator::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation) {
    Vec prevNodes;
    VecDuplicate(nodes_, &prevNodes);
    VecCopy(nodes_, prevNodes);

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    CBStatus rc = CBSolverEquilibrium::SolverStep(time);
    if (rc != CBStatus::SUCCESS) {
        VecDestroy(&prevNodes);
        return rc;
    }
    if (time == timing_.GetStartTime()) {
        VecDestroy(&prevNodes);
        return rc;
    }

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    VecZeroEntries(dist_);
    contact_->Apply(time);
    UpdateMasterNodesOfInterest();
    contact_->GetMasterNodesDistancesToSlaveElements(&dist_);

    Vec subDist;
    PetscScalar norm = 0;
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    VecNorm(subDist, NORM_2, &norm);
    VecRestoreSubVector(dist_, masterNodesIndices_, &subDist);
    DCCtrl::print << "Estimate " << norm << "\n";

    int i = 0;
    TFloat abs = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.AbsTol", 1e-4);
    TFloat rel = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.RelTol", 1e-4);

    while (i < 5 && norm > abs) {
        EstimatorStep(time, i);
        contact_->Apply(time);
        UpdateGhostNodesAndLinkToAdapter();
        CreateNodesJacobianAndLinkToAdapter();

        rc = CBSolverEquilibrium::SolverStep(time);

        UpdateGhostNodesAndLinkToAdapter();
        CreateNodesJacobianAndLinkToAdapter();

        VecZeroEntries(dist_);
        contact_->Apply(time);
        UpdateMasterNodesOfInterest();
        contact_->GetMasterNodesDistancesToSlaveElements(&dist_);
        VecGetSubVector(dist_, masterNodesIndices_, &subDist);
        TFloat lastNorm = norm;
        VecNorm(subDist, NORM_2, &norm);
        VecRestoreSubVector(dist_, masterNodesIndices_, &subDist);

        DCCtrl::print << norm << "\n";

        if (fabs(lastNorm - norm) < rel)
            break;
        else
            i++;

        UpdateGhostNodesAndLinkToAdapter();
        CreateNodesJacobianAndLinkToAdapter();
        WriteToFile(time, i, norm);
    }

    if (i == 0)
        WriteToFile(time, i, norm);

    contact_->Apply(time);
    UpdateMasterNodesOfInterest();
    contact_->GetMasterNodesDistancesToSlaveElements(&dist_);
    ExportNodesVectorData("Dist", dist_);

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    VecDestroy(&prevNodes);

    // Update the tension history (T_i, T_{i-1}, T_{i-2}) used for the temporal regularization terms.
    if (rc == CBStatus::FAILED) {
        PetscScalar *t;
        VecGetArray(ti1_, &t);
        for (auto &it : this->GetSolidElementVector())
            it->GetTensionModel()->SetActiveTensionAtQuadraturePoint(0, -t[it->GetLocalIndex()]);
        VecRestoreArray(ti1_, &t);
        VecCopy(ti1_, ti_);
    } else {
        VecCopy(ti1_, ti2_);
        PetscScalar *t;
        VecGetArray(ti1_, &t);
        for (auto &it : this->GetSolidElementVector())
            t[it->GetLocalIndex()] = -(it->GetTensionModel()->GetActiveTension());
        VecRestoreArray(ti1_, &t);
        VecCopy(ti1_, ti_);
    }

    return rc;
}

CBStatus CBSolverActiveStressEstimator::EstimatorStep(PetscScalar time, int step) {
    Vec subDist;
    UpdateGhostNodesAndLinkToAdapter();
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);

    // df/dtau: sensitivity of nodal forces to element active stress.
    Mat dfdtau;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 3 * numNodes_, GetNumberOfElements(), 0, dfdtauNnz_.data(), &dfdtau);
    adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
    formulation_->CalcNodalForcesActiveStressJacobian();
    MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
    Mat subdfdtau;
    MatCreateSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);

    // inv = (df/dx)^{-1} restricted to the nodes of interest, applied to the master-node unit loads.
    Mat inv;

    if (step == 0 || !(parameters_->Get<bool>("Solver.ActiveStressEstimator.ReUseInverse", false))) {
        CBSolver::CalcNodalForcesJacobian();
        MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);

        Mat subdfdx;
        MatCreateSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdx);
        IS perm, iperm;
        Mat f;
        MatFactorInfo info;
        MatFactorInfoInitialize(&info);
        MatGetOrdering(subdfdx, MATORDERINGND, &perm, &iperm);
        MatGetFactor(subdfdx, MATSOLVERMUMPS, MAT_FACTOR_LU, &f);
        MatLUFactorSymbolic(f, subdfdx, perm, iperm, &info);
        MatLUFactorNumeric(f, subdfdx, &info);

        // Solve for all master-node unit loads at once rather than one at a time: the right-hand
        // sides are the unit vectors at the master DOFs, which MUMPS takes in sparse form. PETSc
        // selects that path when the right-hand side is a virtual transpose of a sparse matrix.
        Mat bt;
        MatCreateSeqAIJ(PETSC_COMM_SELF, numMasterNodesIndices_, numNodesOfInterestIndices_, 1, nullptr, &bt);
        for (int i = 0; i < numMasterNodesIndices_; i++)
            MatSetValue(bt, i, masterNodesIndicesNodesOfInterestMapping_[i], 1, INSERT_VALUES);
        MatAssemblyBegin(bt, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(bt, MAT_FINAL_ASSEMBLY);

        Mat b, x;
        MatCreateTranspose(bt, &b);
        MatCreateSeqDense(PETSC_COMM_SELF, numNodesOfInterestIndices_, numMasterNodesIndices_, nullptr, &x);
        MatMatSolve(f, b, x);

        // The solutions come back as columns of x; inv holds them as rows.
        MatTranspose(x, MAT_INITIAL_MATRIX, &inv);

        MatDestroy(&x);
        MatDestroy(&b);
        MatDestroy(&bt);
        MatDestroy(&f);
        ISDestroy(&iperm);
        ISDestroy(&perm);
        MatDestroy(&subdfdx);

        if (inv_ == nullptr)
            MatDuplicate(inv, MAT_DO_NOT_COPY_VALUES, &inv_);
        MatCopy(inv, inv_, DIFFERENT_NONZERO_PATTERN);
        MatAssemblyBegin(inv_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(inv_, MAT_FINAL_ASSEMBLY);
    } else {
        MatDuplicate(inv_, MAT_COPY_VALUES, &inv);
    }

    // Sensitivity of the master-node positions to the element active stress, and normal equations.
    Mat a;
    MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);
    Mat aTa;
    MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);

    Vec dtau, d, dd, dd2;
    DCPetsc::CreateSeqVector(numElementOfInterestIndices_, &dtau);
    VecDuplicate(dtau, &d);
    VecDuplicate(dtau, &dd);
    VecDuplicate(dtau, &dd2);

    MatMultTranspose(a, subDist, d);

    // Tikhonov regularization terms.
    if (l1_ != 0) {
        VecSet(dd, l1_);
        MatDiagonalSet(aTa, dd, ADD_VALUES);
    }
    if (l2_ != 0) {
        VecSet(dd, l2_);
        MatDiagonalSet(aTa, dd, ADD_VALUES);
        VecAXPY(d, -l2_, ti_);
    }
    auto timeStep = timing_.GetTimeStep();
    if (l3_ != 0) {
        VecSet(dd, l3_ / timeStep);
        MatDiagonalSet(aTa, dd, ADD_VALUES);
        VecAXPY(d, -l3_ / timeStep, ti_);
        VecAXPY(d,  l3_ / timeStep, ti1_);
    }
    if (l4_ != 0) {
        VecSet(dd, l4_ / (timeStep*timeStep));
        MatDiagonalSet(aTa, dd, ADD_VALUES);
        VecAXPY(d,   -l4_ / (timeStep*timeStep), ti_);
        VecAXPY(d,  2*l4_ / (timeStep*timeStep), ti1_);
        VecAXPY(d,   -l4_ / (timeStep*timeStep), ti2_);
    }
    if (l5_ != 0) {
        MatAXPY(aTa, l5_, lTl_, DIFFERENT_NONZERO_PATTERN);
    }
    if (l6_ != 0) {
        MatAXPY(aTa, l6_, lTl_, DIFFERENT_NONZERO_PATTERN);
        MatMult(lTl_, ti_, dd);
        VecAXPY(d, -l6_, dd);
    }
    if (l7_ != 0) {
        MatAXPY(aTa, l7_, lTl_, DIFFERENT_NONZERO_PATTERN);
        VecCopy(ti_, dd2);
        VecAXPY(dd2, -1, ti1_);
        MatMult(lTl_, dd2, dd);
        VecAXPY(d, -l7_, dd);
    }

    KSP ksp;
    KSPCreate(DCPetsc::Comm(), &ksp);
    PC pc;
    KSPGetPC(ksp, &pc);
    KSPSetType(ksp, "preonly");
    KSPSetFromOptions(ksp);
    PCSetType(pc, PCLU);
    KSPSetOperators(ksp, aTa, aTa);
    KSPSetUp(ksp);
    KSPSolve(ksp, d, dtau);

    PetscScalar *t;
    PetscScalar *t2;
    VecGetArray(dtau, &t);
    VecGetArray(ti_, &t2);

    Vec currTension;
    VecDuplicate(dtau, &currTension);
    for (auto &it : this->GetSolidElementVector()) {
        int i = it->GetLocalIndex();
        TFloat valT = 0.0;
        TFloat val = it->GetTensionModel()->GetActiveTension();

        if (t[i] > 10000)
            t[i] = 10000;
        else if (t[i] < -10000)
            t[i] = -10000;

        if ((val - t[i]) < 0) {
            t[i] = val;
            valT = val - t[i];
        } else if ((val - t[i]) > 4e5) {
            t[i] = val - 4e5;
            valT = val - t[i];
        } else {
            valT = val - t[i];
        }

        t2[i] += t[i];
        it->GetTensionModel()->SetActiveTensionAtQuadraturePoint(0, valT);
        VecSetValue(currTension, i, valT, INSERT_VALUES);
    }
    VecAssemblyBegin(currTension);
    VecAssemblyEnd(currTension);

    // L-curve diagnostics: residual and regularization-term norms.
    Vec tmpVector, tmpVector1;
    VecDuplicate(subDist, &tmpVector);
    VecZeroEntries(tmpVector);
    MatMult(a, dtau, tmpVector);
    VecAXPY(tmpVector, -1, subDist);
    VecNorm(tmpVector, NORM_2, &resNorm_);

    VecDuplicate(dtau, &tmpVector1);
    VecZeroEntries(tmpVector1);
    MatMult(elementLaplacian_, currTension, tmpVector1);
    VecNorm(tmpVector1, NORM_2, &regTermNorm_);
    WriteLcurveStuffToFile(time, step, resNorm_, regTermNorm_);

    VecRestoreArray(dtau, &t);
    VecRestoreArray(ti_, &t2);
    VecRestoreSubVector(dist_, masterNodesIndices_, &subDist);

    VecDestroy(&tmpVector);
    VecDestroy(&tmpVector1);
    VecDestroy(&currTension);
    MatDestroy(&nodalForcesJacobian_);
    MatDestroy(&inv);
    MatDestroy(&a);
    MatDestroy(&aTa);
    MatDestroy(&dfdtau);
    MatDestroy(&subdfdtau);
    VecDestroy(&dtau);
    VecDestroy(&d);
    VecDestroy(&dd);
    VecDestroy(&dd2);
    KSPDestroy(&ksp);

    lastEstimatorTimeStep_ = time;
    lastTime_ = time;
    return CBStatus::SUCCESS;
}

void CBSolverActiveStressEstimator::GenerateElementLaplacian() {
    // Element-adjacency (common-face) Laplacian used for spatial smoothing regularization (l5..l7).
    Mat cMat;
    MatCreateAIJ(PETSC_COMM_WORLD, GetNumberOfSolidElements(), GetNumberOfSolidElements(),
                 PETSC_DETERMINE, PETSC_DETERMINE, 100, PETSC_NULLPTR, 100, PETSC_NULLPTR, &cMat);

    Vec diag;
    DCPetsc::CreateVector(GetNumberOfSolidElements(), PETSC_DETERMINE, &diag);
    VecSet(diag, 0);

    for (auto &e : model_->GetSolidElements()) {
        std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(e->GetIndex());
        for (auto &n : neighbors)
            MatSetValue(cMat, e->GetIndex(), n, -1, INSERT_VALUES);
    }

    MatDiagonalSet(cMat, diag, INSERT_VALUES);

    Vec b;
    VecDuplicate(diag, &b);
    VecSet(diag, 1);
    MatMult(cMat, diag, b);
    VecScale(b, -1);
    MatDiagonalSet(cMat, b, INSERT_VALUES);

    MatAssemblyBegin(cMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(cMat, MAT_FINAL_ASSEMBLY);

    Mat subCMat;
    MatCreateSubMatrix(cMat, elementsOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);

    MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &elementLaplacian_);
    MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);
    MatTransposeMatMult(elementLaplacian_, elementLaplacian_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &lTl_);

    MatDestroy(&cMat);
    MatDestroy(&subCMat);
    VecDestroy(&diag);
    VecDestroy(&b);
}

void CBSolverActiveStressEstimator::WriteHeaderToFile(std::string filename) {
    file_.open(filename.c_str());
    if (!file_.good())
        throw std::runtime_error("CBSolverActiveStressEstimator::WriteHeaderToFile: Couldn't create " + filename);
    file_ << std::setw(16) << "time" << std::setw(16) << "NewtonIteration" << std::setw(16) << "DistNorm" << std::endl;
    file_.close();

    fileLcurve_.open(filenameLcurve_.c_str());
    if (!fileLcurve_.good())
        throw std::runtime_error("CBSolverActiveStressEstimator::WriteHeaderToFile: Couldn't create " + filenameLcurve_);
    fileLcurve_ << std::setw(16) << "time" << std::setw(16) << "NewtonIter" << std::setw(16) << "ResidumNorm" << std::setw(16) << "RegTermNorm" << std::endl;
    fileLcurve_.close();
}

void CBSolverActiveStressEstimator::WriteToFile(PetscScalar time, int i, TFloat norm) {
    if (!headerWritten_)
        WriteHeaderToFile(filename_);
    headerWritten_ = true;

    file_.open(filename_.c_str(), std::ios::app);
    if (!file_.good())
        throw std::runtime_error("CBSolverActiveStressEstimator::WriteToFile: Couldn't create " + filename_);
    file_ << std::setprecision(7) << std::fixed << std::setw(16) << time << std::setw(16) << i << std::setw(16) << norm << std::endl;
    file_.close();
}

void CBSolverActiveStressEstimator::WriteLcurveStuffToFile(PetscScalar time, int i, TFloat resNorm, TFloat smoothNorm) {
    if (!headerLcurveWritten_)
        WriteHeaderToFile(filenameLcurve_);
    headerLcurveWritten_ = true;

    fileLcurve_.open(filenameLcurve_.c_str(), std::ios::app);
    if (!fileLcurve_.good())
        throw std::runtime_error("CBSolverActiveStressEstimator::WriteLcurveStuffToFile: Couldn't create " + filenameLcurve_);
    fileLcurve_ << std::setprecision(7) << std::fixed << std::setw(16) << time << std::setw(16) << i << std::setw(16) << resNorm << std::setw(16) << smoothNorm << std::endl;
    fileLcurve_.close();
}
