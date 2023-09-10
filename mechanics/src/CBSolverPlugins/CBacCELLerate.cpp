/* -------------------------------------------------------
 
 CBacCELLerate.cpp
 
 Ver. 1.2.0
 
 Created:       Robin Moss     (19.02.2015)
 Last modified: Tobias Gerach  (05.04.2023)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#include "filesystem.h"
#include "CBacCELLerate.h"
#include "CBDataCtrl.h"
#include "CBSolver.h"
#include "CBElementSolidT4.h"
#include <iostream>
#include <vtkTetra.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include "petscvec.h"
#include <petscviewer.h>
#include "vtksys/SystemTools.hxx"

#include <algorithm>
extern "C" double dgesvd_(const char *, const char *, int *, int *, double *, int *, double *, double *, int *,
                          double *, int *, double *, int *, int *);

std::vector<unsigned char> & split(const std::string &s, char delim, std::vector<unsigned char> &elems) {
    std::stringstream ss(s);
    std::string item;
    
    while (std::getline(ss, item, delim)) {
        elems.push_back((unsigned char)std::stoi(item));
    }
    return elems;
}

std::vector<unsigned char> split(const std::string &s, char delim) {
    std::vector<unsigned char> elems;
    
    split(s, delim, elems);
    return elems;
}

void CBacCELLerate::Init() {
    solidElements_ =  GetAdapter()->GetSolver()->GetSolidElementVector();
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_); CHKERRQ(ierr);
    DCCtrl::debug << "\n--- --- acCELLerate --- ---" << std::endl;
    
    InitParameters();
    
    /// Create accelerate instance and call acCELLerate constructor with the project file
    const char *pf = accprojectFile_.c_str();
    act_ = new acCELLerate();
    DCCtrl::debug << "Load project file ... ";
    InitPvdFile();
    act_->LoadProject(pf);
    DCCtrl::debug << "Done\n";
    
    DCCtrl::debug << "InitMono ... ";
    act_->InitMono();
    
    InitMesh();
    InitPetscVec();
    InitLocalMaps();
    InitMapping();
    
    /// copy sysMatrix and massMatrix allocation from acCELLerate instance
    sysMatrix_ = act_->GetSystemMatrix();
    massMatrix_ = act_->GetMassMatrix();
    DCCtrl::debug << "Done\n";
    
    if (MEF_ == "NONE") {
        /// skipping Prepare() phase
        status_ = CBStatus::DACCORD;
    }
    
    DCCtrl::debug << "\n--- --- acCELLerate --- ---\n" << std::endl;
} // CBacCELLerate::Init

void CBacCELLerate::Prepare() {
    DCCtrl::debug << "\n--- Preparing acCELLerate ---\n" << std::endl;
    
    /// Calculate shape function derivatives with respect to inflated state
    CalcShapeFunctionDeriv();
    
    /// Set new reference configuration
    UpdateNodes(true);
    
    /// not so nice....but vtk arrays have to be filled on process 0 to be correct and F_ was not computed before
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
    Vector3<TFloat> nodesCoords[4];
    for (int i = localElementsFrom_; i <= localElementsTo_; i++) {
        int localEleID = i - localElementsFrom_;
        Q_[localEleID] = GetInitialBasisAtCell(i);
        pointIds = acMesh_->GetCell(i)->GetPointIds();
        for (int j = 0; j < 4; j++) {
            nodesCoords[j] = acMesh_->GetPoint(pointIds->GetId(j));
            nodesCoords[j] /= 1000.0; // adjust mm -> m
        }
        CalcDeformationTensor(localEleID, nodesCoords, F_[localEleID]);
        Matrix3<TFloat> Q = GetBasisAtCell(localEleID);
        
        /// calculate deformed fibers "unloaded state"
        Vector3<TFloat> f = F_[localEleID] * Q.GetCol(0);
        Vector3<TFloat> s = F_[localEleID] * Q.GetCol(1);
        Vector3<TFloat> n = F_[localEleID] * Q.GetCol(2);
        
        /// do Gram Schmidt orthonormalization
        f.Normalize();
        s.Normalize();
        s = s - (f*s) * f;
        n.Normalize();
        n = n - (f*n) * f - (s*n) * s;
        
        /// set new reference basis
        Q_[localEleID] = {f(0), s(0), n(0), f(1), s(1), n(1), f(2), s(2), n(2)};
        
        /// carry F_ as PetscVec for export in P(0)
        PetscScalar F[9] =
        {F_[localEleID].Get(0), F_[localEleID].Get(1), F_[localEleID].Get(2), F_[localEleID].Get(3), F_[localEleID].Get(4),
            F_[localEleID].Get(5), F_[localEleID].Get(6), F_[localEleID].Get(7),
            F_[localEleID].Get(8)};
        PetscInt index[9] = {9*i+0, 9*i+1, 9*i+2, 9*i+3, 9*i+4, 9*i+5, 9*i+6, 9*i+7, 9*i+8};
        VecSetValues(deformation_, 9, index, F, INSERT_VALUES);
    }
    VecAssemblyBegin(deformation_);
    VecAssemblyEnd(deformation_);
    
    VecScatter ScatterDeformation;
    PetscErrorCode ierr;
    Vec localDeformation;
    
    ierr = VecScatterCreateToZero(deformation_, &ScatterDeformation, &localDeformation); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterDeformation, deformation_, localDeformation, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(
                                                                                                                        ierr);
    ierr = VecScatterEnd(ScatterDeformation, deformation_, localDeformation, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    
    PetscScalar *pD;
    ierr = VecGetArray(localDeformation, &pD); CHKERRQ(ierr);
    
    if (mpirank_ == 0) {
        for (int i = 0; i < nCells_; i++) {
            // find new basis
            Matrix3<TFloat> Q = GetInitialBasisAtCell(i);
            Matrix3<TFloat> F =
            {pD[9*i+0], pD[9*i+1], pD[9*i+2], pD[9*i+3], pD[9*i+4], pD[9*i+5], pD[9*i+6], pD[9*i+7], pD[9*i+8]};
            
            /// calculate deformed fibers
            Vector3<TFloat> f = F * Q.GetCol(0);
            Vector3<TFloat> s = F * Q.GetCol(1);
            Vector3<TFloat> n = F * Q.GetCol(2);
            
            /// do Gram Schmidt orthonormalization
            f.Normalize();
            s.Normalize();
            s = s - (f*s) * f;
            n.Normalize();
            n = n - (f*n) * f - (s*n) * s;
            
            acMeshFiberValues_->InsertTuple3(i, f(0), f(1), f(2));
            acMeshSheetValues_->InsertTuple3(i, s(0), s(1), s(2));
            acMeshNormalValues_->InsertTuple3(i, n(0), n(1), n(2));
        }
        
        ierr = VecRestoreArray(localDeformation, &pD); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ScatterDeformation); CHKERRQ(ierr);
        ierr = VecDestroy(&localDeformation); CHKERRQ(ierr);
    }
    
    /// Calculate shape function derivatives with respect to reference coordinates
    CalcShapeFunctionDeriv();
    
    /// Reset accNodes_ to current configuration
    UpdateNodes();
    
    DCCtrl::debug << "\n--- Preparing acCELLerate ---\n" << std::endl;
    
    status_ = CBStatus::DACCORD;
} // CBacCELLerate::Prepare

void CBacCELLerate::Apply(TFloat CMtime) {
    // Initialize all variables
    double t1 = MPI_Wtime();
    float time = CMtime - offsetTime_;
    DCCtrl::debug << "\n ------- acCELLerate -------\n";
    PetscErrorCode ierr;
    
    Vector3<TFloat> ffr;
    
    ffr = Vector3<TFloat>(1.0, 0.0, 0.0);
    
    if (time > stopTime_) {
        // Update Stretch/Velocity
        // Veloctiy is given in 1/s. Some ForceModels e.g. Land/Niederer17 need 1/ms. Check individually what is used in libcell
        ierr = VecCopy(stretchVecF_, velocityVec_); CHKERRQ(ierr);
        UpdateStretch();
        UpdateVelocity();
        
        // update coordinates of the acCELLerate mesh and recalculate the system matrix
        if ((MEF_ == "MINIMAL") || (MEF_ == "FULL")) {
            UpdateNodes();
            AssembleMatrix();
            act_->SetIntraMatrices(sysMatrix_, massMatrix_);
        }
        
        DCCtrl::debug << "Timestep : " << time << endl;
        stopTime_ = time;  // updateTime_;
        stepBackTime_ = time -  adapter_->GetSolver()->GetTiming().GetTimeStep();
        
        // Running acCELLerate
        DCCtrl::debug << "acCELLerate : " << stepBackTime_ << " - " << stopTime_ << "\n";
        
        std::ostringstream ss;
        ss << stopTime_;
        std::string s(ss.str());
        acltTime stopTime = s;
        act_->MonoDomain(stopTime, stretchVecF_, velocityVec_);
        DCCtrl::debug << "acCELLerate Done for : " << time << "\n";
        
        Vec ForceVec = act_->GetForceVec();
        ierr = VecCopy(timestepForce_, stepbackForce_); CHKERRQ(ierr);
        ierr = VecCopy(ForceVec, timestepForce_); CHKERRQ(ierr);
        VecScatter  ScatterForce;
        Vec localForce;
        
        ierr = VecScatterCreateToAll(ForceVec, &ScatterForce, &localForce); CHKERRQ(ierr);
        ierr = VecScatterBegin(ScatterForce, ForceVec, localForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(ScatterForce, ForceVec, localForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        
        PetscScalar *pif;
        ierr = VecGetArray(localForce, &pif); CHKERRQ(ierr);
        
        // Transfer calculated force to CM
        for (auto &e : solidElements_) {
            if (std::find(materialCoupling_.begin(), materialCoupling_.end(),
                          e->GetMaterialIndex()) != materialCoupling_.end()) {
                for (int QPi = 0; QPi < NumQP_; QPi++) {
                    Vector4<TFloat> ForceV4 = {pif[nearP_[e->GetIndex()*NumQP_ + QPi][0]],
                        pif[nearP_[e->GetIndex()*NumQP_ + QPi][1]],
                        pif[nearP_[e->GetIndex()*NumQP_ + QPi][2]],
                        pif[nearP_[e->GetIndex()*NumQP_ + QPi][3]]};
                    TFloat Force = ForceV4* ForceSF_[e->GetIndex()*NumQP_ + QPi];
                    
                    if (Force < 0) {
                        Force = 0; // negative Forces are rarely a nice thing to have
                    } else if (::isnan(Force) || ::isinf(Force) ) {
                        cout << "Force is NaN or inf  --> We might crash soon" << endl;
                    }
                    
                    e->GetTensionModel()->SetfibreRatio(ffr);
                    e->GetTensionModel()->SetActiveTensionAtQuadraturePoint(QPi, Force);
                }
            }
        }
        
        ierr = VecRestoreArray(localForce, &pif); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ScatterForce); CHKERRQ(ierr);
        ierr = VecDestroy(&localForce); CHKERRQ(ierr);
        
        DCCtrl::debug << "\n[" << MPI_Wtime() - t1 << " s]";
    } else if (stepBack_) {
        VecScatter  ScattertsForce, ScattersbForce;
        Vec localtsForce, localsbForce;
        
        ierr = VecScatterCreateToAll(timestepForce_, &ScattertsForce, &localtsForce); CHKERRQ(ierr);
        ierr = VecScatterBegin(ScattertsForce, timestepForce_, localtsForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(ScattertsForce, timestepForce_, localtsForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        
        ierr = VecScatterCreateToAll(stepbackForce_, &ScattersbForce, &localsbForce); CHKERRQ(ierr);
        ierr = VecScatterBegin(ScattersbForce, stepbackForce_, localsbForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(ScattersbForce, stepbackForce_, localsbForce, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        
        PetscScalar *ptsF, *psbF;
        ierr = VecGetArray(localtsForce, &ptsF); CHKERRQ(ierr);
        ierr = VecGetArray(localsbForce, &psbF); CHKERRQ(ierr);
        
        DCCtrl::debug << "\n-------\n Interpolating between timesteps \n-------\n" << endl;
        DCCtrl::debug << "Time: " << time << "\n";
        DCCtrl::debug << "StopTime: " << stopTime_ << "\n";
        DCCtrl::debug << "StepBackTime: " << stepBackTime_ << "\n";
        
        double stepBackFactor = 1 - ((time - stepBackTime_) / (stopTime_ - stepBackTime_));
        if (stepBackFactor == 0) {
            DCCtrl::debug << "StepBackFactor: " << stepBackFactor << "\n";
            DCCtrl::debug << "Doing nothing\n";
        } else {
            DCCtrl::debug << "StepBackFactor: " << stepBackFactor << "\n";
            
            for (auto &e : solidElements_) {
                if (find(materialCoupling_.begin(), materialCoupling_.end(),
                         e->GetMaterialIndex()) != materialCoupling_.end()) {
                    for (int QPi = 0; QPi < NumQP_; QPi++) {
                        Vector4<TFloat> ForceV4sF = {ptsF[nearP_[e->GetIndex()*NumQP_ + QPi][0]],
                            ptsF[nearP_[e->GetIndex()*NumQP_ + QPi][1]],
                            ptsF[nearP_[e->GetIndex()*NumQP_ + QPi][2]],
                            ptsF[nearP_[e->GetIndex()*NumQP_ + QPi][3]]};
                        TFloat timeStepForce = ForceV4sF* ForceSF_[e->GetIndex()*NumQP_ + QPi];
                        
                        Vector4<TFloat> ForceV4bF = {psbF[nearP_[e->GetIndex()*NumQP_ + QPi][0]],
                            psbF[nearP_[e->GetIndex()*NumQP_ + QPi][1]],
                            psbF[nearP_[e->GetIndex()*NumQP_ + QPi][2]],
                            psbF[nearP_[e->GetIndex()*NumQP_ + QPi][3]]};
                        TFloat stepBackForce = ForceV4bF* ForceSF_[e->GetIndex()*NumQP_ + QPi];
                        
                        TFloat Force = timeStepForce - (stepBackFactor * (timeStepForce - stepBackForce));
                        
                        if ((Force < 0) || (time < 0.005)) {
                            Force = 0;
                        } else if (::isnan(Force) || ::isinf(Force)) {
                            cout << "Force is NaN or Inf --> we might crash soon" << endl;
                            Force = 0;
                        }
                        e->GetTensionModel()->SetfibreRatio(ffr);
                        e->GetTensionModel()->SetActiveTensionAtQuadraturePoint(QPi, Force);
                    }
                }
            }
        }
        ierr = VecRestoreArray(localtsForce, &ptsF); CHKERRQ(ierr);
        ierr = VecRestoreArray(localsbForce, &psbF); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ScattersbForce); CHKERRQ(ierr);
        ierr = VecDestroy(&localsbForce); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ScattertsForce); CHKERRQ(ierr);
        ierr = VecDestroy(&localtsForce); CHKERRQ(ierr);
        stepBack_ = false;
    } else {
        DCCtrl::debug << "No Forces added by CBacCELLerate\n" << endl;
    }
    
    DCCtrl::debug << " CBacCELLerate Done" << endl;
    DCCtrl::debug << "\n-----------------------------------\n";
} // CBacCELLerate::Apply

void CBacCELLerate::Export(TFloat time) {
    if (export_ == true) {
        Vec potential;
        Vec calcium;
        
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &potential);
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &calcium);
        
        PetscInt from1, to1;
        PetscInt from2, to2;
        VecGetOwnershipRange(potential, &from1, &to1);
        VecGetOwnershipRange(calcium, &from2, &to2);
        
        VecZeroEntries(potential);
        VecZeroEntries(calcium);
        
        Vec PotentialVec = act_->GetVmVec();
        Vec CalciumVec = act_->GetForceVec();
        
        VecScatter  ScatterCalcium, ScatterPotential;
        Vec localCalcium, localPotential;
        
        VecScatterCreateToAll(CalciumVec, &ScatterCalcium, &localCalcium);
        VecScatterBegin(ScatterCalcium, CalciumVec, localCalcium, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ScatterCalcium, CalciumVec, localCalcium, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterCreateToAll(PotentialVec, &ScatterPotential, &localPotential);
        VecScatterBegin(ScatterPotential, PotentialVec, localPotential, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ScatterPotential, PotentialVec, localPotential, INSERT_VALUES, SCATTER_FORWARD);
        
        PetscScalar *piC, *piV;
        VecGetArray(localCalcium, &piC);
        VecGetArray(localPotential, &piV);
        
        for (auto &e : solidElements_) {
            if (std::find(materialCoupling_.begin(), materialCoupling_.end(),
                          e->GetMaterialIndex()) != materialCoupling_.end()) {
                for (int QPi = 0; QPi < NumQP_; QPi++) {
                    Vector4<TFloat> CalciumV4 = {piC[nearP_[e->GetIndex()*NumQP_ + QPi][0]],
                        piC[nearP_[e->GetIndex()*NumQP_ + QPi][1]],
                        piC[nearP_[e->GetIndex()*NumQP_ + QPi][2]],
                        piC[nearP_[e->GetIndex()*NumQP_ + QPi][3]]};
                    TFloat Calcium = CalciumV4 * ForceSF_[e->GetIndex()*NumQP_ + QPi];
                    
                    Vector4<TFloat> PotentialV4 = {piV[nearP_[e->GetIndex()*NumQP_ + QPi][0]],
                        piV[nearP_[e->GetIndex()*NumQP_ + QPi][1]],
                        piV[nearP_[e->GetIndex()*NumQP_ + QPi][2]],
                        piV[nearP_[e->GetIndex()*NumQP_ + QPi][3]]};
                    TFloat Potential = PotentialV4 * ForceSF_[e->GetIndex()*NumQP_ + QPi] * 1000;
                    
                    VecSetValue(calcium, from2 + e->GetLocalIndex(), Calcium, INSERT_VALUES);
                    VecSetValue(potential, from1 + e->GetLocalIndex(), Potential, INSERT_VALUES);
                }
            }
        }
        
        GetAdapter()->GetSolver()->ExportElementsScalarData("Potential", potential);
        GetAdapter()->GetSolver()->ExportElementsScalarData("Calcium", calcium);
        
        VecRestoreArray(localCalcium, &piC);
        VecRestoreArray(localPotential, &piV);
        
        VecScatterDestroy(&ScatterCalcium);
        VecScatterDestroy(&ScatterPotential);
        
        VecDestroy(&localPotential);
        VecDestroy(&localCalcium);
        VecDestroy(&calcium);
        VecDestroy(&potential);
    }
} // CBacCELLerate::Export

void CBacCELLerate::StepBack() {
    stepBack_ = true;
}

// ----------- Private ------------

void CBacCELLerate::UpdateStretch() {
    PetscErrorCode ierr;
    DCCtrl::debug << "Updating stretch ...";
    Matrix3<TFloat> f[4];
    double alpha = (5 + 3 * sqrt(5)) / 20;
    double beta = (5 - sqrt(5)) / 20;
    Matrix4<TFloat> QPInv = Matrix4<TFloat>(alpha, beta, beta, beta,
                                            beta, alpha, beta, beta,
                                            beta, beta, alpha, beta,
                                            beta, beta, beta, alpha).GetInverse();
    
    ierr = VecSet(stretchVecF_, 0); CHKERRQ(ierr);
    
    for (auto it = eleList_.begin(); it != eleList_.end(); it++) {
        it->second->GetDeformationTensorAtQuadraturePoints(f);
        
        Vector4<TFloat> LocalPos = Vector4<TFloat>(nearC_[it->first][1],
                                                   nearC_[it->first][2],
                                                   nearC_[it->first][3],
                                                   nearC_[it->first][4]);
        Vector4<TFloat> StretchAtQP = Vector4<TFloat>(sqrt(f[0].GetCol(0)*f[0].GetCol(0)),
                                                      sqrt(f[1].GetCol(0)*f[1].GetCol(0)),
                                                      sqrt(f[2].GetCol(0)*f[2].GetCol(0)),
                                                      sqrt(f[3].GetCol(0)*f[3].GetCol(0)));
        
        TFloat StretchVal = LocalPos * (QPInv * StretchAtQP);
        ierr = VecSetValue(stretchVecF_, it->first, StretchVal, INSERT_VALUES);
    }
    ierr = VecAssemblyBegin(stretchVecF_); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(stretchVecF_); CHKERRQ(ierr);
    
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::UpdateStretch

void CBacCELLerate::UpdateVelocity() {
    PetscErrorCode ierr;
    
    DCCtrl::debug << "Updating velocity ...";
    
    if (constStretchRate_) {
        ierr = VecSet(velocityVec_, 0); CHKERRQ(ierr);
    } else {
        ierr = VecAYPX(velocityVec_, -1, stretchVecF_); CHKERRQ(ierr);
        ierr = VecScale(velocityVec_, 1/adapter_->GetSolver()->GetTiming().GetTimeStep()); CHKERRQ(ierr);
    }
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::UpdateVelocity

void CBacCELLerate::UpdateNodes(bool useReferenceNodes) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    PetscErrorCode ierr;
    
    ierr = VecSet(accNodes_, 0); CHKERRQ(ierr);
    Vector4<TFloat> l;
    
    if (useReferenceNodes)
        DCCtrl::debug << "Resetting Nodes to reference configuration ... ";
    else
        DCCtrl::debug << "Updating Nodes ... ";
    
    for (auto it = eleList_.begin(); it != eleList_.end(); it++) {
        CBElementSolid *e = it->second;
        if (e == 0) {
            throw std::runtime_error(
                                     "An element of the list connecting the nodes to the corresponding CM elements seems to be empty");
        }
        PetscInt indices[12];
        PetscScalar coords[12];
        for (unsigned int i = 0; i < 4; i++) {
            indices[3*i]   = 3*(e)->GetNodeIndex(i);
            indices[3*i+1] = 3*(e)->GetNodeIndex(i)+1;
            indices[3*i+2] = 3*(e)->GetNodeIndex(i)+2;
        }
        
        if (useReferenceNodes)
            GetAdapter()->GetRefNodesCoords(12, indices, coords);
        else
            GetAdapter()->GetNodesCoords(12, indices, coords);
        
        Vector3<TFloat> p1(&coords[0]);
        Vector3<TFloat> p2(&coords[3]);
        Vector3<TFloat> p3(&coords[6]);
        Vector3<TFloat> p4(&coords[9]);
        
        Vector3<TFloat> p;
        p(0) = (p1.X() * nearC_[it->first][1] +
                p2.X() * nearC_[it->first][2] +
                p3.X() * nearC_[it->first][3] +
                p4.X() * nearC_[it->first][4]) * 1000;
        
        p(1) = (p1.Y() * nearC_[it->first][1] +
                p2.Y() * nearC_[it->first][2] +
                p3.Y() * nearC_[it->first][3] +
                p4.Y() * nearC_[it->first][4]) * 1000;
        
        p(2) = (p1.Z() * nearC_[it->first][1] +
                p2.Z() * nearC_[it->first][2] +
                p3.Z() * nearC_[it->first][3] +
                p4.Z() * nearC_[it->first][4]) * 1000;
        
        
        ierr = VecSetValue(accNodes_, it->first * 3 + 0, p(0), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(accNodes_, it->first * 3 + 1, p(1), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(accNodes_, it->first * 3 + 2, p(2), INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecAssemblyBegin(accNodes_); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(accNodes_); CHKERRQ(ierr);
    
    
    VecScatter  ScatterElePoints;
    Vec localElePoints;
    
    ierr = VecScatterCreateToAll(accNodes_, &ScatterElePoints, &localElePoints); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterElePoints, accNodes_, localElePoints, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ScatterElePoints, accNodes_, localElePoints, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    PetscScalar *plN;
    ierr = VecGetArray(localElePoints, &plN); CHKERRQ(ierr);
    
    for (int i = 0; i < nPoints_; i++) {
        points->InsertNextPoint(plN[i*3+0], plN[i*3+1], plN[i*3+2]);
    }
    
    ierr = VecRestoreArray(localElePoints, &plN); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ScatterElePoints); CHKERRQ(ierr);
    ierr = VecDestroy(&localElePoints); CHKERRQ(ierr);
    acMesh_->SetPoints(points);
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::UpdateNodes

void CBacCELLerate::CalcShapeFunctionDeriv() {
    /// this function is only to be called during Prepare()
    vtkSmartPointer<vtkIdList> Points = vtkSmartPointer<vtkIdList>::New();
    Vector3<TFloat> nodesCoords[4];
    
    DCCtrl::debug << "Updating shape function derivatives ... ";
    
    for (int CellId = localElementsFrom_; CellId <= localElementsTo_; CellId++) {
        Points = acMesh_->GetCell(CellId)->GetPointIds();
        for (int i = 0; i < 4; i++) {
            nodesCoords[i] = acMesh_->GetPoint(Points->GetId(i));
            nodesCoords[i] /= 1000.0; // adjust mm -> m
        }
        
        TFloat z43 = (nodesCoords[3].Z() - nodesCoords[2].Z());
        TFloat z42 = (nodesCoords[3].Z() - nodesCoords[1].Z());
        TFloat z41 = (nodesCoords[3].Z() - nodesCoords[0].Z());
        TFloat z32 = (nodesCoords[2].Z() - nodesCoords[1].Z());
        TFloat z31 = (nodesCoords[2].Z() - nodesCoords[0].Z());
        TFloat z21 = (nodesCoords[1].Z() - nodesCoords[0].Z());
        
        TFloat y43 = (nodesCoords[3].Y() - nodesCoords[2].Y());
        TFloat y42 = (nodesCoords[3].Y() - nodesCoords[1].Y());
        TFloat y41 = (nodesCoords[3].Y() - nodesCoords[0].Y());
        TFloat y32 = (nodesCoords[2].Y() - nodesCoords[1].Y());
        TFloat y31 = (nodesCoords[2].Y() - nodesCoords[0].Y());
        TFloat y21 = (nodesCoords[1].Y() - nodesCoords[0].Y());
        
        TFloat x41 = (nodesCoords[3].X() - nodesCoords[0].X());
        TFloat x31 = (nodesCoords[2].X() - nodesCoords[0].X());
        TFloat x21 = (nodesCoords[1].X()- nodesCoords[0].X());
        
        TFloat detJ = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
        
        int localEleID = CellId - localElementsFrom_;
        
        dNdX_[localEleID][0] = 1.0 / detJ *(nodesCoords[1].Y() * z43 - nodesCoords[2].Y() * z42 + nodesCoords[3].Y() * z32);
        dNdX_[localEleID][3] = 1.0 / detJ *(-nodesCoords[0].Y()* z43 + nodesCoords[2].Y() * z41 - nodesCoords[3].Y() * z31);
        dNdX_[localEleID][6] = 1.0 / detJ *(nodesCoords[0].Y()* z42 - nodesCoords[1].Y() * z41 + nodesCoords[3].Y() * z21);
        dNdX_[localEleID][9] = 1.0 / detJ *(-nodesCoords[0].Y()* z32 + nodesCoords[1].Y() * z31 - nodesCoords[2].Y() * z21);
        
        dNdX_[localEleID][1] = 1.0 / detJ *(-nodesCoords[1].X()* z43 + nodesCoords[2].X() * z42 - nodesCoords[3].X() * z32);
        dNdX_[localEleID][4] = 1.0 / detJ *(nodesCoords[0].X()* z43 - nodesCoords[2].X() * z41 + nodesCoords[3].X() * z31);
        dNdX_[localEleID][7] = 1.0 / detJ *(-nodesCoords[0].X()* z42 + nodesCoords[1].X()* z41 - nodesCoords[3].X() * z21);
        dNdX_[localEleID][10] = 1.0 / detJ *(nodesCoords[0].X()* z32 - nodesCoords[1].X()* z31 + nodesCoords[2].X() * z21);
        
        dNdX_[localEleID][2] = 1.0 / detJ *(nodesCoords[1].X()* y43 - nodesCoords[2].X() * y42 + nodesCoords[3].X() * y32);
        dNdX_[localEleID][5] = 1.0 / detJ *(-nodesCoords[0].X()* y43 + nodesCoords[2].X() * y41 - nodesCoords[3].X() * y31);
        dNdX_[localEleID][8] = 1.0 / detJ *(nodesCoords[0].X()* y42 - nodesCoords[1].X()* y41 + nodesCoords[3].X() * y21);
        dNdX_[localEleID][11] = 1.0 / detJ *(-nodesCoords[0].X()* y32 + nodesCoords[1].X()* y31 - nodesCoords[2].X() * y21);
        
        // volume of cell
        dNdX_[localEleID][12] = detJ / 6.0;
    }
    
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::CalcShapeFunctionDeriv

void CBacCELLerate::CalcDeformationTensor(vtkIdType cellID, const Vector3<TFloat> *nodesCoords,
                                          Matrix3<TFloat> &deformationTensor) {
    TFloat *f = deformationTensor.GetArray();
    
    for (unsigned int i = 0; i < 3; i++) {
        f[i]   = dNdX_[cellID][i] * nodesCoords[0].X() + dNdX_[cellID][i+3] * nodesCoords[1].X() + dNdX_[cellID][i+6] *
        nodesCoords[2].X() + dNdX_[cellID][i+9] *
        nodesCoords[3].X();
        f[3+i] = dNdX_[cellID][i] * nodesCoords[0].Y() + dNdX_[cellID][i+3] * nodesCoords[1].Y() + dNdX_[cellID][i+6] *
        nodesCoords[2].Y() + dNdX_[cellID][i+9] *
        nodesCoords[3].Y();
        f[6+i] = dNdX_[cellID][i] * nodesCoords[0].Z() + dNdX_[cellID][i+3] * nodesCoords[1].Z() + dNdX_[cellID][i+6] *
        nodesCoords[2].Z() + dNdX_[cellID][i+9] *
        nodesCoords[3].Z();
    }
    deformationTensor = GetBasisAtCell(cellID).GetTranspose() * deformationTensor *
    GetBasisAtCell(cellID).GetInverse().GetTranspose();
}

void CBacCELLerate::AssembleMatrix() {
    PetscErrorCode ierr;
    float AnisotropyX[256];
    float AnisotropyY[256];
    float AnisotropyZ[256];
    int materialPriority[256];
    float K[256];
    static const double Frequency = 0.0;
    
    /// list of active materials in acMesh_
    MaterialListe materialProperties(MaterialFileName_.c_str());
    std::vector<int> nodeMaterial(nPoints_, 0);
    
    /// prioritize smaller anatomical structures over myocardium
    std::vector<unsigned char> priorityVector;
    
    priorityVector = split("72,73,74,75,76,77,78,29,30,31,34,35,32,33,79,80", ',');
    
    DCCtrl::debug << "Assemble matrices ... ";
    
    /// find conductivities in each material
    for (int mat = 0; mat < 256; ++mat) {
        materialPriority[mat] = 0;
        Material *material = materialProperties.Suchen(mat);
        if (material) {
            AnisotropyX[mat] = material->HoleAnisotropyX();
            AnisotropyY[mat] = material->HoleAnisotropyY();
            AnisotropyZ[mat] = material->HoleAnisotropyZ();
            K[mat] = material->LFkappa::Hole(Frequency);
        } else {
            AnisotropyX[mat] = AnisotropyY[mat] = AnisotropyZ[mat] = K[mat] = -1.0;
        }
    }
    
    /// order material by priority
    int prio = 0;
    for (std::vector<unsigned char>::reverse_iterator i = priorityVector.rbegin(); i != priorityVector.rend(); ++i)
        materialPriority[*i] = ++prio;
    
    /// allocate system matrix or set to 0 if it already exists
    if (sysMatrix_) {
        ierr = MatZeroEntries(sysMatrix_);
        CHKERRQ(ierr);
    } else {
        DCCtrl::debug << "Allocate system matrix ... ";
        ierr = MatCreate(PETSC_COMM_WORLD, &sysMatrix_); CHKERRQ(ierr);
        ierr = MatSetType(sysMatrix_, MATMPIAIJ); CHKERRQ(ierr);
        ierr = MatSetSizes(sysMatrix_, PETSC_DECIDE, PETSC_DECIDE, int(nPoints_), int(nPoints_)); CHKERRQ(ierr);
        ierr = MatSetFromOptions(sysMatrix_); CHKERRQ(ierr);
        ierr = MatSetUp(sysMatrix_); CHKERRQ(ierr);
        ierr = MatGetOwnershipRange(sysMatrix_, &Istart, &Iend); CHKERRQ(ierr);
        DCCtrl::debug << "Done\n";
    }
    
    /// allocate mass matrix or set to 0 if it already exists
    if (massMatrix_) {
        ierr = MatZeroEntries(massMatrix_);
        CHKERRQ(ierr);
    } else {
        DCCtrl::debug << "Allocate mass matrix ... ";
        ierr = MatCreate(PETSC_COMM_WORLD, &massMatrix_); CHKERRQ(ierr);
        ierr = MatSetType(massMatrix_, MATMPIAIJ); CHKERRQ(ierr);
        ierr = MatSetSizes(massMatrix_, PETSC_DECIDE, PETSC_DECIDE, int(nPoints_), int(nPoints_)); CHKERRQ(ierr);
        ierr = MatSetFromOptions(massMatrix_); CHKERRQ(ierr);
        ierr = MatSetUp(massMatrix_); CHKERRQ(ierr);
        DCCtrl::debug << "Done\n";
    }
    
    /// loop over cells
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
    Vector3<TFloat> nodesCoords[4];
    for (vtkIdType cellId = localElementsFrom_; cellId <= localElementsTo_; cellId++) {
        int localEleID = cellId - localElementsFrom_;
        int material = 1;
        if (acMeshMaterials_)
            material = acMeshMaterials_->GetValue(cellId);
        if (abs(K[material] - (-1.0)) < 1e-6) {
            throw std::runtime_error("\n\nMaterial " + std::to_string(
                                                                      material) + " does not exist in the material file " + MaterialFileName_ + ".");
        }
        
        /// loop over vertices
        pointIds = acMesh_->GetCell(cellId)->GetPointIds();
        TInt p[4];
        for (int i = 0; i < 4; i++) {
            nodesCoords[i] = acMesh_->GetPoint(pointIds->GetId(i));
            nodesCoords[i] /= 1000.0; // adjust mm -> m
            p[i] = int(pointIds->GetId(i)); // pointIds as int for matrix assembly
                                            /// assign correct material
            if (materialPriority[material] >= materialPriority[nodeMaterial[pointIds->GetId(i)]])
                nodeMaterial[pointIds->GetId(i)] = material;
        } // end loop over vertices
        
        /// calc diffusion tensor
        Matrix3<TFloat> Q = GetBasisAtCell(localEleID);
        Matrix3<TFloat> D = { K[material] * AnisotropyX[material], 0, 0,
            0, K[material] * AnisotropyY[material], 0,
            0, 0, K[material] * AnisotropyZ[material] };
        
        CalcDeformationTensor(localEleID, nodesCoords, F_[localEleID]);
        double J = F_[localEleID].Det();
        
        if (J <= 0) {
            DCCtrl::debug << "\nCorrupt element ID: " << cellId << "\n";
            J = 1;
            F_[localEleID] = Matrix3<TFloat>::Identity();
        }
        
        PetscScalar F[9] =
        {F_[localEleID].Get(0), F_[localEleID].Get(1), F_[localEleID].Get(2), F_[localEleID].Get(3), F_[localEleID].Get(4),
            F_[localEleID].Get(5),
            F_[localEleID].Get(
                               6), F_[localEleID].Get(7), F_[localEleID].Get(8)};
        PetscInt index[9] =
        {9*cellId+0, 9*cellId+1, 9*cellId+2, 9*cellId+3, 9*cellId+4, 9*cellId+5, 9*cellId+6, 9*cellId+7, 9*cellId+8};
        VecSetValues(deformation_, 9, index, F, INSERT_VALUES);
        
        if (MEF_ == "MINIMAL") {
            // no MEF on D
            D = Q * D * Q.GetTranspose();
            D = J * F_[localEleID].GetInverse() * D * F_[localEleID].GetInverse().GetTranspose();
        } else if (MEF_ == "FULL") {
            /// calculate deformed fibers
            Vector3<TFloat> f = F_[localEleID] * Q.GetCol(0);
            Vector3<TFloat> s = F_[localEleID] * Q.GetCol(1);
            Vector3<TFloat> n = F_[localEleID] * Q.GetCol(2);
            
            D =
            D(0,
              0) *
            DyadicProduct(f,
                          f)/(f.Norm()*f.Norm()) +
            D(1, 1) * DyadicProduct(s, s)/(s.Norm()*s.Norm()) + D(2, 2) * DyadicProduct(
                                                                                        n, n)/(n.Norm()*n.Norm());
            
            // rotate
            D = J * F_[localEleID].GetInverse() * D * F_[localEleID].GetInverse().GetTranspose();
        }
        
        /// assemble mass (M) and stiffness (K) matrix
        /// rebuilding M is theoretically not needed, but the one we get from acCELLerate is already scaled and would require changing a lot of code to work
        double M[16], K[16];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Vector3<double> dNdXi(&dNdX_[localEleID][3*i]);
                Vector3<double> dNdXj(&dNdX_[localEleID][3*j]);
                double vol = dNdX_[localEleID][12];
                
                M[4*i+j] = vol * (i == j ? 1.0 / 10.0 : 1.0 / 20.0);
                K[4*i+j] = vol * D * dNdXi * dNdXj;
            }
        } // end assembly loop
        
        ierr = MatSetValues(sysMatrix_, 4, p, 4, p, K, ADD_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(massMatrix_, 4, p, 4, p, M, ADD_VALUES); CHKERRQ(ierr);
    } // end loop over cells
    
    ierr = VecAssemblyBegin(deformation_);
    ierr = VecAssemblyEnd(deformation_);
    ierr = MatAssemblyBegin(sysMatrix_, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(sysMatrix_, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(massMatrix_, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(massMatrix_, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::AssembleMatrix

void CBacCELLerate::InitParameters() {
    DCCtrl::debug << "\nLoading settings ...";
    
    NumQP_ = solidElements_[0]->GetNumberOfQuadraturePoints();
    DCCtrl::debug << "\n Number of quad. points: " << NumQP_;
    
    accprojectFile_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.ProjectFile");
    DCCtrl::debug << "\n Project File: " << accprojectFile_;
    
    MaterialFileName_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.MaterialFile");
    DCCtrl::debug << "\n Material File: " << MaterialFileName_;
    
    resPreFix_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.ResultPreFix");
    DCCtrl::debug << "\n Result Prefix: " << resPreFix_;
    
    resFolder_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.ResultFolder");
    DCCtrl::debug << "\n Result Folder: " << resFolder_;
    
    export_ = GetParameters()->Get<bool>("Plugins.acCELLerate.Export", false);
    DCCtrl::debug << "\n Export: " << export_;
    
    MEF_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.MEF");
    DCCtrl::debug << "\n MEF type: " << MEF_;
    
    pvdFilename_ = GetParameters()->Get<std::string>("Plugins.acCELLerate.PvdFileName");
    DCCtrl::debug << "\n PVD Filename: " << pvdFilename_;
    
    materialCoupling_ =  parameters_->GetArray<TInt>("Plugins.acCELLerate.Material");
    DCCtrl::debug << "\n Active materials: ";
    for (std::vector<TInt>::const_iterator i = materialCoupling_.begin(); i != materialCoupling_.end(); ++i)
        DCCtrl::debug << *i << ", ";
    
    offsetTime_ =   GetParameters()->Get<float>("Plugins.acCELLerate.OffsetTime", 0);
    DCCtrl::debug << "\n Offset Time: " << offsetTime_;
    
    constStretchRate_ = GetParameters()->Get<bool>("Plugins.acCELLerate.constStretchRate", false);
    DCCtrl::debug << "\n Constant stretch-rate: " << constStretchRate_;
    
    if (!(MEF_ == "NONE") && !(MEF_ == "MINIMAL") && !(MEF_ == "FULL"))
        throw std::runtime_error("Plugins.acCELLerate.MEF '" + MEF_ + "' not supported. (NONE, MINIMAL, FULL)");
    
    DCCtrl::debug << "\nLoad settings ... Done";
} // CBacCELLerate::InitParameters

void CBacCELLerate::InitMesh() {
    DCCtrl::debug << "\nLoading mesh for electrophysiology ... ";
    std::string acMeshFilename = GetParameters()->Get<std::string>("Plugins.acCELLerate.acCELLerateMesh");
    
    std::string extension = vtksys::SystemTools::GetFilenameLastExtension(acMeshFilename);
    
    if (extension != ".vtu") {
        throw std::runtime_error(
                                 "Plugins.acCELLerate.acCELLerateMesh is of the wrong type. Please provide the file format .vtu.");
    }
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(acMeshFilename.c_str());
    reader->Update();
    acMesh_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    acMesh_ = reader->GetOutput();
    nPoints_ = acMesh_->GetNumberOfPoints();
    nCells_ = acMesh_->GetNumberOfCells();
    DCCtrl::debug << "\nACMesh Properties: \nPoints: " << nPoints_ << "\nCells: " << nCells_;
    
    if (GetParameters()->Get<bool>("Plugins.acCELLerate.Permute", false)) {
        DCCtrl::debug << "\nPermuting Acc-Mesh ... ";
        ApplySpatialSortPCA();
    }
    
    DCCtrl::debug << "\nLoading fiber orientation ... ";
    acMeshMaterials_ = vtkDoubleArray::SafeDownCast(((acMesh_->GetCellData()->GetArray("Material"))));
    acMeshFiberValues_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Fiber")));
    acMeshSheetValues_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Sheet")));
    acMeshNormalValues_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Sheetnormal")));
    if (!acMeshFiberValues_ || !acMeshSheetValues_ || !acMeshNormalValues_) {
        DCCtrl::debug << "Info: No fiber orientation defined.\n";
        throw std::runtime_error("Automatic mapping of fibers currently not supported.");
        
        //    DCCtrl::debug << "Getting Fibres from CM Mesh ... ";
        //    GetFibersFromMechanicsMesh();
        //    acMeshFiberValues_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Fiber")));
    }
    DCCtrl::debug << "Done\n";
    
    if (!acMeshMaterials_) {
        throw std::runtime_error("CBacCELLerate::InitMesh(): No Material defined.");
    }
    
    DetermineElementRanges();
} // CBacCELLerate::InitMesh

void CBacCELLerate::InitMapping() {
    DCCtrl::debug << "Creating mesh connections ... \n";
    
    /// shape functions for gauss points of tetrahedron
    std::vector<double> ShapeFunVec(20, 0);
    double alpha  = (5 + 3 * sqrt(5)) / 20;
    double beta   = (5 - sqrt(5)) / 20;
    
    ShapeFunVec   = { 0.25, 0.25, 0.25, 0.25,
        alpha, beta, beta, beta,
        beta, alpha, beta, beta,
        beta, beta, alpha, beta,
        beta, beta, beta, alpha };
    
    PetscInt indices[12];
    PetscScalar coords[12];
    for (int i = 0; i < nPoints_; i++) {
        nearC_[i] = {INFINITY, 0, 0, 0, 0, 0};
    }
    int skippedElement = 0;
    vtkSmartPointer<vtkPoints> CenterPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> TempVTK = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkIdList> CellPoints = vtkSmartPointer<vtkIdList>::New();
    for (int CellI = 0; CellI < nCells_; CellI++) {
        CellPoints = acMesh_->GetCell(CellI)->GetPointIds();
        double Center[3] = {(acMesh_->GetPoint(CellPoints->GetId(0))[0] +
                             acMesh_->GetPoint(CellPoints->GetId(1))[0] +
                             acMesh_->GetPoint(CellPoints->GetId(2))[0] +
                             acMesh_->GetPoint(CellPoints->GetId(3))[0])/4,
            (acMesh_->GetPoint(CellPoints->GetId(0))[1] +
             acMesh_->GetPoint(CellPoints->GetId(1))[1] +
             acMesh_->GetPoint(CellPoints->GetId(2))[1] +
             acMesh_->GetPoint(CellPoints->GetId(3))[1])/4,
            (acMesh_->GetPoint(CellPoints->GetId(0))[2] +
             acMesh_->GetPoint(CellPoints->GetId(1))[2] +
             acMesh_->GetPoint(CellPoints->GetId(2))[2] +
             acMesh_->GetPoint(CellPoints->GetId(3))[2])/4};
        CenterPoints->InsertNextPoint(Center);
    }
    TempVTK->SetPoints(CenterPoints);
    
    /// PointLocator contains centroids of acMesh_.Cells
    vtkSmartPointer<vtkPointLocator> PointLocator = vtkSmartPointer<vtkPointLocator>::New();
    PointLocator->SetDataSet(TempVTK);
    PointLocator->BuildLocator();
    
    
    PetscErrorCode ierr;
    
    // For each CM element find the corresponding EP Points and vice versa
    for (auto &e : solidElements_) {
        if (std::find(materialCoupling_.begin(), materialCoupling_.end(),
                      e->GetMaterialIndex()) != materialCoupling_.end()) {
            indices[0] = 3 * (e)->GetNodeIndex(0);
            indices[1] = 3 * (e)->GetNodeIndex(0) + 1;
            indices[2] = 3 * (e)->GetNodeIndex(0) + 2;
            
            indices[3] = 3 * (e)->GetNodeIndex(1);
            indices[4] = 3 * (e)->GetNodeIndex(1) + 1;
            indices[5] = 3 * (e)->GetNodeIndex(1) + 2;
            
            indices[6] = 3 * (e)->GetNodeIndex(2);
            indices[7] = 3 * (e)->GetNodeIndex(2) + 1;
            indices[8] = 3 * (e)->GetNodeIndex(2) + 2;
            
            indices[9] = 3 * (e)->GetNodeIndex(3);
            indices[10] = 3 * (e)->GetNodeIndex(3) + 1;
            indices[11] = 3 * (e)->GetNodeIndex(3) + 2;
            
            GetAdapter()->GetNodesCoords(12, indices, coords);
            
            /// vertices of solid element e
            Vector3<TFloat> p1(&coords[0]);
            Vector3<TFloat> p2(&coords[3]);
            Vector3<TFloat> p3(&coords[6]);
            Vector3<TFloat> p4(&coords[9]);
            
            /// Jacobian matrix of linear tetrahedron
            Matrix4<TFloat> m = {p1.X(), p2.X(), p3.X(), p4.X(),
                p1.Y(), p2.Y(), p3.Y(), p4.Y(),
                p1.Z(), p2.Z(), p3.Z(), p4.Z(),
                1,      1,      1,      1};
            m.Invert();
            
            for (PetscInt i = 0; i < nPoints_; i++) {
                double p[3] = {acMesh_->GetPoint(i)[0], acMesh_->GetPoint(i)[1], acMesh_->GetPoint(i)[2]};
                
                /// l: point p expressed with shape functions of solidElement e
                Vector4<TFloat> l = m * Vector4<TFloat>(p[0] / 1000, p[1] / 1000, p[2] / 1000, 1);
                
                /// since the centroid of the linear tetrahedron is expressed as l={0.25,0.25,0.25,0.25}, distance to the centroid becomes minimal if l.max() - l.min() approaches 0
                double max = l.Max();
                double min = l.Min();
                double CDist = max - min;
                
                /// if point p is inside solidElement e, update map
                if ((min >= 0) && (max <= 1)) {
                    nearC_[i][0] = CDist;
                    nearC_[i][1] = l(0);
                    nearC_[i][2] = l(1);
                    nearC_[i][3] = l(2);
                    nearC_[i][4] = l(3);
                    nearC_[i][5] = 1; // this value marks the point as mapped
                    
                    /// map current acc node i to solid element e
                    eleList_[i] = (e);
                    
                    ierr =
                    VecSetValue(ClosestEleShapeFun_, mpirank_* int(nPoints_)+ i,  CDist, INSERT_VALUES);
                    CHKERRQ(ierr);
                    
                    /// points that are not inside a tetrahedron are mapped to the closest centroid
                } else if ((min >= -0.1) && (max <= 1.1) && (nearC_[i][0] == INFINITY)) {
                    nearC_[i][0] = CDist;
                    nearC_[i][1] = l(0);
                    nearC_[i][2] = l(1);
                    nearC_[i][3] = l(2);
                    nearC_[i][4] = l(3);
                    
                    /// map current acc node i to solid element e
                    eleList_[i] = (e);
                    
                    ierr =
                    VecSetValue(ClosestEleShapeFun_, mpirank_* int(nPoints_)+ i,  CDist, INSERT_VALUES);
                    CHKERRQ(ierr);
                }
            }
            
            /// build gauss point i of solidElement_ e
            /// NumQP_ = 1: centroid
            /// NumQP_ = 5: centroid, QP1, QP2, QP3, QP4
            for (int QPi = 0; QPi < NumQP_; QPi++) { // NumQP_ = 1 or 5
                TFloat QP[3];
                QP[0] = ((p1.Get(0) * ShapeFunVec[QPi*NumQP_ + 0] +
                          p2.Get(0) * ShapeFunVec[QPi*NumQP_ + 1] +
                          p3.Get(0) * ShapeFunVec[QPi*NumQP_ + 2] +
                          p4.Get(0) * ShapeFunVec[QPi*NumQP_ + 3])) * 1000;
                QP[1] = ((p1.Get(1) * ShapeFunVec[QPi*NumQP_ + 0] +
                          p2.Get(1) * ShapeFunVec[QPi*NumQP_ + 1] +
                          p3.Get(1) * ShapeFunVec[QPi*NumQP_ + 2] +
                          p4.Get(1) * ShapeFunVec[QPi*NumQP_ + 3])) * 1000;
                QP[2] = ((p1.Get(2) * ShapeFunVec[QPi*NumQP_ + 0] +
                          p2.Get(2) * ShapeFunVec[QPi*NumQP_ + 1] +
                          p3.Get(2) * ShapeFunVec[QPi*NumQP_ + 2] +
                          p4.Get(2) * ShapeFunVec[QPi*NumQP_ + 3])) * 1000;
                
                /// Find closest Cell in accMesh to later interpolate Force/Cai from EP to CM
                vtkSmartPointer<vtkIdList> Points = vtkSmartPointer<vtkIdList>::New();
                vtkIdType CellId = PointLocator->FindClosestPoint(QP);
                Points = acMesh_->GetCell(CellId)->GetPointIds();
                Vector3<TFloat> acCellPoints[4];
                
                /// assign acMesh_ nodes to QP
                for (int i = 0; i < 4; i++) {
                    acCellPoints[i] = acMesh_->GetPoint(Points->GetId(i));
                    nearP_[e->GetIndex()*NumQP_ + QPi].push_back(PetscInt(Points->GetId(i)));
                }
                
                /// Determine Shape fun to interpolate force from Acc points to QP later on
                Matrix4<TFloat> mAcc = {acCellPoints[0].X(), acCellPoints[1].X(), acCellPoints[2].X(), acCellPoints[3].X(),
                    acCellPoints[0].Y(), acCellPoints[1].Y(), acCellPoints[2].Y(), acCellPoints[3].Y(),
                    acCellPoints[0].Z(), acCellPoints[1].Z(), acCellPoints[2].Z(), acCellPoints[3].Z(),
                    1,           1,           1,          1};
                mAcc.Invert();
                
                /// Gauss point of solid element e expressed with shape functions of acc element
                ForceSF_[e->GetIndex()*NumQP_+ QPi] = mAcc * Vector4<TFloat>(QP[0],  QP[1], QP[2], 1);
            } // end loop over QPs
        } else {
            skippedElement++;
        }
    } // end loop over solidElements_
    
    DCCtrl::debug << "Distribute mapping to processes...\n";
    ierr = VecAssemblyBegin(ClosestEleShapeFun_); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ClosestEleShapeFun_); CHKERRQ(ierr);
    
    VecScatter  ScatterEleDist;
    Vec localEleDist;
    
    ierr = VecScatterCreateToAll(ClosestEleShapeFun_, &ScatterEleDist, &localEleDist); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterEleDist, ClosestEleShapeFun_, localEleDist, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(
                                                                                                                       ierr);
    ierr = VecScatterEnd(ScatterEleDist, ClosestEleShapeFun_, localEleDist, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    
    PetscScalar *pCESF;
    ierr = VecGetArray(localEleDist, &pCESF); CHKERRQ(ierr);
    
    double closestPc = 1e6;
    TInt closestPci = -1;
    for (auto it = eleList_.begin(); it != eleList_.end();) {
        closestPc = 1e6;
        closestPci = -1;
        
        /// determine vertice distribution based on solidElements_
        for (TInt Pci = 0; Pci < mpisize_; Pci++) {
            if (pCESF[Pci*nPoints_+it->first] < closestPc) {
                closestPci = Pci;
                closestPc = pCESF[Pci*nPoints_+it->first];
            }
        }
        
        /// remove elements/mapping if it is not need by current process
        if (closestPci != mpirank_) {
            nearC_.erase(it->first);
            it = eleList_.erase(it);
        } else {
            it++;
        }
    }
    
    DCCtrl::debug << "Done\n";
    
    ierr = VecRestoreArray(localEleDist, &pCESF); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ScatterEleDist); CHKERRQ(ierr);
    ierr = VecDestroy(&localEleDist); CHKERRQ(ierr);
    ierr = VecDestroy(&ClosestEleShapeFun_); CHKERRQ(ierr);
} // CBacCELLerate::InitMapping

void CBacCELLerate::InitPetscVec() {
    PetscErrorCode ierr;
    
    if (DCCtrl::IsParallel())
        VecCreateMPI(Petsc::Comm(), 9 * numLocalElements_, PETSC_DETERMINE, &deformation_);
    else
        VecCreateSeq(PETSC_COMM_SELF, 9 * nCells_, &deformation_);
    
    VecZeroEntries(deformation_);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &stretchVecF_); CHKERRQ(ierr);
    ierr = VecSetSizes(stretchVecF_, PETSC_DECIDE,  int(nPoints_)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(stretchVecF_); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &accNodes_); CHKERRQ(ierr);
    ierr = VecSetSizes(accNodes_, PETSC_DECIDE,  int(nPoints_*3)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(accNodes_); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &ClosestEleShapeFun_); CHKERRQ(ierr);
    ierr = VecSetSizes(ClosestEleShapeFun_, PETSC_DECIDE, int(nPoints_ * (mpisize_))); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ClosestEleShapeFun_); CHKERRQ(ierr);
    
    ierr = VecDuplicate(stretchVecF_, &velocityVec_); CHKERRQ(ierr);
    ierr = VecDuplicate(stretchVecF_, &timestepForce_); CHKERRQ(ierr);
    ierr = VecDuplicate(stretchVecF_, &stepbackForce_); CHKERRQ(ierr);
    ierr = VecSet(stretchVecF_, 1); CHKERRQ(ierr);
    ierr = VecSet(velocityVec_, 0); CHKERRQ(ierr);
    ierr = VecSet(timestepForce_, 0); CHKERRQ(ierr);
    ierr = VecSet(stepbackForce_, 0); CHKERRQ(ierr);
    ierr = VecSet(ClosestEleShapeFun_, 1e6); CHKERRQ(ierr);
} // CBacCELLerate::InitPetscVec

void CBacCELLerate::InitPvdFile() {
    std::string timeStepsDir = resFolder_ + "/" + pvdFilename_ + "_vtuData";
    
    if (!frizzle::filesystem::CreateDirectory(resFolder_)) {
        throw(std::string("void CBModelExporterVTK::InitPvdFile(): Path: " + resFolder_ +
                          " exists but is not a directory"));
    }
    
    //  if (!frizzle::filesystem::CreateDirectory(timeStepsDir)) {
    //    throw std::runtime_error(
    //            "void CBModelExporterVTK::InitPvdFile(): Path: " + timeStepsDir + " exists but is not a directory");
    //  }
} // CBacCELLerate::InitPvdFile

void CBacCELLerate::InitLocalMaps() {
    /// populate these maps only on local processes to save memory
    for (int globalID = localElementsFrom_; globalID <= localElementsTo_; globalID++) {
        int localID = globalID - localElementsFrom_;
        
        dNdX_[localID] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        F_[localID] = Matrix3<TFloat>::Identity();
        Q_[localID] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    }
}

void CBacCELLerate::ApplySpatialSortPCA() {
    /// Adapted from the PCA sorting done in CM
    /// Most likely not compatible when used in combination with CBloadUnloadedState
    int             mDim = 3, nDim = int(nPoints_);
    int             lda  = mDim, ldu = mDim, ldvt = nDim, info, lwork;
    double          wkopt;
    double *work;
    double *sVec = new double[nDim];
    double *uVec = new double[lda * mDim];
    double *a = new double[lda * nDim];
    Vector3<TFloat> mean(0, 0, 0);
    
    backwardMapping_.resize(nPoints_);
    forwardMapping_.resize(nPoints_);
    
    for (int i = 0; i < nPoints_; i++) {
        forwardMapping_.at(i) = i;
        backwardMapping_.at(i) = i;
    }
    
    for (TInt i = 0; i < nPoints_; i++) {
        mean += acMesh_->GetPoint(i);
    }
    mean /= nPoints_;
    
    for (TInt i = 0; i < nPoints_; i++) {
        Vector3<TFloat> node = acMesh_->GetPoint(i);
        node -= mean;
        TFloat *n = node.GetArray();
        for (int k = 0; k < 3; k++) {
            a[3 * i + k] = n[k];
        }
    }
    
    lwork = -1;
    dgesvd_("S", "N", &mDim, &nDim, a, &lda, sVec, uVec, &ldu, 0, &ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work  = new double[lwork];
    dgesvd_("S", "N", &mDim, &nDim, a, &lda, sVec, uVec, &ldu, 0, &ldvt, work, &lwork, &info);
    if (info > 0) {
        throw std::runtime_error("CBModel::ApplySpatialSortPCA(): The algorithm computing SVD failed to converge");
    }
    
    double *score = new double[nPoints_];
    
    for (TInt i = 0; i < nPoints_; i++) {
        double          s    = 0;
        Vector3<TFloat> node = acMesh_->GetPoint(i);
        node -= mean;
        
        TFloat *n = node.GetArray();
        
        for (int k = 0; k < 3; k++) {
            s += n[k] * uVec[k];
        }
        score[i] = s;
    }
    
    std::vector<std::pair<TInt, double>> nodesIndexes;
    nodesIndexes.reserve(nPoints_);
    
    for (TInt i = 0; i < nPoints_; i++) {
        std::pair<TInt, double> p;
        p.first  = i;
        p.second = score[i];
        nodesIndexes.push_back(p);
    }
    
    std::sort(nodesIndexes.begin(), nodesIndexes.end(), [](std::pair<TInt, double> f, std::pair<TInt, double> b) {
        return f.second < b.second;
    });
    
    TInt *mapping = new TInt[nPoints_];
    
    std::vector<Vector3<TFloat>> sortedNodes;
    std::vector<TInt>             sortedNodesComponentsBoundaryConditions;
    
    sortedNodes.reserve(nPoints_);
    sortedNodesComponentsBoundaryConditions.reserve(nPoints_);
    
    std::vector<TInt> tmpMapping = backwardMapping_;
    
    vtkSmartPointer<vtkPoints> sortedpoints = vtkSmartPointer<vtkPoints>::New();
    
    for (TInt i = 0; i < nPoints_; i++) {
        sortedpoints->InsertNextPoint(acMesh_->GetPoint(nodesIndexes.at(i).first));
        forwardMapping_.at(tmpMapping.at(nodesIndexes.at(i).first)) = i; // using tmpMapping to map from current indices to original to update forwardMapping
        backwardMapping_.at(i) = tmpMapping.at(nodesIndexes.at(i).first); // and vice versa
        
        mapping[nodesIndexes.at(i).first] = i;
    }
    
    vtkSmartPointer<vtkTetra> tetra1 = vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType i = 0; i < nCells_; ++i) {
        vtkSmartPointer<vtkCell> c = acMesh_->GetCell(i);
        for (int k = 0; k < c->GetNumberOfPoints(); ++k) {
            TInt n = int((c)->GetPointId(k));
            tetra1->GetPointIds()->SetId(k, mapping[n]);
        }
        cellArray->InsertNextCell(tetra1);
    }
    acMesh_->SetPoints(sortedpoints);
    acMesh_->SetCells(VTK_TETRA, cellArray);
    delete[] work;
    delete[] score;
    delete[] sVec;
    delete[] uVec;
    
    delete[] a;
} // CBacCELLerate::ApplySpatialSortPCA

void CBacCELLerate::GetFibersFromMechanicsMesh() {
    PetscErrorCode ierr;
    
    Vec FiberVecX, FiberVecY, FiberVecZ;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &FiberVecX); CHKERRQ(ierr);
    ierr = VecSetSizes(FiberVecX, PETSC_DECIDE, int(nPoints_)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(FiberVecX); CHKERRQ(ierr);
    ierr = VecSet(FiberVecX, 0); CHKERRQ(ierr);
    
    ierr = VecDuplicate(FiberVecX, &FiberVecY); CHKERRQ(ierr);
    ierr = VecDuplicate(FiberVecX, &FiberVecZ); CHKERRQ(ierr);
    
    // Get fiber information for each AccNode from CMElement Centroid via eleList.
    for (auto it = eleList_.begin(); it != eleList_.end(); it++) {
        CBElementSolid *e = it->second;
        PetscInt ix = {it->first};
        PetscScalar x = {e->GetBasis()->GetCol(0)(0)};
        PetscScalar y = {e->GetBasis()->GetCol(0)(1)};
        PetscScalar z = {e->GetBasis()->GetCol(0)(2)};
        ierr = VecSetValue(FiberVecX, ix, x, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(FiberVecY, ix, y, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(FiberVecZ, ix, z, INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecAssemblyBegin(FiberVecX);
    ierr = VecAssemblyEnd(FiberVecX);
    
    ierr = VecAssemblyBegin(FiberVecY);
    ierr = VecAssemblyEnd(FiberVecY);
    
    ierr = VecAssemblyBegin(FiberVecZ);
    ierr = VecAssemblyEnd(FiberVecZ);
    
    // Create sequential Vector to write to FiberXYZArray
    VecScatter  ScatterFiberVecX, ScatterFiberVecY, ScatterFiberVecZ;
    Vec localFiberVecX, localFiberVecY, localFiberVecZ;
    
    ierr = VecScatterCreateToAll(FiberVecX, &ScatterFiberVecX, &localFiberVecX); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterFiberVecX, FiberVecX, localFiberVecX, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ScatterFiberVecX, FiberVecX, localFiberVecX, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    ierr = VecScatterCreateToAll(FiberVecY, &ScatterFiberVecY, &localFiberVecY); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterFiberVecY, FiberVecY, localFiberVecY, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ScatterFiberVecY, FiberVecY, localFiberVecY, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    ierr = VecScatterCreateToAll(FiberVecZ, &ScatterFiberVecZ, &localFiberVecZ); CHKERRQ(ierr);
    ierr = VecScatterBegin(ScatterFiberVecZ, FiberVecZ, localFiberVecZ, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ScatterFiberVecZ, FiberVecZ, localFiberVecZ, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    PetscScalar *pFiberX, *pFiberY, *pFiberZ;
    ierr = VecGetArray(localFiberVecX, &pFiberX); CHKERRQ(ierr);
    ierr = VecGetArray(localFiberVecY, &pFiberY); CHKERRQ(ierr);
    ierr = VecGetArray(localFiberVecZ, &pFiberZ); CHKERRQ(ierr);
    
    vtkSmartPointer<vtkDoubleArray> FiberXYZArray = vtkSmartPointer<vtkDoubleArray>::New();
    FiberXYZArray->SetNumberOfComponents(3);
    FiberXYZArray->Allocate(nCells_);
    FiberXYZArray->SetName("Fiber");
    
    // Set FiberOrientation for each AccCell from AccNode with ID 0 within that Cell.
    // For now, we think this is a better solution than averaging over all nodes within the Cell since this will lead to
    // unwanted effects at Cells were the FiberOrientation differs greatly.
    for (int i = 0; i < int(nCells_); i++) {
        FiberXYZArray->InsertTuple3(i, pFiberX[acMesh_->GetCell(i)->GetPointId(0)],
                                    pFiberY[acMesh_->GetCell(i)->GetPointId(0)],
                                    pFiberZ[acMesh_->GetCell(i)->GetPointId(0)]);
    }
    
    ierr = VecRestoreArray(localFiberVecX, &pFiberX); CHKERRQ(ierr);
    ierr = VecRestoreArray(localFiberVecY, &pFiberY); CHKERRQ(ierr);
    ierr = VecRestoreArray(localFiberVecZ, &pFiberZ); CHKERRQ(ierr);
    
    if (FiberXYZArray->GetNumberOfTuples() != nCells_) {
        cout << "ArrayLength: " << FiberXYZArray->GetNumberOfTuples() << "\n Number of Cells: " <<
        nCells_ << endl;
        throw std::runtime_error("ArrayLength != #Cells");
    }
    acMesh_->GetCellData()->AddArray(FiberXYZArray);
    
    ierr = VecScatterDestroy(&ScatterFiberVecX); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ScatterFiberVecY); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ScatterFiberVecZ); CHKERRQ(ierr);
    ierr = VecDestroy(&localFiberVecX); CHKERRQ(ierr);
    ierr = VecDestroy(&localFiberVecY); CHKERRQ(ierr);
    ierr = VecDestroy(&localFiberVecZ); CHKERRQ(ierr);
    ierr = VecDestroy(&FiberVecX); CHKERRQ(ierr);
    ierr = VecDestroy(&FiberVecY); CHKERRQ(ierr);
    ierr = VecDestroy(&FiberVecZ); CHKERRQ(ierr);
} // CBacCELLerate::GetFibersfromMechanicsMesh

Matrix3<TFloat> CBacCELLerate::GetInitialBasisAtCell(vtkIdType cellID) {
    Matrix3<TFloat> basis =
    { acMeshFiberValues_->GetComponent(cellID, 0), acMeshSheetValues_->GetComponent(cellID, 0),
        acMeshNormalValues_->GetComponent(cellID, 0),
        acMeshFiberValues_->GetComponent(cellID, 1), acMeshSheetValues_->GetComponent(cellID, 1),
        acMeshNormalValues_->GetComponent(cellID, 1),
        acMeshFiberValues_->GetComponent(cellID, 2), acMeshSheetValues_->GetComponent(cellID, 2),
        acMeshNormalValues_->GetComponent(cellID, 2)};
    
    return basis;
}

void CBacCELLerate::DetermineElementRanges() {
    /// This function creates a parallel layout for the acMesh_ elements
    elementRanges_.resize(DCCtrl::GetNumberOfProcesses() + 1);
    elementRanges_[0] = 0;
    PetscInt n = int(nCells_);
    
    for (unsigned int i = 1; i < DCCtrl::GetNumberOfProcesses(); i++) {
        PetscInt m = n / ((DCCtrl::GetNumberOfProcesses() - i) + 1);
        elementRanges_[i] = elementRanges_[i - 1] + m;
        n -= m;
    }
    
    elementRanges_[DCCtrl::GetNumberOfProcesses()] = int(nCells_);
    
    localElementsFrom_        = elementRanges_[DCCtrl::GetProcessID()];
    localElementsTo_          = elementRanges_[DCCtrl::GetProcessID() + 1] - 1;
    numLocalElements_         = (localElementsTo_ - localElementsFrom_) + 1;
}
