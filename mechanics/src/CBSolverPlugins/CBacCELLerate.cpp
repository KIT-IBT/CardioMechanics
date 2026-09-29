/*
 * File: CBacCELLerate.cpp
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


#include "filesystem.h"
#include "CBacCELLerate.h"
#include "CBDataCtrl.h"
#include "CBSolver.h"
#include "CBElementSolidT4.h"
#include <iostream>
#include <vtkTetra.h>
#include <vtkCellTypes.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkStaticPointLocator.h>
#include "petscvec.h"
#include <petscviewer.h>
#include "vtksys/SystemTools.hxx"

#include <algorithm>
#include <limits>
#include <unordered_map>
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
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_); CHKERRQ(ierr);
    DCCtrl::debug << "\n--- --- acCELLerate --- ---" << std::endl;
    
    InitParameters();
    
    /// Create accelerate instance and call acCELLerate constructor with the project file
    const char *pf = accprojectFile_.c_str();
    act_ = std::make_unique<acCELLerate>();
    DCCtrl::debug << "\nLoad project file ... ";
    InitPvdFile();
    act_->LoadProject(pf);
    DCCtrl::debug << "Done\n";
    
    DCCtrl::debug << "InitMono ... ";
    act_->InitMono();
    
    InitMesh();
    InitPetscVec();
    InitMapping();
    InitLocalMaps();
    
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
    
    Vector3<TFloat> nodesCoords[4];
    Matrix3<TFloat> deformationTensor;
    for (int localEleID = 0; localEleID < numLocalElements_; localEleID++) {
        for (int j = 0; j < 4; j++) {
            nodesCoords[j] = localCoords_[localCells_[localEleID][j]];
            nodesCoords[j] /= 1000.0; // adjust mm -> m
        }
        CalcDeformationTensor(localEleID, nodesCoords, deformationTensor);
        Matrix3<TFloat> Q = GetBasisAtCell(localEleID);
        
        /// calculate deformed fibers "unloaded state"
        Vector3<TFloat> f = deformationTensor * Q.GetCol(0);
        Vector3<TFloat> s = deformationTensor * Q.GetCol(1);
        Vector3<TFloat> n = deformationTensor * Q.GetCol(2);
        
        /// do Gram Schmidt orthonormalization
        f.Normalize();
        s.Normalize();
        s = s - (f*s) * f;
        n.Normalize();
        n = n - (f*n) * f - (s*n) * s;
        
        /// set new reference basis
        Q_[localEleID] = {f(0), s(0), n(0), f(1), s(1), n(1), f(2), s(2), n(2)};
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
        for (size_t eIdx = 0; eIdx < solidElements_.size(); eIdx++) {
            CBElementSolid *e = solidElements_[eIdx];
            if (std::find(materialCoupling_.begin(), materialCoupling_.end(),
                          e->GetMaterialIndex()) != materialCoupling_.end()) {
                for (int QPi = 0; QPi < NumQP_; QPi++) {
                    const QPMapping &qpMapping = qpMappings_[eIdx*NumQP_ + QPi];
                    Vector4<TFloat> ForceV4 = {pif[qpMapping.points[0]],
                        pif[qpMapping.points[1]],
                        pif[qpMapping.points[2]],
                        pif[qpMapping.points[3]]};
                    TFloat Force = ForceV4* qpMapping.shapeFun;
                    
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
            
            for (size_t eIdx = 0; eIdx < solidElements_.size(); eIdx++) {
                CBElementSolid *e = solidElements_[eIdx];
                if (find(materialCoupling_.begin(), materialCoupling_.end(),
                         e->GetMaterialIndex()) != materialCoupling_.end()) {
                    for (int QPi = 0; QPi < NumQP_; QPi++) {
                        const QPMapping &qpMapping = qpMappings_[eIdx*NumQP_ + QPi];
                        Vector4<TFloat> ForceV4sF = {ptsF[qpMapping.points[0]],
                            ptsF[qpMapping.points[1]],
                            ptsF[qpMapping.points[2]],
                            ptsF[qpMapping.points[3]]};
                        TFloat timeStepForce = ForceV4sF* qpMapping.shapeFun;
                        
                        Vector4<TFloat> ForceV4bF = {psbF[qpMapping.points[0]],
                            psbF[qpMapping.points[1]],
                            psbF[qpMapping.points[2]],
                            psbF[qpMapping.points[3]]};
                        TFloat stepBackForce = ForceV4bF* qpMapping.shapeFun;
                        
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
        
        DCPetsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &potential);
        DCPetsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &calcium);
        
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
        
        for (size_t eIdx = 0; eIdx < solidElements_.size(); eIdx++) {
            CBElementSolid *e = solidElements_[eIdx];
            if (std::find(materialCoupling_.begin(), materialCoupling_.end(),
                          e->GetMaterialIndex()) != materialCoupling_.end()) {
                for (int QPi = 0; QPi < NumQP_; QPi++) {
                    const QPMapping &qpMapping = qpMappings_[eIdx*NumQP_ + QPi];
                    Vector4<TFloat> CalciumV4 = {piC[qpMapping.points[0]],
                        piC[qpMapping.points[1]],
                        piC[qpMapping.points[2]],
                        piC[qpMapping.points[3]]};
                    TFloat Calcium = CalciumV4 * qpMapping.shapeFun;
                    
                    Vector4<TFloat> PotentialV4 = {piV[qpMapping.points[0]],
                        piV[qpMapping.points[1]],
                        piV[qpMapping.points[2]],
                        piV[qpMapping.points[3]]};
                    TFloat Potential = PotentialV4 * qpMapping.shapeFun * 1000;
                    
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
    
    for (auto &pointMapping : pointMappings_) {
        pointMapping.element->GetDeformationTensorAtQuadraturePoints(f);
        
        const Vector4<TFloat> &LocalPos = pointMapping.shapeFun;
        Vector4<TFloat> StretchAtQP = Vector4<TFloat>(sqrt(f[0].GetCol(0)*f[0].GetCol(0)),
                                                      sqrt(f[1].GetCol(0)*f[1].GetCol(0)),
                                                      sqrt(f[2].GetCol(0)*f[2].GetCol(0)),
                                                      sqrt(f[3].GetCol(0)*f[3].GetCol(0)));
        
        TFloat StretchVal = LocalPos * (QPInv * StretchAtQP);
        ierr = VecSetValue(stretchVecF_, pointMapping.point, StretchVal, INSERT_VALUES);
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
    PetscErrorCode ierr;
    
    ierr = VecSet(accNodes_, 0); CHKERRQ(ierr);
    Vector4<TFloat> l;
    
    if (useReferenceNodes)
        DCCtrl::debug << "Resetting Nodes to reference configuration ... ";
    else
        DCCtrl::debug << "Updating Nodes ... ";
    
    for (auto &pointMapping : pointMappings_) {
        CBElementSolid *e = pointMapping.element;
        const Vector4<TFloat> &shapeFun = pointMapping.shapeFun;
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
        p(0) = (p1.X() * shapeFun(0) +
                p2.X() * shapeFun(1) +
                p3.X() * shapeFun(2) +
                p4.X() * shapeFun(3)) * 1000;
        
        p(1) = (p1.Y() * shapeFun(0) +
                p2.Y() * shapeFun(1) +
                p3.Y() * shapeFun(2) +
                p4.Y() * shapeFun(3)) * 1000;
        
        p(2) = (p1.Z() * shapeFun(0) +
                p2.Z() * shapeFun(1) +
                p3.Z() * shapeFun(2) +
                p4.Z() * shapeFun(3)) * 1000;
        
        
        ierr = VecSetValue(accNodes_, pointMapping.point * 3 + 0, p(0), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(accNodes_, pointMapping.point * 3 + 1, p(1), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(accNodes_, pointMapping.point * 3 + 2, p(2), INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecAssemblyBegin(accNodes_); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(accNodes_); CHKERRQ(ierr);
    
    ierr = VecScatterBegin(localNodesScatter_, accNodes_, localNodes_, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(localNodesScatter_, accNodes_, localNodes_, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    const PetscScalar *plN;
    ierr = VecGetArrayRead(localNodes_, &plN); CHKERRQ(ierr);
    
    /// the deformed EP geometry is held in single precision; the reference results of the coupled model depend on it
    for (size_t k = 0; k < localPoints_.size(); k++)
        localCoords_[k] = Vector3<TFloat>(float(plN[3*k+0]), float(plN[3*k+1]), float(plN[3*k+2]));
    
    ierr = VecRestoreArrayRead(localNodes_, &plN); CHKERRQ(ierr);
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::UpdateNodes

void CBacCELLerate::CalcShapeFunctionDeriv() {
    /// this function is only to be called during Prepare()
    Vector3<TFloat> nodesCoords[4];
    
    DCCtrl::debug << "Updating shape function derivatives ... ";
    
    for (int localEleID = 0; localEleID < numLocalElements_; localEleID++) {
        for (int i = 0; i < 4; i++) {
            nodesCoords[i] = localCoords_[localCells_[localEleID][i]];
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
    float K[256];
    static const double Frequency = 0.0;
    
    /// list of active materials in acMesh_
    MaterialListe materialProperties(MaterialFileName_.c_str());
    
    DCCtrl::debug << "Assemble matrices ... ";
    
    /// find conductivities in each material
    for (int mat = 0; mat < 256; ++mat) {
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
    Vector3<TFloat> nodesCoords[4];
    for (int localEleID = 0; localEleID < numLocalElements_; localEleID++) {
        int material = localMaterials_[localEleID];
        if (abs(K[material] - (-1.0)) < 1e-6) {
            throw std::runtime_error("\n\nMaterial " + std::to_string(
                                                                      material) + " does not exist in the material file " + MaterialFileName_ + ".");
        }
        
        /// loop over vertices
        TInt p[4];
        for (int i = 0; i < 4; i++) {
            nodesCoords[i] = localCoords_[localCells_[localEleID][i]];
            nodesCoords[i] /= 1000.0; // adjust mm -> m
            p[i] = int(localPoints_[localCells_[localEleID][i]]); // global point index for matrix assembly
        } // end loop over vertices
        
        /// calc diffusion tensor
        Matrix3<TFloat> Q = GetBasisAtCell(localEleID);
        Matrix3<TFloat> D = { K[material] * AnisotropyX[material], 0, 0,
            0, K[material] * AnisotropyY[material], 0,
            0, 0, K[material] * AnisotropyZ[material] };
        
        Matrix3<TFloat> F;
        CalcDeformationTensor(localEleID, nodesCoords, F);
        double J = F.Det();
        
        if (J <= 0) {
            DCCtrl::debug << "\nCorrupt element ID: " << localElementsFrom_ + localEleID << "\n";
            J = 1;
            F = Matrix3<TFloat>::Identity();
        }
        
        if (MEF_ == "MINIMAL") {
            // no MEF on D
            D = Q * D * Q.GetTranspose();
            D = J * F.GetInverse() * D * F.GetInverse().GetTranspose();
        } else if (MEF_ == "FULL") {
            /// calculate deformed fibers
            Vector3<TFloat> f = F * Q.GetCol(0);
            Vector3<TFloat> s = F * Q.GetCol(1);
            Vector3<TFloat> n = F * Q.GetCol(2);
            
            D =
            D(0,
              0) *
            DyadicProduct(f,
                          f)/(f.Norm()*f.Norm()) +
            D(1, 1) * DyadicProduct(s, s)/(s.Norm()*s.Norm()) + D(2, 2) * DyadicProduct(
                                                                                        n, n)/(n.Norm()*n.Norm());
            
            // rotate
            D = J * F.GetInverse() * D * F.GetInverse().GetTranspose();
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

    priorityVector_ = parameters_->GetArray<TInt>("Plugins.acCELLerate.TissuePriority", priorityVector_);
    if (priorityVector_.empty()) {
        for (TInt i = 0; i < 256; i++)
            priorityVector_.push_back(i);
        DCCtrl::debug << "\n Tissue Priority: 0, 1, 2, 3, 4, ... (default)";
    } else {
        DCCtrl::debug << "\n Tissue Priority: ";
        for (std::vector<TInt>::const_iterator i = priorityVector_.begin(); i != priorityVector_.end(); ++i)
            DCCtrl::debug << *i << ", ";
    }
    
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
    
    /// only process 0 reads the mesh, so that memory does not grow with the number of processes;
    /// InitMapping() and InitLocalMaps() distribute what the other processes need
    std::string error;
    if (mpirank_ == 0) {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(acMeshFilename.c_str());
        reader->Update();
        acMesh_ = reader->GetOutput();
        nPoints_ = acMesh_->GetNumberOfPoints();
        nCells_ = acMesh_->GetNumberOfCells();
        
        acMeshFiberValues_ = acMesh_->GetCellData()->GetArray("Fiber");
        acMeshSheetValues_ = acMesh_->GetCellData()->GetArray("Sheet");
        acMeshNormalValues_ = acMesh_->GetCellData()->GetArray("Sheetnormal");
        acMeshMaterials_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Material")));
        
        /// the coupling reads every cell as a linear tetrahedron, so any other cell would be silently misread
        vtkIdType cellId = 0;
        while (cellId < nCells_ && acMesh_->GetCellType(cellId) == VTK_TETRA)
            cellId++;
        
        if (cellId < nCells_) {
            error = "CBacCELLerate::InitMesh(): " + acMeshFilename + " contains a cell of type " +
                    vtkCellTypes::GetClassNameFromTypeId(acMesh_->GetCellType(cellId)) +
                    " (cell " + std::to_string(cellId) + "), but only linear tetrahedra (vtkTetra) are supported.";
        } else if (!acMeshFiberValues_ || !acMeshSheetValues_ || !acMeshNormalValues_) {
            error = "Automatic mapping of fibers currently not supported.";

            //    DCCtrl::debug << "Getting Fibres from CM Mesh ... ";
            //    GetFibersFromMechanicsMesh();
            //    acMeshFiberValues_ = vtkDoubleArray::SafeDownCast((acMesh_->GetCellData()->GetArray("Fiber")));
        } else if (!acMeshMaterials_) {
            error = "CBacCELLerate::InitMesh(): No Material defined.";
        } else if (GetParameters()->Get<bool>("Plugins.acCELLerate.Permute", false)) {
            DCCtrl::debug << "\nPermuting Acc-Mesh ... ";
            ApplySpatialSortPCA();
        }
    }
    
    /// every process raises the error, otherwise the others would wait for process 0 forever
    long long header[3] = {nPoints_, nCells_, !error.empty()};
    PetscErrorCode ierr = MPI_Bcast(header, 3, MPI_LONG_LONG, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (header[2])
        throw std::runtime_error(mpirank_ == 0 ? error : "CBacCELLerate::InitMesh(): process 0 failed to load the mesh.");
    nPoints_ = header[0];
    nCells_ = header[1];
    DCCtrl::debug << "\nACMesh Properties: \nPoints: " << nPoints_ << "\nCells: " << nCells_ << "\nDone\n";
    
    DetermineElementRanges();
} // CBacCELLerate::InitMesh

template<typename T>
static MPI_Datatype CommitBytesType() {
    MPI_Datatype type;
    MPI_Type_contiguous(int(sizeof(T)), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    return type;
}

void CBacCELLerate::InitMapping() {
    DCCtrl::debug << "Creating mesh connections ... \n";
    
    /// the mapping is computed on process 0, which holds acMesh_
    std::vector<MappingElement> elements;
    for (size_t eIdx = 0; eIdx < solidElements_.size(); eIdx++) {
        CBElementSolid *e = solidElements_[eIdx];
        if (std::find(materialCoupling_.begin(), materialCoupling_.end(), e->GetMaterialIndex()) == materialCoupling_.end())
            continue;
        MappingElement element;
        element.eIdx = PetscInt(eIdx);
        PetscInt indices[12];
        for (int i = 0; i < 4; i++)
            for (int c = 0; c < 3; c++)
                indices[3*i + c] = 3 * e->GetNodeIndex(i) + c;
        GetAdapter()->GetNodesCoords(12, indices, element.coords);
        elements.push_back(element);
    }
    
    MPI_Datatype elementType = CommitBytesType<MappingElement>();
    MPI_Datatype pointType = CommitBytesType<MappingPoint>();
    MPI_Datatype qpType = CommitBytesType<MappingQP>();
    
    int numElements = int(elements.size());
    std::vector<int> elementCounts(mpisize_), elementDispls(mpisize_ + 1, 0);
    PetscErrorCode ierr = MPI_Gather(&numElements, 1, MPI_INT, elementCounts.data(), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    CHKERRQ(ierr);
    for (int r = 0; r < mpisize_; r++)
        elementDispls[r + 1] = elementDispls[r] + elementCounts[r];
    
    std::vector<MappingElement> allElements(mpirank_ == 0 ? elementDispls[mpisize_] : 0);
    ierr = MPI_Gatherv(elements.data(), numElements, elementType, allElements.data(), elementCounts.data(),
                       elementDispls.data(), elementType, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    std::vector<int> pointCounts(mpisize_, 0), pointDispls(mpisize_ + 1, 0);
    std::vector<int> qpCounts(mpisize_), qpDispls(mpisize_ + 1);
    for (int r = 0; r <= mpisize_; r++) {
        if (r < mpisize_)
            qpCounts[r] = elementCounts[r] * NumQP_;
        qpDispls[r] = elementDispls[r] * NumQP_;
    }
    std::vector<MappingPoint> allPoints;
    std::vector<MappingQP> allQPs;
    if (mpirank_ == 0)
        MapElements(allElements, elementDispls, allPoints, pointCounts, allQPs);
    for (int r = 0; r < mpisize_; r++)
        pointDispls[r + 1] = pointDispls[r] + pointCounts[r];

    /// UpdateNodes() places a point that no solid element maps at the origin, which would silently distort its cells
    long long numUnmapped = mpirank_ == 0 ? nPoints_ - pointDispls[mpisize_] : 0;
    ierr = MPI_Bcast(&numUnmapped, 1, MPI_LONG_LONG, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (numUnmapped > 0)
        throw std::runtime_error("CBacCELLerate::InitMapping(): " + std::to_string(numUnmapped) + " of " +
                                 std::to_string(nPoints_) + " acCELLerate mesh points lie outside all solid elements "
                                 "of the coupled materials (Plugins.acCELLerate.Material).");

    int numPoints;
    ierr = MPI_Scatter(pointCounts.data(), 1, MPI_INT, &numPoints, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    std::vector<MappingPoint> points(numPoints);
    ierr = MPI_Scatterv(allPoints.data(), pointCounts.data(), pointDispls.data(), pointType,
                        points.data(), numPoints, pointType, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    std::vector<MappingQP> qps(numElements * NumQP_);
    ierr = MPI_Scatterv(allQPs.data(), qpCounts.data(), qpDispls.data(), qpType,
                        qps.data(), int(qps.size()), qpType, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    MPI_Type_free(&elementType);
    MPI_Type_free(&pointType);
    MPI_Type_free(&qpType);
    
    pointMappings_.reserve(points.size());
    for (auto &p : points)
        pointMappings_.push_back({p.point, solidElements_[p.eIdx], Vector4<TFloat>(p.shapeFun)});
    
    qpMappings_.resize(solidElements_.size() * NumQP_);
    for (size_t k = 0; k < elements.size(); k++) {
        for (int QPi = 0; QPi < NumQP_; QPi++) {
            const MappingQP &qp = qps[k*NumQP_ + QPi];
            QPMapping &qpMapping = qpMappings_[elements[k].eIdx*NumQP_ + QPi];
            std::copy(qp.points, qp.points + 4, qpMapping.points);
            qpMapping.shapeFun = Vector4<TFloat>(qp.shapeFun);
        }
    }
    
    DCCtrl::debug << "Done\n";
} // CBacCELLerate::InitMapping

void CBacCELLerate::MapElements(const std::vector<MappingElement> &elements, const std::vector<int> &elementDispls,
                                std::vector<MappingPoint> &points, std::vector<int> &pointCounts,
                                std::vector<MappingQP> &qps) {
    /// shape functions for gauss points of tetrahedron
    std::vector<double> ShapeFunVec(20, 0);
    double alpha  = (5 + 3 * sqrt(5)) / 20;
    double beta   = (5 - sqrt(5)) / 20;
    
    ShapeFunVec   = { 0.25, 0.25, 0.25, 0.25,
        alpha, beta, beta, beta,
        beta, alpha, beta, beta,
        beta, beta, alpha, beta,
        beta, beta, beta, alpha };
    
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
    
    /// locates the acMesh_ points that are candidates for a solid element
    vtkSmartPointer<vtkStaticPointLocator> acPointLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
    acPointLocator->SetDataSet(acMesh_);
    acPointLocator->BuildLocator();
    vtkSmartPointer<vtkIdList> nearPoints = vtkSmartPointer<vtkIdList>::New();
    
    /// candidate of a process for a point, with the distance to the element centroid used to pick a single owner
    struct Candidate {
        double dist;
        int element;  // index into elements
        Vector4<TFloat> shapeFun;
    };
    
    /// each point is owned by the process with the closest element centroid; ties resolve to the lowest rank
    struct Owner {
        double dist = std::numeric_limits<double>::infinity();
        int rank = -1;
        int element = -1;
        Vector4<TFloat> shapeFun;
    };
    std::vector<Owner> owners(nPoints_);
    qps.resize(elements.size() * NumQP_);
    
    for (int r = 0; r < mpisize_; r++) {
        /// candidates of process r; the result depends on the order of the elements within a process
        std::unordered_map<PetscInt, Candidate> candidates;
        for (int k = elementDispls[r]; k < elementDispls[r + 1]; k++) {
            /// vertices of solid element k
            Vector3<TFloat> p1(&elements[k].coords[0]);
            Vector3<TFloat> p2(&elements[k].coords[3]);
            Vector3<TFloat> p3(&elements[k].coords[6]);
            Vector3<TFloat> p4(&elements[k].coords[9]);
            
            /// Jacobian matrix of linear tetrahedron
            Matrix4<TFloat> m = {p1.X(), p2.X(), p3.X(), p4.X(),
                p1.Y(), p2.Y(), p3.Y(), p4.Y(),
                p1.Z(), p2.Z(), p3.Z(), p4.Z(),
                1,      1,      1,      1};
            m.Invert();
            
            /// points with all shape functions >= -0.1 lie in the tetrahedron scaled by 1.4 about its centroid,
            /// which is enclosed by a sphere of 1.4 times the largest centroid-vertex distance
            Vector3<TFloat> centroid = (p1 + p2 + p3 + p4) / 4;
            double radius = std::max({(p1 - centroid).Norm(), (p2 - centroid).Norm(),
                                      (p3 - centroid).Norm(), (p4 - centroid).Norm()});
            double center[3] = {centroid.X() * 1000, centroid.Y() * 1000, centroid.Z() * 1000};
            acPointLocator->FindPointsWithinRadius(1.4 * radius * 1000 * (1 + 1e-3), center, nearPoints);
            
            for (vtkIdType n = 0; n < nearPoints->GetNumberOfIds(); n++) {
                PetscInt i = PetscInt(nearPoints->GetId(n));
                double p[3] = {acMesh_->GetPoint(i)[0], acMesh_->GetPoint(i)[1], acMesh_->GetPoint(i)[2]};
                
                /// l: point p expressed with shape functions of solidElement e
                Vector4<TFloat> l = m * Vector4<TFloat>(p[0] / 1000, p[1] / 1000, p[2] / 1000, 1);
                
                /// since the centroid of the linear tetrahedron is expressed as l={0.25,0.25,0.25,0.25}, distance to the centroid becomes minimal if l.max() - l.min() approaches 0
                double max = l.Max();
                double min = l.Min();
                double CDist = max - min;
                
                /// if point p is inside solidElement e, update map
                if ((min >= 0) && (max <= 1)) {
                    candidates[i] = {CDist, k, l};
                    
                    /// points that are not inside a tetrahedron are mapped to the closest centroid
                } else if ((min >= -0.1) && (max <= 1.1)) {
                    candidates.try_emplace(i, Candidate{CDist, k, l});
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
                vtkIdType CellId = PointLocator->FindClosestPoint(QP);
                vtkIdList *Points = acMesh_->GetCell(CellId)->GetPointIds();
                Vector3<TFloat> acCellPoints[4];
                MappingQP &qp = qps[k*NumQP_ + QPi];
                
                /// assign acMesh_ nodes to QP
                for (int i = 0; i < 4; i++) {
                    acCellPoints[i] = acMesh_->GetPoint(Points->GetId(i));
                    qp.points[i] = PetscInt(Points->GetId(i));
                }
                
                /// Determine Shape fun to interpolate force from Acc points to QP later on
                Matrix4<TFloat> mAcc = {acCellPoints[0].X(), acCellPoints[1].X(), acCellPoints[2].X(), acCellPoints[3].X(),
                    acCellPoints[0].Y(), acCellPoints[1].Y(), acCellPoints[2].Y(), acCellPoints[3].Y(),
                    acCellPoints[0].Z(), acCellPoints[1].Z(), acCellPoints[2].Z(), acCellPoints[3].Z(),
                    1,           1,           1,          1};
                mAcc.Invert();
                
                /// Gauss point of solid element e expressed with shape functions of acc element
                Vector4<TFloat> shapeFun = mAcc * Vector4<TFloat>(QP[0],  QP[1], QP[2], 1);
                for (int i = 0; i < 4; i++)
                    qp.shapeFun[i] = shapeFun(i);
            } // end loop over QPs
        } // end loop over elements of process r
        
        for (auto &c : candidates) {
            Owner &owner = owners[c.first];
            if (c.second.dist < owner.dist)
                owner = {c.second.dist, r, c.second.element, c.second.shapeFun};
        }
    } // end loop over processes
    
    /// grouped by process and sorted by point within a process
    for (auto &owner : owners)
        if (owner.rank >= 0)
            pointCounts[owner.rank]++;
    std::vector<int> next(mpisize_, 0);
    for (int r = 1; r < mpisize_; r++)
        next[r] = next[r - 1] + pointCounts[r - 1];
    points.resize(next[mpisize_ - 1] + pointCounts[mpisize_ - 1]);
    for (PetscInt i = 0; i < nPoints_; i++) {
        const Owner &owner = owners[i];
        if (owner.rank < 0)
            continue;
        MappingPoint &point = points[next[owner.rank]++];
        point.point = i;
        point.eIdx = elements[owner.element].eIdx;
        for (int j = 0; j < 4; j++)
            point.shapeFun[j] = owner.shapeFun(j);
    }
} // CBacCELLerate::MapElements

void CBacCELLerate::InitPetscVec() {
    PetscErrorCode ierr;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &stretchVecF_); CHKERRQ(ierr);
    ierr = VecSetSizes(stretchVecF_, PETSC_DECIDE,  int(nPoints_)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(stretchVecF_); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &accNodes_); CHKERRQ(ierr);
    ierr = VecSetSizes(accNodes_, PETSC_DECIDE,  int(nPoints_*3)); CHKERRQ(ierr);
    ierr = VecSetFromOptions(accNodes_); CHKERRQ(ierr);
    
    ierr = VecDuplicate(stretchVecF_, &velocityVec_); CHKERRQ(ierr);
    ierr = VecDuplicate(stretchVecF_, &timestepForce_); CHKERRQ(ierr);
    ierr = VecDuplicate(stretchVecF_, &stepbackForce_); CHKERRQ(ierr);
    ierr = VecSet(stretchVecF_, 1); CHKERRQ(ierr);
    ierr = VecSet(velocityVec_, 0); CHKERRQ(ierr);
    ierr = VecSet(timestepForce_, 0); CHKERRQ(ierr);
    ierr = VecSet(stepbackForce_, 0); CHKERRQ(ierr);
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
    /// process 0 sends every process its cells, one process at a time so that it holds a single extra block
    struct CellData {
        PetscInt points[4];
        double material;
        double basis[9];  // row-major, columns fiber, sheet, sheet normal
    };
    MPI_Datatype cellType = CommitBytesType<CellData>();
    PetscErrorCode ierr;
    
    auto packCell = [&](PetscInt cellId) {
        CellData cell;
        vtkIdList *pointIds = acMesh_->GetCell(cellId)->GetPointIds();
        for (int i = 0; i < 4; i++)
            cell.points[i] = PetscInt(pointIds->GetId(i));
        cell.material = acMeshMaterials_->GetValue(cellId);
        for (int i = 0; i < 3; i++) {
            cell.basis[3*i + 0] = acMeshFiberValues_->GetComponent(cellId, i);
            cell.basis[3*i + 1] = acMeshSheetValues_->GetComponent(cellId, i);
            cell.basis[3*i + 2] = acMeshNormalValues_->GetComponent(cellId, i);
        }
        return cell;
    };
    
    /// sorted global ids of the points used by cells [from, to), and their coordinates [mm]
    auto packPoints = [&](PetscInt from, PetscInt to, std::vector<PetscInt> &points, std::vector<double> &coords) {
        points.clear();
        points.reserve(4 * (to - from));
        for (PetscInt cellId = from; cellId < to; cellId++) {
            vtkIdList *pointIds = acMesh_->GetCell(cellId)->GetPointIds();
            for (int i = 0; i < 4; i++)
                points.push_back(PetscInt(pointIds->GetId(i)));
        }
        std::sort(points.begin(), points.end());
        points.erase(std::unique(points.begin(), points.end()), points.end());
        points.shrink_to_fit();
        coords.resize(3 * points.size());
        for (size_t k = 0; k < points.size(); k++)
            std::copy_n(acMesh_->GetPoint(points[k]), 3, &coords[3*k]);
    };
    
    /// process 0 reads its own cells directly from acMesh_
    std::vector<CellData> cells;
    std::vector<double> coords;
    if (mpirank_ == 0) {
        std::vector<PetscInt> points;
        for (int r = 1; r < mpisize_; r++) {
            std::vector<CellData> block(elementRanges_[r + 1] - elementRanges_[r]);
            for (size_t k = 0; k < block.size(); k++)
                block[k] = packCell(elementRanges_[r] + PetscInt(k));
            ierr = MPI_Send(block.data(), int(block.size()), cellType, r, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            block = {};
            packPoints(elementRanges_[r], elementRanges_[r + 1], points, coords);
            ierr = MPI_Send(points.data(), int(points.size()), MPIU_INT, r, 1, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Send(coords.data(), int(coords.size()), MPI_DOUBLE, r, 2, PETSC_COMM_WORLD); CHKERRQ(ierr);
        }
        packPoints(localElementsFrom_, localElementsTo_ + 1, localPoints_, coords);
    } else {
        cells.resize(numLocalElements_);
        ierr = MPI_Recv(cells.data(), int(numLocalElements_), cellType, 0, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
        CHKERRQ(ierr);
        MPI_Status status;
        int numPoints;
        ierr = MPI_Probe(0, 1, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);
        ierr = MPI_Get_count(&status, MPIU_INT, &numPoints); CHKERRQ(ierr);
        localPoints_.resize(numPoints);
        coords.resize(3 * size_t(numPoints));
        ierr = MPI_Recv(localPoints_.data(), numPoints, MPIU_INT, 0, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE); CHKERRQ(ierr);
        ierr = MPI_Recv(coords.data(), 3 * numPoints, MPI_DOUBLE, 0, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE); CHKERRQ(ierr);
    }
    MPI_Type_free(&cellType);
    auto localCell = [&](PetscInt localEleID) {
        return mpirank_ == 0 ? packCell(localElementsFrom_ + localEleID) : cells[localEleID];
    };
    
    dNdX_.assign(numLocalElements_, {});
    Q_.resize(numLocalElements_);
    localCells_.resize(numLocalElements_);
    localMaterials_.resize(numLocalElements_);
    for (PetscInt localEleID = 0; localEleID < numLocalElements_; localEleID++) {
        CellData cell = localCell(localEleID);
        for (int i = 0; i < 4; i++)
            localCells_[localEleID][i] = PetscInt(std::lower_bound(localPoints_.begin(), localPoints_.end(),
                                                                   cell.points[i]) - localPoints_.begin());
        localMaterials_[localEleID] = int(cell.material);
        Q_[localEleID] = Matrix3<TFloat>(cell.basis);
    }
    cells = {};
    
    std::vector<PetscInt> nodeIndices(3 * localPoints_.size());
    for (size_t k = 0; k < localPoints_.size(); k++)
        for (int c = 0; c < 3; c++)
            nodeIndices[3*k + c] = 3 * localPoints_[k] + c;
    IS nodeIS;
    ierr = ISCreateGeneral(PETSC_COMM_SELF, PetscInt(nodeIndices.size()), nodeIndices.data(), PETSC_COPY_VALUES, &nodeIS);
    CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, PetscInt(nodeIndices.size()), &localNodes_); CHKERRQ(ierr);
    ierr = VecScatterCreate(accNodes_, nodeIS, localNodes_, NULL, &localNodesScatter_); CHKERRQ(ierr);
    ierr = ISDestroy(&nodeIS); CHKERRQ(ierr);
    
    localCoords_.resize(localPoints_.size());
    for (size_t k = 0; k < localPoints_.size(); k++)
        localCoords_[k] = Vector3<TFloat>(&coords[3*k]);
    
    acMeshMaterials_ = nullptr;
    acMeshFiberValues_ = nullptr;
    acMeshSheetValues_ = nullptr;
    acMeshNormalValues_ = nullptr;
    acMesh_ = nullptr;
} // CBacCELLerate::InitLocalMaps

void CBacCELLerate::ApplySpatialSortPCA() {
    /// Adapted from the PCA sorting done in CM
    /// Most likely not compatible when used in combination with CBloadUnloadedState
    int             mDim = 3, nDim = int(nPoints_);
    int             lda  = mDim, ldu = mDim, ldvt = nDim, info, lwork;
    double          wkopt;
    std::vector<double> sVec(nDim), uVec(lda * mDim), a(lda * nDim);
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
    dgesvd_("S", "N", &mDim, &nDim, a.data(), &lda, sVec.data(), uVec.data(), &ldu, 0, &ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    std::vector<double> work(lwork);
    dgesvd_("S", "N", &mDim, &nDim, a.data(), &lda, sVec.data(), uVec.data(), &ldu, 0, &ldvt, work.data(), &lwork,
            &info);
    if (info > 0) {
        throw std::runtime_error("CBacCELLerate::ApplySpatialSortPCA(): The algorithm computing SVD failed to converge");
    }
    
    std::vector<double> score(nPoints_);
    
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
    
    std::vector<TInt> mapping(nPoints_);
    
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
    
    // Get fiber information for each AccNode from CMElement Centroid via pointMappings_.
    for (auto &pointMapping : pointMappings_) {
        CBElementSolid *e = pointMapping.element;
        PetscInt ix = {pointMapping.point};
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
