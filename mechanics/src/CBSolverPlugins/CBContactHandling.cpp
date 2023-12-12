/*
 * File: CBContactHandling.cpp
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


#include <limits>

#include <petscksp.h>

#include "filesystem.h"

#include "Vector3.h"
#include "CBContactHandling.h"
#include "CBSolver.h"
#include <algorithm>


// ----------- Public ------------

CBContactHandling::CBContactHandling() {}

void CBContactHandling::Prepare() {
    if (parameters_->Get<bool>("Plugins.ContactHandling.Initialize", true))
        status_ = CBStatus::PREPARING_SIMULATION;
    else
        status_ = CBStatus::DACCORD;
}

void CBContactHandling::Init() {
    if (parameters_->IsAvailable("Plugins.MortarContact.Alpha") ||
        parameters_->IsAvailable("Plugins.MortarContact.Beta") || parameters_->IsAvailable("Plugins.MortarContact.Steps"))
        throw std::runtime_error("CBContactHandling::Init(): MortarContact is a legacy, use ContactHandling instead");
    
    useContactHandlingToFitPeri_ = parameters_->Get<bool>("Plugins.ContactHandling.useContactHandlingToFitPeri", false);
    
    startTime_ = parameters_->Get<TFloat>("Plugins.ContactHandling.StartTime", std::numeric_limits<double>::lowest());
    if (startTime_ != std::numeric_limits<double>::lowest())
        status_ = CBStatus::DACCORD; // disable initialization
    alpha_  = parameters_->Get<TFloat>("Plugins.ContactHandling.Alpha", 5e7);       // penalizes large distance, spring
                                                                                    // stiffness Pa/m
    beta_   = parameters_->Get<TFloat>("Plugins.ContactHandling.Beta", 1);          // multiplication factor for alpha
                                                                                    // during initialization
    cntMax_ = parameters_->Get<TFloat>("Plugins.ContactHandling.Steps", 5);         // number of initialization steps
    bool flipSlaveNormals = parameters_->Get<bool>("Plugins.ContactHandling.FlipSlaveNormals", false);
    if (flipSlaveNormals)
        normalVectorSign_ = -1;
    oneWayForce_ = parameters_->Get<TInt>("Plugins.ContactHandling.OneWayForce", 0);  // 0 -> contact force in both
                                                                                      // directions, 1 -> only positive
                                                                                      // force, -1 -> only negative force
    if (std::abs(oneWayForce_) > 1) {
        throw std::runtime_error("CBContactHandling::Init(): OneWayForce '" + std::to_string(
                                                                                             oneWayForce_) + "' not valid. Has to be 0, 1 or -1.");
    }
    maxDistanceToSlave_     = parameters_->Get<TFloat>("Plugins.ContactHandling.MaxDistance", 5e-3);            // keep
                                                                                                                // contact
                                                                                                                // up to
                                                                                                                // this
                                                                                                                // distance,
                                                                                                                // m
    transitionDistance_     = parameters_->Get<TFloat>("Plugins.ContactHandling.TransitionDistance", 1e-4);     // F(d)
                                                                                                                // behaves
                                                                                                                // quadratically
                                                                                                                // until
                                                                                                                // this
                                                                                                                // distance,
                                                                                                                // linear
                                                                                                                // elsewhere,
                                                                                                                // m
    surfaceNormalDirection_ = parameters_->Get<std::string>("Plugins.ContactHandling.SurfaceNormalDirection",
                                                            "bidirectional");  // unidirectional if you only want to use
                                                                               // master surface normal
    if ((surfaceNormalDirection_ != "bidirectional") && (surfaceNormalDirection_ != "unidirectional") ) {
        throw std::runtime_error(
                                 "CBContactHandling::Init: Plugins.ContactHandling.SurfaceNormalDirection parameter not supported. Choose between unidirectional and bidirectional (default).");
    }
    TFloat maxAngleDeg = parameters_->Get<TFloat>("Plugins.ContactHandling.MaxAngle", 90.0);                // max angle
                                                                                                            // between
                                                                                                            // master and
                                                                                                            // slave
                                                                                                            // normals
                                                                                                            // that still
                                                                                                            // gets
                                                                                                            // considered,
                                                                                                            // in deg
    maxAngle_ = maxAngleDeg*M_PI/180.0;
    cnt_      = 1;
    
    initType_ = parameters_->Get<std::string>("Plugins.ContactHandling.InitType", "Exponential");
    if ((initType_ != "Exponential") && (initType_ != "Linear") ) {
        throw std::runtime_error(
                                 "CBContactHandling::Init(): InitType " + initType_ +
                                 "not found. InitType has to be one of: Exponential, Linear");
    }
    
    // if startTime_ is set in the XML, initialization is disabled. Therefore, the following code should not be executed
    // or contactPressure is 0 all the time
    if ((initType_ == "Linear") && (status_ != CBStatus::DACCORD) ) {
        maxAlpha_ = alpha_;
        alpha_    = 0;
    }
    
    filename_ = parameters_->Get<std::string>("Plugins.ContactHandling.ExportPath", "");
    export_ = parameters_->Get<bool>("Plugins.ContactHandling.Export", false);
    
    if (filename_ != "") {
        size_t pos = filename_.find_last_of("/\\");
        if (pos != std::string::npos) {
            std::string dir = filename_.substr(0, pos);
            frizzle::filesystem::CreateDirectory(dir);
        }
        
        std::ofstream file(filename_.c_str());
        if (!file.good())
            throw std::runtime_error("CBContactHandling::Init: Couldn't create " + filename_ + ".");
        file << "Time\tContactDistance\tContactPressure\n";
        file.close();
    } else {
        throw std::runtime_error("CBContactHandling::Init: Filename is empty.");
    }
    
    LoadSlaveElements();
    LoadMasterElements();
    
    TInt from = 0;
    
    // Create vector for nodes coordinates of the slave elements, these are needed on all processes -> VecScatter ...
    if (DCCtrl::IsParallel()) {
        VecCreateMPI(Petsc::Comm(), 9 * slaveElements_.size(), PETSC_DETERMINE, &slaveElementsNodes_);
        VecAssemblyBegin(slaveElementsNodes_);
        VecAssemblyEnd(slaveElementsNodes_);
        
        VecGetSize(slaveElementsNodes_, &numGlobalSlaveElements_);
        numGlobalSlaveElements_ /= 9;
        VecGetOwnershipRange(slaveElementsNodes_, &from, PETSC_NULL);
        from /= 9;
        
        VecCreateSeq(PETSC_COMM_SELF, 9 * numGlobalSlaveElements_, &slaveElementsNodesSeq_);
        VecScatterCreateToAll(slaveElementsNodes_, &scatter_, PETSC_NULL);
    } else {
        VecCreateSeq(PETSC_COMM_SELF, 9 * slaveElements_.size(), &slaveElementsNodesSeq_);
        numGlobalSlaveElements_ = slaveElements_.size();
        from                    = 0;
    }
    
    // Reindexing the slave element to a global index
    for (auto e : slaveElements_) {
        e->SetLocalIndex(from);
        from++;
    }
    
    if (DCCtrl::IsParallel()) {
        TInt *indices = new TInt[3*slaveElements_.size()];
        
        int s = 0;
        
        for (auto it : slaveElements_) {
            for (int i = 0; i < 3; i++)
                indices[3*s+i] = it->GetNodeIndex(i);
            s++;
            slaveElementsSurfaceIndices_.push_back(it->GetSurfaceIndex());
        }
        
        std::vector<TInt> tmp(9*slaveElements_.size());
        for (int i = 0; i < 3*slaveElements_.size(); i++) {
            tmp[3*i]   = 3*indices[i];
            tmp[3*i+1] = 3*indices[i]+1;
            tmp[3*i+2] = 3*indices[i]+2;
        }
        adapter_->ApplyLocalToGlobalMapping(&tmp[0], 9*slaveElements_.size());
        
        for (int i = 0; i < 3*slaveElements_.size(); i++)
            slaveElementsNodesIndicesGlobal_.push_back(tmp[3*i] / 3);
        
        
        Petsc::GatherToZero(slaveElementsNodesIndicesGlobal_);
        Petsc::CopyFromZeroToAll(slaveElementsNodesIndicesGlobal_);
        
        
        Petsc::GatherToZero(slaveElementsSurfaceIndices_);
        Petsc::CopyFromZeroToAll(slaveElementsSurfaceIndices_);
        
        if (slaveElementsSurfaceIndices_.size() != numGlobalSlaveElements_) {
            throw std::runtime_error(
                                     "void CBContactHandling::Init(): Something very bad is happening here, check the source code");
        }
        
        delete[] indices;
    } else {
        for (auto e : slaveElements_) {
            for (int i = 0; i < 3; i++)
                slaveElementsNodesIndicesGlobal_.push_back(e->GetNodeIndex(i));
            slaveElementsSurfaceIndices_.push_back(e->GetSurfaceIndex());
        }
    }
    double t1 = MPI_Wtime();
    DetermineSlaveElementsAtGaussPoints();
    DCCtrl::debug << "\nCBContactHandling::DetermineSlaveElementsAtGaussPoints [" << MPI_Wtime() - t1 << " s]\n";
    
    // ----- find neighbors of slave elements
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    
    for (int j = 0; j < numGlobalSlaveElements_; j++) {
        // DCCtrl::print << "Neighbors for " << j << ":\t";
        TFloat i10 = slaveNodes[9*j];
        TFloat i11 = slaveNodes[9*j+1];
        TFloat i12 = slaveNodes[9*j+2];
        TFloat i13 = slaveNodes[9*j+3];
        TFloat i14 = slaveNodes[9*j+4];
        TFloat i15 = slaveNodes[9*j+5];
        TFloat i16 = slaveNodes[9*j+6];
        TFloat i17 = slaveNodes[9*j+7];
        TFloat i18 = slaveNodes[9*j+8];
        
        std::vector<int> *neighbors = new std::vector<int>();
        
        for (int k = 0; k < numGlobalSlaveElements_; k++) {
            if (j == k)
                continue;
            
            TFloat i20 = slaveNodes[9*k];
            TFloat i21 = slaveNodes[9*k+1];
            TFloat i22 = slaveNodes[9*k+2];
            TFloat i23 = slaveNodes[9*k+3];
            TFloat i24 = slaveNodes[9*k+4];
            TFloat i25 = slaveNodes[9*k+5];
            TFloat i26 = slaveNodes[9*k+6];
            TFloat i27 = slaveNodes[9*k+7];
            TFloat i28 = slaveNodes[9*k+8];
            
            TFloat a = 1e-12;
            if (((std::abs(i10-i20) < a) && (std::abs(i11-i21) < a) && (std::abs(i12-i22) < a)) ||
                ((std::abs(i10-i23) < a) && (std::abs(i11-i24) < a) && (std::abs(i12-i25) < a)) ||
                ((std::abs(i10-i26) < a) && (std::abs(i11-i27) < a) && (std::abs(i12-i28) < a)) ||
                ((std::abs(i13-i20) < a) && (std::abs(i14-i21) < a) && (std::abs(i15-i22) < a)) ||
                ((std::abs(i13-i23) < a) && (std::abs(i14-i24) < a) && (std::abs(i15-i25) < a)) ||
                ((std::abs(i13-i26) < a) && (std::abs(i14-i27) < a) && (std::abs(i15-i28) < a)) ||
                ((std::abs(i16-i20) < a) && (std::abs(i17-i21) < a) && (std::abs(i18-i22) < a)) ||
                ((std::abs(i16-i23) < a) && (std::abs(i17-i24) < a) && (std::abs(i18-i25) < a)) ||
                ((std::abs(i16-i26) < a) && (std::abs(i17-i27) < a) && (std::abs(i18-i28) < a))) {
                neighbors->push_back(k);
                
                // DCCtrl::print << k << ", ";
            }
        }
        slaveNeighbors_.insert(std::pair<int, std::vector<int> *>(j, neighbors));
        
        // DCCtrl::print << std::endl;
    }
    
    // -----
    
    //    //----- find neighbors of slave elements
    //    for(auto e1 : slaveElements_)
    //    {
    //        //DCCtrl::print << "Neighbors for " << e1->GetLocalIndex() << ":\t";
    //        TInt i10 = e1->GetNodeIndex(0);
    //        TInt i11 = e1->GetNodeIndex(1);
    //        TInt i12 = e1->GetNodeIndex(2);
    //
    //        std::vector<int>* neighbors = new std::vector<int>();
    //
    //        for(auto e2 : slaveElements_)
    //        {
    //            if(e1 == e2)
    //                continue;
    //
    //            TInt i20 = e2->GetNodeIndex(0);
    //            TInt i21 = e2->GetNodeIndex(1);
    //            TInt i22 = e2->GetNodeIndex(2);
    //
    //            if(i10 == i20 || i10 == i21 || i10 == i22 ||
    //               i11 == i20 || i11 == i21 || i11 == i22 ||
    //               i12 == i20 || i12 == i21 || i12 == i22)
    //            {
    //                neighbors->push_back(e2->GetLocalIndex());
    //                //DCCtrl::print << e2->GetLocalIndex() << ", ";
    //            }
    //        }
    //        slaveNeighbors_.insert(std::pair<int, std::vector<int>*>(e1->GetLocalIndex(), neighbors));
    //        //DCCtrl::print << std::endl;
    //    }
    //    //-----
}  // CBContactHandling::Init

std::set<TInt> CBContactHandling::GetMasterNodesLocalIndices() {
    if (masterNodesLocalIndices_.size() == 0) {
        throw std::runtime_error(
                                 "void CBContactHandling::GetMasterNodesLocalIndices(): Something very bad is happening here, check the source code");
    } else {
        return masterNodesLocalIndices_;
    }
}

std::set<TInt> CBContactHandling::GetMasterWithSlaveNodesLocalIndices() {
    DetermineSlaveElementsAtVertices();
    std::set<TInt> masterNodesWithSlaveNodesLocalIndices;
    for (auto e : masterElements_)
        for (int i = 0; i < e->GetNumberOfNodesIndices(); i++)
            if (e->GetSlaveAtVertex(i) != -1)
                masterNodesWithSlaveNodesLocalIndices.insert(e->GetNodeIndex(i));
    
    return masterNodesWithSlaveNodesLocalIndices;
}

void CBContactHandling::StepBack() {
    stepBack_ = true;
}

void CBContactHandling::Apply(TFloat time) {
    if (!hasStarted_) {
        if (time > startTime_)
            hasStarted_ = true;
        else
            return;
    }
    
    prevDt_ = dt_;
    dt_     = Base::adapter_->GetSolver()->GetTiming().GetTimeStep();
    
    if (!stepBack_ || isFirstStep_) {
        double t1 = MPI_Wtime();
        DetermineSlaveElementsAtGaussPoints();
        DCCtrl::debug << "\nCBContactHandling::DetermineSlaveElementsAtGaussPoints [" << MPI_Wtime() - t1 << " s]\n";
        
        //        t1 = MPI_Wtime();
        //        DetermineSlaveElementsAtVertices();
        //        DCCtrl::debug << "CBContactHandling::DetermineSlaveElementsAtVertices [" << MPI_Wtime() - t1 << " s]\n";
    }
    
    if (status_ == CBStatus::PREPARING_SIMULATION) {
        if (time >= lastTime_) {
            if (cnt_ < cntMax_) {
                alpha_ *= beta_;
                if (initType_ == "Linear") {
                    alpha_ = maxAlpha_*(cnt_+1)/cntMax_;
                }
                
                // this line does not work well together with the LoadUnloadedState Plugin [lb451]
                // GetAdapter()->GetSolver()->RelaxElementsAndBases();
            } else {
                alpha_  = maxAlpha_;
                status_ = CBStatus::DACCORD;
            }
            cnt_++;
        }
        
        // step back case
        else if (cnt_ < cntMax_) {
            alpha_ /= beta_;
            if (initType_ == "Linear") {
                alpha_ = maxAlpha_*(cnt_-1)/cntMax_;
            }
            cnt_--;
        }
    }
    
    lastTime_    = time;
    isFirstStep_ = false;
    stepBack_    = false;
}  // CBContactHandling::Apply

double CBContactHandling::GetPreparationProgress() {
    return double(cnt_)/cntMax_;
}

void CBContactHandling::ApplyToNodalForces() {
    if (!IsActive() || !hasStarted_)
        return;
    
    DetermineSlaveNodes();
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    
    
    for (int i = 0; i < masterElements_.size(); i++) {
        auto e = masterElements_.at(i);
        
        TInt pos[18];
        bool bc[18]            = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        TFloat nodalForces[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        TFloat distances[3]    = {0, 0, 0};
        TFloat nodesCoords[9];
        e->GetNodesCoords(nodesCoords);
        TFloat scaling = e->GetSurfaceTractionScaling();
        
        Triangle<TFloat> masterTriangle(nodesCoords);
        
        masterContactForces_.at(i) = Vector3<TFloat>(0, 0, 0);
        Vector3<TFloat> slaveNormal(0, 0, 0);
        
        masterContactDistances_.at(i) = Vector3<TFloat>(0, 0, 0);
        
        for (int j = 0; j < 3; j++) {
            TInt slave = e->GetSlaveAtGaussPoint(j);
            if (slave == -1)
                continue;
            else
                masterCorrespondingSlaveFound_.at(i) = 1;
            
            Triangle<TFloat> slaveTriangle(&slaveNodes[9*slave]);
            slaveNormal += slaveTriangle.GetNormalVector();
            
            CalcContributionToContactForceAtGaussPoint(masterTriangle, slaveTriangle, j, nodalForces, distances, scaling);
            masterContactForces_.at(i)    += Vector3<TFloat>(nodalForces[3*j], nodalForces[3*j+1], nodalForces[3*j+2]);
            masterContactDistances_.at(i) += Vector3<TFloat>(distances[0], distances[1], distances[2]);
            
            for (int k = 0; k < 3; k++) {
                pos[3*k]   = 3.0 * e->GetNodeIndex(k);
                pos[3*k+1] = 3.0 * e->GetNodeIndex(k) + 1;
                pos[3*k+2] = 3.0 * e->GetNodeIndex(k) + 2;
            }
            adapter_->ApplyLocalToGlobalMapping(pos, 9);
            
            for (int k = 3; k < 6; k++) {
                pos[3*k]   = 3 * slaveElementsNodesIndicesGlobal_.at(3*slave+(k-3));
                pos[3*k+1] = 3 * slaveElementsNodesIndicesGlobal_.at(3*slave+(k-3)) + 1;
                pos[3*k+2] = 3 * slaveElementsNodesIndicesGlobal_.at(3*slave+(k-3)) + 2;
            }
            
            GetAdapter()->GetNodesComponentsBoundaryConditionsGlobal(18, pos, bc);
            
            for (int k = 0; k < 6; k++) {
                if (bc[3*k] != 0)
                    nodalForces[3*k] = 0;
                if (bc[3*k+1] != 0)
                    nodalForces[3*k+1] = 0;
                if (bc[3*k+2] != 0)
                    nodalForces[3*k+2] = 0;
            }
            
            Base::GetAdapter()->AddNodalForcesComponentsGlobal(18, pos, nodalForces);
        }
        
        if (slaveNormal.Norm() != 0)
            slaveNormal.Normalize();
        masterCorrespondingSlaveNormal_.at(i) = slaveNormal;
    }
    VecRestoreArray(slaveElementsNodesSeq_, &slaveNodes);
}  // CBContactHandling::ApplyToNodalForces

void CBContactHandling::GetMasterNodesDistancesToSlaveElements(Vec *distances) {
    if (!IsActive())
        return;
    
    DetermineSlaveNodes();
    DetermineSlaveElementsAtVertices();
    
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    Vec cnt;
    VecDuplicate(*distances, &cnt);
    
    VecZeroEntries(cnt);
    VecZeroEntries(*distances);
    
    for (auto e : masterElements_) {
        for (int i = 0; i < 3; i++) {
            if (e->GetSlaveAtGaussPoint(i) == -1)
                continue;
            PetscInt nodesIndices[3] = {3*e->GetNodeIndex(i), 3*e->GetNodeIndex(i) + 1, 3*e->GetNodeIndex(i) + 2};
            PetscScalar d[3]         =
            {e->GetDistanceVectorToSlave(i)(0), e->GetDistanceVectorToSlave(i)(1), e->GetDistanceVectorToSlave(i)(2)};
            VecSetValues(*distances, 3, nodesIndices, d, ADD_VALUES);
            PetscScalar c[3] = {1, 1, 1};
            VecSetValues(cnt, 3, nodesIndices, c, ADD_VALUES);
        }
    }
    
    VecAssemblyBegin(*distances);
    VecAssemblyEnd(*distances);
    
    VecAssemblyBegin(cnt);
    VecAssemblyEnd(cnt);
    
    PetscScalar *c;
    PetscInt l = 0;
    VecGetArray(cnt, &c);
    
    VecGetLocalSize(cnt, &l);
    
    for (int i = 0; i < l; i++)
        if (c[i] == 0)
            c[i] = 1;
    VecRestoreArray(cnt, &c);
    
    VecPointwiseDivide(*distances, *distances, cnt);
    VecAssemblyBegin(*distances);
    VecAssemblyEnd(*distances);
    VecDestroy(&cnt);
    VecRestoreArray(slaveElementsNodesSeq_, &slaveNodes);
}  // CBContactHandling::GetMasterNodesDistancesToSlaveElements

void CBContactHandling::ApplyToNodalForcesJacobian() {
    if (!IsActive() || !hasStarted_)
        return;
    
    DetermineSlaveNodes();
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    for (auto e : masterElements_) {
        TInt pos[18];
        bool bc[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        TFloat nodalForces[18];
        TFloat nodesCoords[9];
        TFloat distances[3] = {0, 0, 0};
        e->GetNodesCoords(nodesCoords);
        Triangle<TFloat> masterTriangle(nodesCoords);
        TFloat scaling = e->GetSurfaceTractionScaling();
        
        for (int i = 0; i < 3; i++) {
            TInt slave = e->GetSlaveAtGaussPoint(i);
            
            if (slave == -1)
                continue;
            
            Triangle<TFloat> slaveTriangle(&slaveNodes[9*slave]);
            
            CalcContributionToContactForceAtGaussPoint(masterTriangle, slaveTriangle, i, nodalForces, distances, scaling);
            
            for (int j = 0; j < 3; j++) {
                pos[3*j]   = 3.0 * e->GetNodeIndex(j);
                pos[3*j+1] = 3.0 * e->GetNodeIndex(j) + 1;
                pos[3*j+2] = 3.0 * e->GetNodeIndex(j) + 2;
            }
            
            adapter_->ApplyLocalToGlobalMapping(pos, 9);
            
            for (int j = 3; j < 6; j++) {
                pos[3*j]   = 3 * slaveElementsNodesIndicesGlobal_[3*slave+(j-3)];
                pos[3*j+1] = 3 * slaveElementsNodesIndicesGlobal_[3*slave+(j-3)] + 1;
                pos[3*j+2] = 3 * slaveElementsNodesIndicesGlobal_[3*slave+(j-3)] + 2;
            }
            
            GetAdapter()->GetNodesComponentsBoundaryConditionsGlobal(18, pos, bc);
            
            
            TFloat epsilon = Base::GetAdapter()->GetFiniteDifferencesEpsilon();
            TFloat nodalForcesJacobian[18*18];
            TFloat nodalForces2[18];
            TFloat slaveNodesTmp[9];
            
            
            for (int j = 0; j < 18; j++) {
                e->GetNodesCoords(nodesCoords);
                memcpy(slaveNodesTmp, &slaveNodes[9*slave], 9*sizeof(TFloat));
                
                if (j < 9)
                    nodesCoords[j] += epsilon;
                else
                    slaveNodesTmp[j-9] += epsilon;
                
                masterTriangle.SetNodes(nodesCoords);
                slaveTriangle.SetNodes(slaveNodesTmp);
                
                CalcContributionToContactForceAtGaussPoint(masterTriangle, slaveTriangle, i, nodalForces2, distances, scaling);
                
                for (int k = 0; k < 18; k++) {
                    if ((bc[k] == 0) && (bc[j] == 0))
                        nodalForcesJacobian[18*k + j] = (nodalForces2[k] - nodalForces[k]) / epsilon;
                    else
                        nodalForcesJacobian[18*k + j] = 0;
                }
            }
            Base::GetAdapter()->AddNodalForcesJacobianEntriesGlobal(18, pos, 18, pos, nodalForcesJacobian);
        }
    }
    VecRestoreArray(slaveElementsNodesSeq_, &slaveNodes);
}  // CBContactHandling::ApplyToNodalForcesJacobian

void CBContactHandling::Export(TFloat time) {
    if (export_) {
        if (DCCtrl::IsParallel()) {
            DCCtrl::WeightedAverage(averageDist_, double(masterElements_.size()), globalAverageDist_);
        } else {
            globalAverageDist_ = averageDist_;
        }
        
        Vec contactForce;
        Vec contactPressure;
        Vec contactSlaveFound;
        
        Petsc::CreateVector(3*GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &contactForce);
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &contactPressure);
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &contactSlaveFound);
        
        /// Try to export better distances: 1 CellData array per QP
        Vec distance;
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &distance);
        
        PetscInt from4, to4;
        VecGetOwnershipRange(distance, &from4, &to4);
        
        PetscInt from1, to1;
        PetscInt from2, to2;
        PetscInt from3, to3;
        VecGetOwnershipRange(contactForce, &from1, &to1);
        VecGetOwnershipRange(contactPressure, &from2, &to2);
        VecGetOwnershipRange(contactSlaveFound, &from3, &to3);
        
        VecZeroEntries(contactPressure);
        VecZeroEntries(contactForce);
        VecZeroEntries(contactSlaveFound);
        
        averageContactPressure_ = 0;
        
        for (int i = 0; i < masterElements_.size(); i++) {
            auto e = masterElements_.at(i);
            Vector3<TFloat> sn = masterCorrespondingSlaveNormal_.at(i);
            Vector3<TFloat> cf = masterContactForces_.at(i);
            
            Vector3<TFloat> dist = masterContactDistances_.at(i);
            VecSetValue(distance, from4 + e->GetLocalIndex(), (dist(0)+dist(1)+dist(2))/3.0, INSERT_VALUES);
            
            PetscInt indices[3] =
            {from1 + 3*e->GetLocalIndex(), from1 +  3*e->GetLocalIndex()+1, from1 + 3*e->GetLocalIndex()+2};
            PetscScalar f[3] = {cf(0), cf(1), cf(2)};
            
            VecSetValues(contactForce, 3, indices, f, INSERT_VALUES);
            
            PetscScalar p = cf.Norm() / e->GetArea();
            
            if (cf * sn < 0)
                p *= -1;
            
            averageContactPressure_ += p;
            VecSetValue(contactPressure, from2 + e->GetLocalIndex(), p, INSERT_VALUES);
            
            PetscScalar slaveFound = masterCorrespondingSlaveFound_.at(i);  // alternative: masterCorrespondingSlaveFound_.at(i)
                                                                            // * e->GetArea()
            VecSetValue(contactSlaveFound, from3 + e->GetLocalIndex(), slaveFound, INSERT_VALUES);
        }
        if (masterElements_.size() != 0)
            averageContactPressure_ /= masterElements_.size();
        
        if (DCCtrl::IsParallel()) {
            DCCtrl::WeightedAverage(averageContactPressure_, double(masterElements_.size()), globalAverageContactPressure_);
        } else {
            globalAverageContactPressure_ = averageContactPressure_;
        }
        
        GetAdapter()->GetSolver()->ExportElementsVectorData("ContactForce", contactForce);
        GetAdapter()->GetSolver()->ExportElementsScalarData("ContactPressure", contactPressure);
        GetAdapter()->GetSolver()->ExportElementsScalarData("ContactSlaveFound", contactSlaveFound);
        GetAdapter()->GetSolver()->ExportElementsScalarData("ContactDistance", distance);
        
        // free memory
        VecDestroy(&contactPressure);
        VecDestroy(&contactForce);
        VecDestroy(&contactSlaveFound);
        VecDestroy(&distance);
    }
}  // CBContactHandling::Export

// void CBContactHandling::GetMasterNodesDistancesToSlaveElements(Vec* d)
// {
//    Vec cnt;
//    VecDuplicate(*d, &cnt);
//    DetermineSlaveNodes();
//    DetermineSlaveElementsAtGaussPoints();
//    TFloat* slaveNodes;
//
//    VecZeroEntries(*d);
//
//    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
//    for(auto e : masterElements_)
//    {
//        TFloat nodesCoords[9];
//
//        e->GetNodesCoords(nodesCoords);
//        Triangle<TFloat> masterTriangle(nodesCoords);
//
//        for(int i = 0; i < 3; i++)
//        {
//            TInt slave = e->GetSlaveAtGaussPoint(i);
//            if(slave == -1)
//                continue;
//
//            Triangle<TFloat> slaveTriangle(&slaveNodes[9*slave]);
//
//            Vector3<TFloat> v = masterTriangle.GetNode(i);
//            Vector3<TFloat> nv =  masterTriangle.GetNormalVector();
//
//            TFloat dist = slaveTriangle.GetDistanceTo(v);
//
//            if(slaveTriangle.GetDistanceTo(v + dist*nv) > dist)
//                nv *= -1;
//
//
//            Vector3<TFloat> x = slaveTriangle.CalcIntersectionPoint(v,nv) - v;
//
//            TInt indices[3] = {3*e->GetNodeIndex(i),3*e->GetNodeIndex(i)+1,3*e->GetNodeIndex(i)+2};
//            bool* bc      = new bool[3];
//
//            GetAdapter()->GetNodesComponentsBoundaryConditions(3,indices, bc);
//
//            TFloat t[3] = {bc[0]?0:x(0),bc[1]?0:x(1),bc[2]?0:x(2)};
//
//            VecSetValues(*d, 3, indices, t, ADD_VALUES);
//            TFloat n[3] = {1.0,1.0,1.0};
//            VecSetValues(cnt, 3, indices, n, ADD_VALUES);
//
//            delete bc;
//        }
//
//    }
//    PetscScalar* c;
//    PetscInt l=0;
//    VecGetArray(cnt, &c);
//
//    VecGetLocalSize(cnt, &l);
//
//    for(int i=0; i < l; i++)
//        if(c[i]==0)
//            c[i]=1;
//    VecRestoreArray(cnt, &c);
//
//    VecPointwiseDivide(*d,*d,cnt);
//    VecAssemblyBegin(*d);
//    VecAssemblyEnd(*d);
//    VecDestroy(&cnt);

// }

void CBContactHandling::WriteToFile(TFloat time) {
    if (filename_ != "") {
        std::ofstream file(filename_.c_str(), std::ios::app);
        file << time  << " " << globalAverageDist_ << " " << globalAverageContactPressure_ << std::endl;
        file.close();
    }
}

// ----------- Private ------------

void CBContactHandling::DetermineSlaveNodes() {
    for (auto e : slaveElements_) {
        TInt pos[9];
        TFloat nodesCoords[9];
        
        for (int i = 0; i < 9; i++)
            pos[i] = 9* e->GetLocalIndex()+i;
        
        e->GetNodesCoords(nodesCoords);
        
        if (DCCtrl::IsParallel())
            VecSetValues(slaveElementsNodes_, 9, pos, nodesCoords, INSERT_VALUES);
        else
            VecSetValues(slaveElementsNodesSeq_, 9, pos, nodesCoords, INSERT_VALUES);
    }
    if (DCCtrl::IsParallel()) {
        VecAssemblyBegin(slaveElementsNodes_);
        VecAssemblyEnd(slaveElementsNodes_);
        VecScatterBegin(scatter_, slaveElementsNodes_, slaveElementsNodesSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter_, slaveElementsNodes_, slaveElementsNodesSeq_, INSERT_VALUES, SCATTER_FORWARD);
    }
}

bool CBContactHandling::CheckIfSlave(TFloat *slaveNodes, int slaveInd, Vector3<TFloat> *p, Vector3<TFloat> *nv,
                                     TFloat &dist) {
    Triangle<TFloat> slaveTriangle(&slaveNodes[9*slaveInd]);
    Vector3<TFloat>  ip = slaveTriangle.CalcIntersectionPoint(*p, *nv);
    
    if (slaveTriangle.IsPointWithinTriangle(ip)) {
        dist = (ip-*p).Norm();
        if (dist <= maxDistanceToSlave_) {
            Vector3<TFloat> snv = slaveTriangle.GetNormalVector();
            if (*nv * snv * normalVectorSign_ > 0)
                return true;
        }
    }
    return false;
}

int CBContactHandling::SearchForSlave(TFloat *slaveNodes, int oldSlave, Vector3<TFloat> *p, Vector3<TFloat> *nv,
                                      TFloat &dist) {
    bool isNewSlave;
    
    // search all slave elements if no old slave is available
    if (oldSlave == -1) {
        for (int j = 0; j < numGlobalSlaveElements_; j++) {
            isNewSlave = CheckIfSlave(slaveNodes, j, p, nv, dist);
            if (isNewSlave)
                return j;
        }
        return -1;
    }
    
    // check for the old slave first
    isNewSlave = CheckIfSlave(slaveNodes, oldSlave, p, nv, dist);
    if (isNewSlave)
        return oldSlave;
    
    // search by moving from the old slave outwards in "rings of elements"
    std::vector<int> prevStartSlaves;
    std::vector<int> startSlaves;
    std::vector<int> checkedSlaves;
    
    startSlaves.push_back(oldSlave);
    
    int depth        = 1;
    int numNeighbors = 1;
    
    while (numNeighbors > 0 && depth <= maxDepth_) {
        // std::cout << std::endl << oldSlave << "\t\tdepth " << depth << ": ";
        
        numNeighbors = 0;
        
        for (auto startSlave : startSlaves) {
            for (auto neighbor : *(slaveNeighbors_.at(startSlave))) {
                if ((std::find(prevStartSlaves.begin(), prevStartSlaves.end(), neighbor) == prevStartSlaves.end()) &&
                    (std::find(startSlaves.begin(), startSlaves.end(), neighbor) == startSlaves.end()) &&
                    (std::find(checkedSlaves.begin(), checkedSlaves.end(), neighbor) == checkedSlaves.end())) {
                    // std::cout << neighbor << " " << std::flush;
                    
                    checkedSlaves.push_back(neighbor);
                    isNewSlave = CheckIfSlave(slaveNodes, neighbor, p, nv, dist);
                    if (isNewSlave) {  // check if neighbor is slave of master
                                       // if(depth > 1)
                                       //    DCCtrl::print << "\nFound new slave at depth " << depth << "\n";
                        return neighbor;
                    }
                    
                    numNeighbors++;
                }
            }
        }
        
        prevStartSlaves = startSlaves;
        startSlaves     = checkedSlaves;
        checkedSlaves.clear();
        
        depth++;
    }
    
    // DCCtrl::print << "\nNo new slave found!\n";
    return -1;
}  // CBContactHandling::SearchForSlave

void CBContactHandling::DetermineSlaveElementsAtGaussPoints() {
    DetermineSlaveNodes();
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    
    averageDist_ = 0;
    TInt numContacts = 0;
    
    for (auto e : masterElements_) {
        TFloat nodesCoords[9];
        
        e->GetNodesCoords(nodesCoords);
        
        Triangle<TFloat> masterTriangle(nodesCoords);
        Vector3<TFloat>  nv =  masterTriangle.GetNormalVector();
        nv.Normalize();
        
        for (int i = 0; i < 3; i++) {
            Vector3<TFloat> gp = masterTriangle.GetGaussPoint3(i);
            int oldSlave       = e->GetSlaveAtGaussPoint(i);
            TFloat dist;
            int slave = SearchForSlave(slaveNodes, oldSlave, &gp, &nv, dist);
            
            if (slave != -1) {
                e->SetDistanceVectorToSlave(i, nv*dist);
                averageDist_ += dist;
                numContacts++;
            } else {
                e->SetDistanceVectorToSlave(i, Vector3<TFloat>(0, 0, 0));
            }
            
            e->SetSlaveAtGaussPoint(i, slave);
        }
    }
    if (numContacts != 0)
        averageDist_ /= numContacts;
    DCCtrl::WeightedAverage(averageDist_, double(masterElements_.size()), globalAverageDist_);
    
    if (status_ == CBStatus::PREPARING_SIMULATION) {
        DCCtrl::debug << "\t--- --- Contact Handling --- ---" << std::endl;
        DCCtrl::debug << "\tContact Force Paramter: " << alpha_ << std::endl;
        DCCtrl::debug << "\tAvg. gap = " << globalAverageDist_ << std::endl;
        DCCtrl::debug << "\t--- --- Contact Handling --- ---" << std::endl;
    }
    VecRestoreArray(slaveElementsNodesSeq_, &slaveNodes);
}  // CBContactHandling::DetermineSlaveElementsAtGaussPoints

void CBContactHandling::DetermineSlaveElementsAtVertices() {
    DetermineSlaveNodes();
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    
    averageDist_ = 0;
    TInt numContacts = 0;
    
    for (auto e : masterElements_) {
        TFloat nodesCoords[9];
        
        e->GetNodesCoords(nodesCoords);
        
        Triangle<TFloat> masterTriangle(nodesCoords);
        Vector3<TFloat>  nv =  masterTriangle.GetNormalVector();
        nv.Normalize();
        
        for (int i = 0; i < 3; i++) {
            Vector3<TFloat> vertex = masterTriangle.GetNode(i);
            int oldSlave           = e->GetSlaveAtVertex(i);
            TFloat dist;
            int slave = SearchForSlave(slaveNodes, oldSlave, &vertex, &nv, dist);
            
            if (slave != -1) {
                e->SetDistanceVectorToSlave(i, nv*dist);
                averageDist_ += dist;
                numContacts++;
            } else {
                e->SetDistanceVectorToSlave(i, Vector3<TFloat>(0, 0, 0));
            }
            
            e->SetSlaveAtVertex(i, slave);
        }
    }
    if (numContacts != 0)
        averageDist_ /= numContacts;
    DCCtrl::WeightedAverage(averageDist_, double(masterElements_.size()), globalAverageDist_);
    
    if (status_ == CBStatus::PREPARING_SIMULATION) {
        DCCtrl::debug << "\t--- --- Contact Handling --- ---" << std::endl;
        DCCtrl::debug << "\tContact Force Paramter: " << alpha_ << std::endl;
        DCCtrl::debug << "\tAvg. gap = " << globalAverageDist_ << std::endl;
        DCCtrl::debug << "\t--- --- Contact Handling --- ---" << std::endl;
    }
    VecRestoreArray(slaveElementsNodesSeq_, &slaveNodes);
}  // CBContactHandling::DetermineSlaveElementsAtVertices

void CBContactHandling::LoadMasterElements() {
    for (auto &it : Base::GetAdapter()->GetElementVector()) {
        CBElementContactMaster *element = dynamic_cast<CBElementContactMaster *>(it);
        if (element != 0) {
            int *indices = new int[3*element->GetNumberOfNodesIndices()];
            bool *bc     = new bool[3*element->GetNumberOfNodesIndices()];
            
            for (int i = 0; i < element->GetNumberOfNodesIndices(); i++) {
                auto index     = it->GetNodeIndex(i);
                indices[3*i]   = 3*index;
                indices[3*i+1] = 3*index+1;
                indices[3*i+2] = 3*index+2;
                masterNodesLocalIndices_.insert(it->GetNodeIndex(i));
            }
            
            GetAdapter()->GetNodesComponentsBoundaryConditions(3*element->GetNumberOfNodesIndices(), indices, bc);
            
            bool discard = false;
            if (!useContactHandlingToFitPeri_) {
                for (int i = 0; i < 3*element->GetNumberOfNodesIndices(); i++) {
                    if (bc[i] == true)
                        discard = true;
                }
            }
            
            if (!discard) {
                masterElements_.push_back(element);
                masterContactForces_.push_back(Vector3<TFloat>(0, 0, 0));
                masterCorrespondingSlaveNormal_.push_back(Vector3<TFloat>(0, 0, 0));
                masterCorrespondingSlaveFound_.push_back(0);
                masterContactDistances_.push_back(Vector3<TFloat>(0, 0, 0));
            }
            
            delete[] indices;
            delete[] bc;
        }
    }
}  // CBContactHandling::LoadMasterElements

void CBContactHandling::LoadSlaveElements() {
    for (auto &it : Base::GetAdapter()->GetElementVector()) {
        CBElementContactSlave *element = dynamic_cast<CBElementContactSlave *>(it);
        if (element)
            slaveElements_.push_back(element);
    }
}

void CBContactHandling::UpdateDistancesMasterSlave() {
    DetermineSlaveNodes();
    TFloat *slaveNodes;
    VecGetArray(slaveElementsNodesSeq_, &slaveNodes);
    
    for (auto e : masterElements_) {
        TFloat nodesCoords[9];
        
        e->GetNodesCoords(nodesCoords);
        Triangle<TFloat> masterTriangle(nodesCoords);
        TFloat dist = 0;
        
        for (int i = 0; i < 3; i++) {
            TInt slave = e->GetSlaveAtGaussPoint(i);
            if (slave == -1)
                continue;
            
            Triangle<TFloat> slaveTriangle(&slaveNodes[9*slave]);
            
            Vector3<TFloat> gp = masterTriangle.GetGaussPoint3(i);
            dist += slaveTriangle.GetDistanceTo(gp);
        }
        
        dist /= 3.0;
        e->SetDistanceToSlave(dist);
    }
}  // CBContactHandling::UpdateDistancesMasterSlave

void CBContactHandling::CalcContributionToContactForceAtGaussPoint(const Triangle<TFloat> &masterTriangle,
                                                                   const Triangle<TFloat> &slaveTriangle,
                                                                   TInt gaussPointIndex, TFloat *nodalForces,
                                                                   TFloat *distances, TFloat scaling) {
    Vector3<TFloat> gp   = masterTriangle.GetGaussPoint3(gaussPointIndex);
    TFloat          dist = slaveTriangle.GetDistanceTo(gp);
    
    Vector3<TFloat> mnv = masterTriangle.GetNormalVector();
    Vector3<TFloat> snv = slaveTriangle.GetNormalVector();
    
    mnv.Normalize();
    snv.Normalize();
    
    if (slaveTriangle.GetDistanceTo(gp + dist*mnv) > dist)
        mnv *= -1;
    
    if (acos(mnv*snv) > M_PI_2)
        snv *= -1;
    
    TFloat normalAngle             = 0;
    Vector3<TFloat> forceDirection = mnv;
    if (surfaceNormalDirection_ == "bidirectional") {
        forceDirection += snv;
        forceDirection.Normalize();
        
        normalAngle = acos(std::abs(mnv*snv));  // angular difference between the two normal vectors
        if (normalAngle > maxAngle_)
            normalAngle = maxAngle_;
    }
    Vector3<TFloat> ip = slaveTriangle.CalcIntersectionPoint(gp, forceDirection);
    
    // For a contact distance d < transitionDistance_, forceMagnitude depends quadratically on d
    // For d >= transitionDistance_, forceMagnitude depends linearly on d
    // This yields a smooth abs function
    TFloat d = (ip-gp).Norm();
    
    if (((oneWayForce_ == 1) && ((ip-gp)*masterTriangle.GetNormalVector() > 0)) ||
        ((oneWayForce_ == -1) && ((ip-gp)*masterTriangle.GetNormalVector() < 0)))
        d = 0;
    
    TFloat forceMagnitude;
    if (d < transitionDistance_)
        forceMagnitude = d*d / (2*transitionDistance_);
    else
        forceMagnitude = d - transitionDistance_/2;
    
    forceMagnitude *= scaling * alpha_;
    
    TFloat masterArea = masterTriangle.GetArea();
    TFloat slaveArea  = slaveTriangle.GetArea();
    
    // force contribution at QP of element e pointing in direction N (average direction between master surface normal and
    // slave surface normal)
    // f_i = - A_e * sum[ W * N_i(l1, l2, l3) * g * n ]
    // W = 1/3
    Vector3<TFloat> contactForce = -1.0 / 3.0 * masterArea * forceDirection * forceMagnitude * (maxAngle_-normalAngle)/
    maxAngle_;
    
    // This loop adds the shape functions for each QP to the nodal forces: sum[N_i] = N_1(2/3,1/6,1/6) + N_2(1/6,2/3,1/6)
    // + N_3(1/6,1/6,2/3)
    for (int i = 0; i < 3; i++) {
        if (i == gaussPointIndex) {
            nodalForces[3*i]   = 2.0 / 3.0 * contactForce(0);
            nodalForces[3*i+1] = 2.0 / 3.0 * contactForce(1);
            nodalForces[3*i+2] = 2.0 / 3.0 * contactForce(2);
            distances[i]       = d;
        } else {
            nodalForces[3*i]   = 1.0 / 6.0 * contactForce(0);
            nodalForces[3*i+1] = 1.0 / 6.0 * contactForce(1);
            nodalForces[3*i+2] = 1.0 / 6.0 * contactForce(2);
            distances[i]       = 0.0;
        }
    }
    
    // contact force of slave elements
    for (int i = 3; i < 6; i++) {
        TFloat w = slaveTriangle.GetSubArea(ip, i-3) / slaveArea;
        nodalForces[3*i]   = -w *contactForce(0);
        nodalForces[3*i+1] = -w *contactForce(1);
        nodalForces[3*i+2] = -w *contactForce(2);
    }
}  // CBContactHandling::CalcContributionToContactForceAtGaussPoint

Vector3<TFloat> CBContactHandling::CalculateDistanceAtGaussPoint(const Triangle<TFloat> &masterTriangle,
                                                                 const Triangle<TFloat> &slaveTriangle,
                                                                 TInt                    gaussPointIndex) {
    Vector3<TFloat> gp   = masterTriangle.GetGaussPoint3(gaussPointIndex);
    TFloat          dist = slaveTriangle.GetDistanceTo(gp);
    
    Vector3<TFloat> nv =   masterTriangle.GetNormalVector();  // + slaveTriangle.GetNormalVector();
    
    if (nv.Norm() != 0)
        nv.Normalize();
    
    if (slaveTriangle.GetDistanceTo(gp + dist*nv) > dist)
        nv *= -1;
    
    Vector3<TFloat> ip = slaveTriangle.CalcIntersectionPoint(gp, nv);
    
    return ip-gp;
}

Vector3<TFloat> CBContactHandling::CalculateDistanceAtVertex(const Triangle<TFloat> &masterTriangle,
                                                             const Triangle<TFloat> &slaveTriangle, TInt index) {
    Vector3<TFloat> vertex = masterTriangle.GetNode(index);
    TFloat          dist   = slaveTriangle.GetDistanceTo(vertex);
    
    Vector3<TFloat> nv =   masterTriangle.GetNormalVector();  // + slaveTriangle.GetNormalVector();
    
    if (nv.Norm() != 0)
        nv.Normalize();
    
    if (slaveTriangle.GetDistanceTo(vertex + dist*nv) > dist)
        nv *= -1;
    
    Vector3<TFloat> ip = slaveTriangle.CalcIntersectionPoint(vertex, nv);
    
    return ip-vertex;
}
