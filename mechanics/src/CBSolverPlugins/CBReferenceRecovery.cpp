/* -------------------------------------------------------
 
 CBReferenceRecovery.cpp
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach  (25.05.2023)
 Last modified: Tobias Gerach  (25.05.2023)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#include "DCCtrl.h"
#include "DCCtrlPETSc.h"

#include "CBSolver.h"
#include <math.h>
#include <vector>
#include <limits>
#include "petscvec.h"

#include "CBReferenceRecovery.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"


/// basic structure - functions called by the solver

/// main initialization function
void CBReferenceRecovery::Init() {
    // init cavities_ using inherited init function
    CBSolverPluginCavities::InitCavities("Plugins.ReferenceRecovery");
    
    InitParamsFromXML();
    InitPetscVec();
    
    for (auto s : surfaces_) {
        S_.currentPressures_.insert(std::pair<TInt, TFloat>(s, 0.0));
        S_.currentVolumes_.insert(std::pair<TInt, TFloat>(s, cavities_.at(s)->CalcVolume()));
    }
    initialVolumes_     = S_.currentVolumes_;
    S_.unloadedVolumes_ = S_.currentVolumes_;
    
    InitTargetCoordsFromCurrentCoords();
    
    InitExportPressureVolumeInfo();
    InitExportCycleInfo();
    
    // init the state variable objects: copy data to last guess (S_prev_) and last successful guess (S_curr_)
    S_prev_ = S_;
    S_curr_ = S_;
}  // CBReferenceRecovery::Init

/// apply/main function of this plugin, gets called at every time step
void CBReferenceRecovery::Apply(TFloat time) {
    // check if time proceeded and the last known correct state needs to be updated
    SaveStateVariablesAsNeeded(time);
    
    S_.time_ = time;
    
    status_ = CBStatus::SUCCESS;
    
    if (time < startTime_)
        return;
    
    DCCtrl::print << "\t--- --- Unloading --- ---" << std::endl;
    
    if ((time >= inflationStartTime_) && (time < inflationEndTime_) ) {
        DCCtrl::print << "\tInflating from time " << inflationStartTime_ << " to " << inflationEndTime_;
        CalcPressures(time);
    } else if (!S_.inflationFinished_) {
        DCCtrl::print << "\tFinishing inflation";
        CalcPressures(inflationEndTime_);
        S_.inflationFinished_ = true;
    }
    
    DCCtrl::print << "\tCurrent pressures:";
    char str[10];
    for (auto s : surfaces_) {
        sprintf(str, "%.1f", S_.currentPressures_.at(s));
        DCCtrl::print << " " << str;
    }
    DCCtrl::print << std::endl;
    
    // this should be the last line to store all computed variables
    S_curr_ = S_;
} // CBReferenceRecovery::Apply

void CBReferenceRecovery::ApplyToNodalForces() {
    for (auto s : surfaces_)
        cavities_.at(s)->ApplyPressure(S_.currentPressures_.at(s));
}

void CBReferenceRecovery::ApplyToNodalForcesJacobian() {
    //    for (auto s : surfaces_)
    //        cavities_.at(s)->ApplyPressureToJacobian(S_.currentPressures_.at(s));
}

void CBReferenceRecovery::Export(TFloat time) {
    GetAdapter()->GetSolver()->ExportNodesVectorData("Residual", S_.residual_);
}

void CBReferenceRecovery::StepBack() {
    DCCtrl::print << "\tReferenceRecovery - Step back" << std::endl;
    
    S_ = S_prev_;  // reset state variable counter
                   // S_curr_ = S_prev_; // discard guessed variables
}

void CBReferenceRecovery::AnalyzeResults() {
    for (auto s : surfaces_)
        S_.currentVolumes_.at(s) = cavities_.at(s)->CalcVolume();
    
    if (S_.inflationFinished_) {
        DCCtrl::print << "\tFinishing cycle " << S_.cycleCount_ << std::endl;
        S_.cycleFinishedTime_ = S_.time_;
        CalcResidual();
        S_.residualNorm_   = CalcResidualNorm(NORM_INFINITY);
        S_.residualNormL2_ = CalcResidualNorm(NORM_2);
        DCCtrl::print << "\tResidual norm: " << S_.residualNorm_ << std::endl;
        DCCtrl::print << "\tUnloaded volumes:";
        char str[6];
        for (auto s : surfaces_) {
            sprintf(str, "%.1f", 100*S_.unloadedVolumes_.at(s)/initialVolumes_.at(s));
            DCCtrl::print << " " << str << "%";
        }
        DCCtrl::print << std::endl;
        
        ExportCycleInfo();
        
        if (S_.incrementCount_ != numIncrements_) {
            if ((S_.residualNorm_ < tolerance_) || (S_.cycleCount_ >= cycleMax_)) {
                ExportCoordsAsNodeFile(S_.unloadedCoords_);
                ExportCoordsAsNodeFile(S_.loadedCoords_, true);
                DCCtrl::print << "\tApplying Backward Displacement" << std::endl;
                if (augmentation_ && (S_.cycleCount_ > 0) ) {
                    displacementFactor_ = UpdateAugmentationParameter(displacementFactor_);
                }
                DCCtrl::print << "\tDisplacement factor: " << displacementFactor_ << std::endl;
                ApplyBackwardDisplacement();
                ExportFibersAsBasesFile();
                inflationStartTime_   = S_.time_;
                inflationEndTime_     = inflationStartTime_ + inflationDuration_;
                S_.inflationFinished_ = false;
                S_.cycleCount_        = 0;
                S_.incrementCount_++;
            } else {
                DCCtrl::print << "\tApplying Backward Displacement" << std::endl;
                if (augmentation_ && (S_.cycleCount_ > 0) ) {
                    displacementFactor_ = UpdateAugmentationParameter(displacementFactor_);
                }
                DCCtrl::print << "\tDisplacement factor: " << displacementFactor_ << std::endl;
                ApplyBackwardDisplacement();
                inflationStartTime_   = S_.time_;
                inflationEndTime_     = inflationStartTime_ + inflationDuration_;
                S_.inflationFinished_ = false;
                S_.cycleCount_++;
            }
        } else {
            if ((S_.residualNorm_ < tolerance_) || (S_.cycleCount_ >= cycleMax_)) {
                shallExit_ = true;
                ExportCoordsAsNodeFile(S_.unloadedCoords_);
                ExportCoordsAsNodeFile(S_.loadedCoords_, true);
            } else {
                DCCtrl::print << "\tApplying Backward Displacement" << std::endl;
                if (augmentation_ && (S_.cycleCount_ > 0) ) {
                    displacementFactor_ = UpdateAugmentationParameter(displacementFactor_);
                }
                DCCtrl::print << "\tDisplacement factor: " << displacementFactor_ << std::endl;
                ApplyBackwardDisplacement();
                ExportFibersAsBasesFile();
                inflationStartTime_   = S_.time_;
                inflationEndTime_     = inflationStartTime_ + inflationDuration_;
                S_.inflationFinished_ = false;
                S_.cycleCount_++;
            }
        }
        
        S_curr_ = S_;
    }
} // CBReferenceRecovery::AnalyzeResults

void CBReferenceRecovery::WriteToFile(TFloat time) {
    ExportPressureVolumeInfo(time);
}

bool CBReferenceRecovery::ExitCheck(TFloat time) {
    return shallExit_;
}

/// private helper functions

/// convert xml -> values for class member variables, call this at initialization
void CBReferenceRecovery::InitParamsFromXML() {
    exportDir_                  =
    parameters_->Get<std::string>("Plugins.ReferenceRecovery.ExportDir", "ReferenceRecovery");
    pressureVolumeInfoFilename_ = exportDir_ + "/PressureVolumeInfo.dat";
    cycleInfoFilename_          = exportDir_ + "/CycleInfo.dat";
    
    surfaces_ = parameters_->GetArray<TInt>("Plugins.ReferenceRecovery.Surfaces");
    if (surfaces_.size() == 0)
        throw std::runtime_error("CBReferenceRecovery: At least one surface must be given.");
    for (auto s : surfaces_) {
        if (cavities_.find(s) == cavities_.end()) {
            throw std::runtime_error("CBReferenceRecovery: No cavity with surface index " + std::to_string(
                                                                                                           s) + " found. Maybe you did not use CAVITY surface elements?");
        }
    }
    
    std::vector<TFloat> pressures = parameters_->GetArray<TFloat>("Plugins.ReferenceRecovery.Pressures");
    if (pressures.size() != surfaces_.size()) {
        throw std::runtime_error(
                                 "CBReferenceRecovery: Expected the two vectors in xml ReferenceRecovery.Surfaces and ReferenceRecovery.Pressures to be of the same size.");
    }
    
    for (TInt i = 0; i < surfaces_.size(); i++) {
        targetPressures_.insert(std::pair<TInt, TFloat>(surfaces_.at(i), pressures.at(i)));
    }
    
    numIncrements_ = parameters_->Get<TInt>("Plugins.ReferenceRecovery.NumberOfIncrements", 1);
    if (numIncrements_ < 1)
        throw std::runtime_error("CBReferenceRecovery: NumberOfIncrements must be greater than or equal 1.");
    cycleMax_ = parameters_->Get<TInt>("Plugins.ReferenceRecovery.CycleMax", std::numeric_limits<TInt>::max());
    
    startTime_            = parameters_->Get<TFloat>("Plugins.ReferenceRecovery.StartTime", 0.0);
    inflationStartTime_   = startTime_;
    inflationDuration_    = parameters_->Get<TFloat>("Plugins.ReferenceRecovery.InflationDuration", 0.1);
    inflationEndTime_     = inflationStartTime_ + inflationDuration_;
    tolerance_            = parameters_->Get<TFloat>("Plugins.ReferenceRecovery.Tolerance", 1e-3);
    displacementFactor_   = parameters_->Get<TFloat>("Plugins.ReferenceRecovery.DisplacementFactor", 1.0);
    
    /// Applies the Augmentation Method proposed by M. Rausch et al. (2017): An augmented iterative method for identifying
    /// a stress-free reference configuration in image-based biomechanical modeling
    augmentation_ = parameters_->Get<bool>("Plugins.ReferenceRecovery.Augmentation", false);
} // CBReferenceRecovery::InitParamsFromXML

void CBReferenceRecovery::InitPetscVec() {
    VecCreate(Petsc::Comm(), &targetCoords_);
    VecSetSizes(targetCoords_, 3*adapter_->GetSolver()->GetNumberOfLocalNodes(), PETSC_DETERMINE);
    VecSetFromOptions(targetCoords_);
    VecZeroEntries(targetCoords_);
    
    VecScatterCreateToZero(targetCoords_, &coordsScatter_, &coordsSeq_);
    
    // init petsc vectors for the state variables objects, needs targetCoords_ to determine size of vectors
    S_.Init(targetCoords_);
    S_prev_.Init(targetCoords_);
    S_curr_.Init(targetCoords_);
}

void CBReferenceRecovery::InitTargetCoordsFromCurrentCoords() {
    adapter_->GetSolver()->GetNodeCoordinates(targetCoords_);
    VecCopy(targetCoords_, S_.unloadedCoords_);
    VecCopy(targetCoords_, S_.prevUnloadedCoords_);
    VecCopy(targetCoords_, S_.prevLoadedCoords_);
}

void CBReferenceRecovery::InitExportPressureVolumeInfo() {
    if (DCCtrl::IsProcessZero()) {
        if (!frizzle::filesystem::CreateDirectory(exportDir_)) {
            throw std::runtime_error(
                                     "CBReferenceRecovery::InitExportPressureVolumeInfo: Path: " + exportDir_ +
                                     " exists but is not a directory.");
        }
        
        
        pressureVolumeInfoFile_.open(pressureVolumeInfoFilename_.c_str());
        if (!pressureVolumeInfoFile_.good()) {
            throw std::runtime_error(
                                     "CBReferenceRecovery::InitExportPressureVolumeInfo: Couldn't create " + pressureVolumeInfoFilename_ +
                                     ".");
        }
        
        pressureVolumeInfoFile_ << "Time";
        for (auto s : surfaces_)
            pressureVolumeInfoFile_ << "    Pressure" << s << "    Volume" << s;
        pressureVolumeInfoFile_ << std::endl;
        
        pressureVolumeInfoFile_.close();
    }
}

void CBReferenceRecovery::InitExportCycleInfo() {
    if (DCCtrl::IsProcessZero()) {
        if (!frizzle::filesystem::CreateDirectory(exportDir_)) {
            throw std::runtime_error(
                                     "CBReferenceRecovery::InitExportCycleInfo: Path: " + exportDir_ + " exists but is not a directory.");
        }
        
        cycleInfoFile_.open(cycleInfoFilename_.c_str());
        if (!cycleInfoFile_.good())
            throw std::runtime_error("CBReferenceRecovery::InitExportCycleInfo: Couldn't create " + cycleInfoFilename_ + ".");
        
        cycleInfoFile_ << "Cycle    Increment    Time    ResidualNorm    ResidualNormL2";
        for (auto s : surfaces_)
            cycleInfoFile_ << "    UnloadedVolume" << s << "    CurrentPressure" << s;
        cycleInfoFile_ << std::endl;
        
        cycleInfoFile_.close();
    }
}

void CBReferenceRecovery::CalcPressures(TFloat time) {
    TFloat pressureFactor = (time - inflationStartTime_) / (inflationEndTime_ - inflationStartTime_);
    
    std::string incrementProgress = "";
    
    if (numIncrements_ != 1)
        incrementProgress = " of increment " + std::to_string(S_.incrementCount_) + "/" + std::to_string(numIncrements_);
    DCCtrl::print << " - " << 100*pressureFactor << "%" << incrementProgress << std::endl;
    
    for (auto s : surfaces_)
        S_.currentPressures_.at(s) = pressureFactor * targetPressures_.at(s) / numIncrements_ * S_.incrementCount_;
}

void CBReferenceRecovery::CalcResidual() {
    adapter_->GetSolver()->GetNodeCoordinates(S_.loadedCoords_);
    
    // residual = loadedCoords_ - targetCoords_
    VecCopy(S_.loadedCoords_, S_.residual_);
    VecAXPY(S_.residual_, -1, targetCoords_);
    
    VecAssemblyBegin(S_.residual_);
    VecAssemblyEnd(S_.residual_);
}

double CBReferenceRecovery::UpdateAugmentationParameter(TFloat augmentationParameter) {
    Vec residualDiff;
    PetscScalar a, b;
    
    VecDuplicate(S_.residual_, &residualDiff);
    VecCopy(S_.residual_, residualDiff);
    VecAXPY(residualDiff, -1, S_prev_.residual_);
    VecDot(S_prev_.residual_, residualDiff, &a);
    VecDot(residualDiff, residualDiff, &b);
    augmentationParameter *= -a/b;
    
    VecAssemblyBegin(residualDiff);
    VecAssemblyEnd(residualDiff);
    
    VecDestroy(&residualDiff);
    return augmentationParameter;
}

TFloat CBReferenceRecovery::CalcResidualNorm(NormType normType) {
    TFloat residualNorm;
    
    VecNorm(S_.residual_, normType, &residualNorm);
    return residualNorm;
}

void CBReferenceRecovery::ApplyBackwardDisplacement() {
    VecCopy(S_.unloadedCoords_, S_.prevUnloadedCoords_);
    
    // apply backward displacement: unloadedCoords_ = unloadedCoords_ - displacementFactor_ * residual_
    VecAXPY(S_.unloadedCoords_, -displacementFactor_, S_.residual_);
    
    VecAssemblyBegin(S_.unloadedCoords_);
    VecAssemblyEnd(S_.unloadedCoords_);
    
    adapter_->GetSolver()->SetNodeCoordinates(S_.unloadedCoords_);
    
    for (auto s : surfaces_)
        S_.unloadedVolumes_.at(s) = cavities_.at(s)->CalcVolume();
    
    adapter_->GetSolver()->RelaxElementsAndBases();
    adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
} // CBReferenceRecovery::ApplyBackwardDisplacement

void CBReferenceRecovery::ExportPressureVolumeInfo(TFloat time) {
    if (DCCtrl::IsProcessZero()) {
        pressureVolumeInfoFile_.open(pressureVolumeInfoFilename_.c_str(), std::ios::app);
        if (!pressureVolumeInfoFile_.good()) {
            throw std::runtime_error(
                                     "CBReferenceRecovery::ExportPressureVolumeInfo: Couldn't create " + pressureVolumeInfoFilename_ + ".");
        }
        pressureVolumeInfoFile_ << time;
        
        for (auto s : surfaces_)
            pressureVolumeInfoFile_ << "    " << S_.currentPressures_.at(s) << "    " << S_.currentVolumes_.at(s) * 1e6;
        
        pressureVolumeInfoFile_ << std::endl;
        
        pressureVolumeInfoFile_.close();
    }
}

void CBReferenceRecovery::ExportCycleInfo() {
    if (DCCtrl::IsProcessZero()) {
        cycleInfoFile_.open(cycleInfoFilename_.c_str(), std::ios::app);
        if (!cycleInfoFile_.good())
            throw std::runtime_error("CBReferenceRecovery::ExportCycleInfo: Couldn't create " + cycleInfoFilename_ + ".");
        cycleInfoFile_ << S_.cycleCount_ << "    " << S_.incrementCount_ << "    " << S_.cycleFinishedTime_ << "    " <<
        S_.residualNorm_ << "    " << S_.residualNormL2_;
        
        for (auto s : surfaces_)
            cycleInfoFile_ << "    " << S_.unloadedVolumes_.at(s) * 1e6 << "    " << S_.currentPressures_.at(s);
        
        cycleInfoFile_ << std::endl;
        
        cycleInfoFile_.close();
    }
}

void CBReferenceRecovery::ExportCoordsAsNodeFile(Vec coords, bool isFinalCycle) {
    if (DCCtrl::IsParallel()) {
        // Gather to process zero
        VecScatterBegin(coordsScatter_, coords, coordsSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(coordsScatter_, coords, coordsSeq_, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    if (DCCtrl::IsProcessZero()) {
        TInt numNodes  = adapter_->GetSolver()->GetNumberOfNodes();
        TInt numCoords = 3*numNodes;
        
        ////// Get node coords //////
        
        std::vector<PetscInt> pos(numCoords);
        for (PetscInt i = 0; i < numCoords; i++)
            pos[i] = i;
        
        std::vector<PetscScalar> values(numCoords);
        if (DCCtrl::IsParallel())
            VecGetValues(coordsSeq_, numCoords, pos.data(), values.data());
        else
            VecGetValues(coords, numCoords, pos.data(), values.data());
        
        ////// Get boundary conditions //////
        
        bool *bc = new bool[numCoords];
        adapter_->GetNodesComponentsBoundaryConditionsGlobal(numCoords, pos.data(), bc);
        
        ////// Write node file //////
        
        std::string nodeFilename = exportDir_;
        if (isFinalCycle)
            nodeFilename += "/InflatedState_Incr" + std::to_string(
                                                                   S_.incrementCount_) + ".node";
        else
            nodeFilename += "/UnloadedState_Incr" + std::to_string(
                                                                   S_.incrementCount_) + ".node";
        
        nodeFile_.open(nodeFilename.c_str());
        if (!nodeFile_.good())
            throw std::runtime_error("CBReferenceRecovery::ExportCoordsAsNodeFile: Couldn't create " + nodeFilename + ".");
        
        nodeFile_ << numNodes << " 3 1 0" << std::endl;
        
        for (PetscInt i = 0; i < numNodes; i++) {
            TFloat x = 1e3*values[3*i];
            TFloat y = 1e3*values[3*i+1];
            TFloat z = 1e3*values[3*i+2];
            TInt   b = bc[3*i] + (bc[3*i+1] << 1) + (bc[3*i+2] << 2); // 001 (1): x fixed; 010 (2): y fixed; 100 (4): z fixed;
                                                                      // 111 (7): x,y,z fixed
            nodeFile_ << i+1 << " " << x << " " << y << " " << z << " " << b << std::endl;
        }
        
        delete[] bc;
        nodeFile_.close();
    }
} // CBReferenceRecovery::ExportCoordsAsNodeFile

/// blatant copy of CBSolver, we are only interested in bases after ApplyBackwardsDisplacement()
void CBReferenceRecovery::ExportFibersAsBasesFile(bool isFinalCycle) {
    std::vector<CBElementSolid *> solidElements = GetSolver()->GetSolidElementVector();
    TInt nQP = solidElements.at(0)->GetNumberOfQuadraturePoints();
    
    for (auto &e : solidElements) {
        if (e->GetNumberOfQuadraturePoints() != nQP)
            throw std::runtime_error("ERROR: Tetgen export is currently only implemented for all T4 or all T10 elements.");
    }
    
    /// gather to zero first...
    Vec fVec, sVec, snVec;
    
    if (DCCtrl::IsParallel())
        VecCreateMPI(Petsc::Comm(), 3 * solidElements.size() * nQP, PETSC_DETERMINE, &fVec);
    else
        VecCreateSeq(PETSC_COMM_SELF, 3 * solidElements.size() * nQP, &fVec);
    
    VecDuplicate(fVec, &sVec);
    VecDuplicate(fVec, &snVec);
    
    VecZeroEntries(fVec);
    VecZeroEntries(sVec);
    VecZeroEntries(snVec);
    
    for (auto &it : solidElements) {
        for (TInt q = 0; q < nQP; q++) {
            Vector3<TFloat> f  = it->GetBasisAtQuadraturePoint(q)->GetCol(0);
            Vector3<TFloat> s  = it->GetBasisAtQuadraturePoint(q)->GetCol(1);
            Vector3<TFloat> sn = it->GetBasisAtQuadraturePoint(q)->GetCol(2);
            
            // rescale f to length = 1
            f.Normalize();
            
            // Gram Schmidt orthogonalize s wrt. f
            s.Normalize();
            s = s - (f*s) * f;
            
            // sn needs to be orthogonal to the others
            sn.Normalize();
            sn = CrossProduct(f, s);
            
            int idx = it->GetIndex();
            PetscInt index1[3] = {nQP *3*idx + 3*q, nQP *3*idx + 3*q + 1, nQP *3*idx + 3*q + 2};
            VecSetValues(fVec, 3, index1, f.GetArray(), INSERT_VALUES);
            PetscInt index2[3] = {nQP *3*idx + 3*q, nQP *3*idx + 3*q + 1, nQP *3*idx + 3*q + 2};
            VecSetValues(sVec, 3, index2, s.GetArray(), INSERT_VALUES);
            PetscInt index3[3] = {nQP *3*idx + 3*q, nQP *3*idx + 3*q + 1, nQP *3*idx + 3*q + 2};
            VecSetValues(snVec, 3, index3, sn.GetArray(), INSERT_VALUES);
        } // end nQP
    } // end solidElements
    
    VecAssemblyBegin(fVec);  VecAssemblyEnd(fVec);
    VecAssemblyBegin(sVec);  VecAssemblyEnd(sVec);
    VecAssemblyBegin(snVec); VecAssemblyEnd(snVec);
    
    VecScatter fScatter, sScatter, nScatter;
    Vec fLocal, sLocal, nLocal;
    
    VecScatterCreateToZero(fVec, &fScatter, &fLocal);
    VecScatterCreateToZero(sVec, &sScatter, &sLocal);
    VecScatterCreateToZero(snVec, &nScatter, &nLocal);
    
    VecScatterBegin(fScatter, fVec, fLocal, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(fScatter, fVec, fLocal, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBegin(sScatter, sVec, sLocal, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(sScatter, sVec, sLocal, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBegin(nScatter, snVec, nLocal, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(nScatter, snVec, nLocal, INSERT_VALUES, SCATTER_FORWARD);
    
    PetscScalar *pf, *ps, *pn;
    
    VecGetArray(fLocal, &pf);
    VecGetArray(sLocal, &ps);
    VecGetArray(nLocal, &pn);
    
    if (DCCtrl::IsProcessZero()) {
        /// Write bases file
        std::string basesFilename = exportDir_;
        if (isFinalCycle)
            basesFilename += "/UnloadedState_Final.bases";
        else
            basesFilename += "/UnloadedState_Incr" + std::to_string(S_.incrementCount_) + ".bases";
        
        basesFile_.open(basesFilename);
        if (!basesFile_.good())
            throw std::runtime_error("CBModelExporterTetgen::ExportModel: Couldn't create " + basesFilename + ".");
        
        size_t numElements = GetSolver()->GetModel()->GetSolidElements().size();
        basesFile_ << numElements << " " << nQP << std::endl; // HEADER
        
        for (TInt e = 0; e < numElements; e++) {
            basesFile_ << e+1;
            for (TInt i = 0; i < nQP; i++) {
                Vector3<TFloat> f  = {pf[nQP*3*e + 3*i], pf[nQP*3*e + 3*i + 1], pf[nQP*3*e + 3*i + 2]};
                Vector3<TFloat> s  = {ps[nQP*3*e + 3*i], ps[nQP*3*e + 3*i + 1], ps[nQP*3*e + 3*i + 2]};
                Vector3<TFloat> sn = {pn[nQP*3*e + 3*i], pn[nQP*3*e + 3*i + 1], pn[nQP*3*e + 3*i + 2]};
                
                // rescale f to length = 1
                f.Normalize();
                
                // Gram Schmidt orthogonalize s wrt. f
                s.Normalize();
                s = s - (f*s) * f;
                
                // sn needs to be orthogonal to the others
                sn.Normalize();
                sn = CrossProduct(f, s);
                
                basesFile_ << " ";
                basesFile_ << " " <<  f(0) << " " <<  f(1) << " " <<  f(2);
                basesFile_ << " " <<  s(0) << " " <<  s(1) << " " <<  s(2);
                basesFile_ << " " << sn(0) << " " << sn(1) << " " << sn(2);
            }
            basesFile_ << std::endl;
        }
        basesFile_.close();
    }
    
    VecRestoreArray(fLocal, &pf);
    VecRestoreArray(sLocal, &ps);
    VecRestoreArray(nLocal, &pn);
    
    VecScatterDestroy(&fScatter);
    VecDestroy(&fLocal);
    VecScatterDestroy(&sScatter);
    VecDestroy(&sLocal);
    VecScatterDestroy(&nScatter);
    VecDestroy(&nLocal);
    
    VecDestroy(&fVec);
    VecDestroy(&sVec);
    VecDestroy(&snVec);
} // CBReferenceRecovery::ExportFibersAsBasesFile

void CBReferenceRecovery::SaveStateVariablesAsNeeded(TFloat time) {
    // time is bigger than in the last call? -> the last call was successfull, store S_curr to S_prev
    // time is smaller than in last call? -> previous run was non-successfull, discard last computed values S_curr
    // time is roughly the same? -> update the current guess S_curr (with what?) ...
    if (time - S_curr_.time_ > 0) {
        S_prev_ = S_curr_;
    } else if (time - S_prev_.time_ < 0) {
        S_curr_ = S_prev_;
    } else {
        S_curr_ = S_prev_;
    }
    
    // reset state variables to the last known correct state
    S_ = S_prev_;
}
