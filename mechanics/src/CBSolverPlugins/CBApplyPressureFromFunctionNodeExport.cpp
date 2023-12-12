/*
 * File: CBApplyPressureFromFunctionNodeExport.cpp
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

#include "CBApplyPressureFromFunctionNodeExport.h"
#include "CBSolver.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"

void CBApplyPressureFromFunctionNodeExport::Init() {
    // init filenames
    std::string defaultFilename = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir() +
    "/ApplyPressureFromFunctionNodeExport.dat";
    
    filename_ = parameters_->Get<std::string>("Plugins.ApplyPressureFromFunctionNodeExport.ExportFile", defaultFilename);
    
    // init cavities map using function from base class
    CBSolverPluginCavities::InitCavities("Plugins.ApplyPressureFromFunctionNodeExport");
    
    // init groups and intervals
    std::vector<std::string> groupKeys = parameters_->GetChildNodes("Plugins.ApplyPressureFromFunctionNodeExport.Groups");
    
    for (auto &groupKey : groupKeys) {
        std::vector<TInt> surfaces = parameters_->GetArray<TInt>(groupKey + ".Surfaces");
        
        std::vector<std::string> intervalKeys = parameters_->GetChildNodes(groupKey + ".Intervals");
        if (intervalKeys.size() == 0)
            throw std::runtime_error(
                                     "void CBApplyPressureFromFunctionNodeExport::Init(): No Intervals specified in " + groupKey);
        
        TFloat prevStopTime = 0.0;
        std::vector<IntervalStruct> intervals;
        for (auto &intervalKey : intervalKeys) {
            IntervalStruct interval;
            
            interval.function  = std::make_shared<CBDataFromFunction>(parameters_, intervalKey);
            interval.relaxElementsAtStart = parameters_->Get<bool>(intervalKey + ".RelaxElementsAtStart", false);
            if (interval.relaxElementsAtStart)
                interval.materialsToRelax = parameters_->GetArray<TInt>(intervalKey + ".MaterialsToRelax");
            interval.invert    = parameters_->Get<bool>(intervalKey + ".Invert", false);
            interval.startTime = parameters_->Get<TFloat>(intervalKey + ".StartTime", 0.0);
            if (interval.startTime < prevStopTime)
                throw std::runtime_error(
                                         "void CBApplyPressureFromFunctionNodeExport::Init(): StartTime of " + intervalKey +
                                         " is smaller than StopTime of previous Interval");
            interval.stopTime  = parameters_->Get<TFloat>(intervalKey + ".StopTime", 0.0);
            prevStopTime = interval.stopTime;
            interval.offset    = parameters_->Get<TFloat>(intervalKey + ".Offset", 0.0);
            interval.amplitude = parameters_->Get<TFloat>(intervalKey + ".Amplitude", 0.0);
            
            intervals.push_back(interval);
        }
        
        for (auto &surface : surfaces) {
            if (cavities_.find(surface) == cavities_.end())
                throw std::runtime_error("void CBApplyPressureFromFunctionNodeExport::Init(): No cavity with surface index " + std::to_string(
                                                                                                                                              surface) + " found.");
            
            groups_.insert(std::pair<TInt, std::vector<IntervalStruct>>(surface, intervals));
            ValueStruct vs;
            vs.pressure = 0.0;
            vs.volume = 0.0;
            valuesStructs_.insert(std::pair<TInt, ValueStruct>(surface, vs));
        }
    }
    
    stopSurface_ = parameters_->Get<TInt>("Plugins.ApplyPressureFromFunctionNodeExport.StopSurface", -1);
    if (stopSurface_ != -1)
        stopVolume_ = parameters_->Get<TFloat>("Plugins.ApplyPressureFromFunctionNodeExport.StopVolume");
    else
        nodeExportTime_ = parameters_->Get<TFloat>("Plugins.ApplyPressureFromFunctionNodeExport.NodeExportTime");
    nodeFilename_ = parameters_->Get<std::string>("Plugins.ApplyPressureFromFunctionNodeExport.NodeExportFile");
    InitPetscVectors();
} // CBApplyPressureFromFunctionNodeExport::Init

void CBApplyPressureFromFunctionNodeExport::Apply(TFloat time) {
    currentTime_ = time;
    auto &timing = adapter_->GetSolver()->GetTiming();
    for (groupsIt_ = groups_.begin(); groupsIt_ != groups_.end(); groupsIt_++) {
        for (auto &interval : groupsIt_->second) {
            TFloat dt = timing.GetTimeStep();
            if ((time >= interval.startTime) && (time < interval.stopTime+0.5*dt) ) {
                if (interval.relaxElementsAtStart) {
                    for (auto &mat : interval.materialsToRelax) {
                        materialsRelaxedIt_ = materialsRelaxed_.find(mat);
                        if ((materialsRelaxedIt_ == materialsRelaxed_.end()) ||
                            (materialsRelaxedIt_->second != interval.startTime) ) {
                            DCCtrl::print << "Relaxing elements with material " << mat << std::endl;
                            Base::adapter_->GetSolver()->RelaxElementsAndBasesByMaterialIndex(mat);
                            materialsRelaxed_.insert(std::pair<TInt, TFloat>(mat, interval.startTime));
                        }
                    }
                    interval.relaxElementsAtStart = false;
                }
                
                TFloat f = interval.function->Get(time, 0);
                if (interval.invert)
                    valuesStructs_[groupsIt_->first].pressure = interval.offset + interval.amplitude * (1.0-f);
                else
                    valuesStructs_[groupsIt_->first].pressure = interval.offset + interval.amplitude * f;
                
                break;
            }
        }
    }
} // CBApplyPressureFromFunctionNodeExport::Apply

void CBApplyPressureFromFunctionNodeExport::ApplyToNodalForces() {
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
        cavities_[valuesIt_->first]->ApplyPressure(valuesIt_->second.pressure);
}

void CBApplyPressureFromFunctionNodeExport::ApplyToNodalForcesJacobian() {
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
        cavities_[valuesIt_->first]->ApplyPressureToJacobian(valuesIt_->second.pressure);
}

void CBApplyPressureFromFunctionNodeExport::AnalyzeResults() {
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++) {
        valuesIt_->second.volume = cavities_[valuesIt_->first]->CalcVolume();
        
        if ( ((stopSurface_ == -1) && (std::abs(currentTime_-nodeExportTime_) < 1e-5) ) ||
            ((valuesIt_->first == stopSurface_) && (valuesIt_->second.volume*1e6 <= stopVolume_) ) ) {
            ExportNodeFile();
            shallExit_ = true;
        }
    }
}

void CBApplyPressureFromFunctionNodeExport::WriteHeaderToFile() {
    file_.open(filename_.c_str());
    if (!file_.good())
        throw std::runtime_error(
                                 "void CBApplyPressureFromFunctionNodeExport::WriteToFile(TFloat time): Couldn't create " + filename_);
    file_ << std::setw(16) << "time";
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++) {
        std::ostringstream p, v;
        p << "pressure" << valuesIt_->first;
        v << "volume" << valuesIt_->first;
        file_ << std::setw(16) << p.str() << std::setw(16) << v.str();
    }
    file_ << std::endl;
    file_.close();
}

void CBApplyPressureFromFunctionNodeExport::WriteToFile(TFloat time) {
    // write header on first call
    if (!headerWritten_)
        WriteHeaderToFile();
    headerWritten_ = true;
    
    // write plugin data
    file_.open(filename_.c_str(), std::ios::app);
    if (!file_.good())
        throw std::runtime_error(
                                 "void CBApplyPressureFromFunctionNodeExport::WriteToFile(TFloat time): Couldn't create " + filename_);
    file_ << std::setprecision(7) << std::fixed << std::setw(16) << time;
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++) {
        file_ << std::setw(16) << valuesIt_->second.pressure/133.322 << std::setw(16) << (valuesIt_->second).volume*1e6;
        std::cout << std::setw(16) << valuesIt_->second.pressure/133.322 << std::setw(16) << (valuesIt_->second).volume*1e6;
    }
    std::cout << std::endl;
    file_ << std::endl;
    file_.close();
}

//////////////////// Functions needed for nodes export ////////////////////

void CBApplyPressureFromFunctionNodeExport::InitPetscVectors() {
    DCCtrlPETSc::CreateVector(3*adapter_->GetSolver()->GetNumberOfLocalNodes(), PETSC_DETERMINE, &coords_);
    
    if (DCCtrl::IsParallel()) {
        VecCreateSeq(PETSC_COMM_SELF, 3*adapter_->GetSolver()->GetNumberOfNodes(), &coordsSeq_);
        VecScatterCreateToZero(coords_, &coordsScatter_, &coordsSeq_);
    }
}

void CBApplyPressureFromFunctionNodeExport::ExportNodeFile() {
    adapter_->GetSolver()->GetNodeCoordinates(coords_);
    
    if (DCCtrl::IsParallel()) {
        // Gather to process zero
        VecScatterBegin(coordsScatter_, coords_, coordsSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(coordsScatter_, coords_, coordsSeq_, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    if (DCCtrl::IsProcessZero()) {
        TInt numNodes = adapter_->GetSolver()->GetNumberOfNodes();
        TInt numCoords = 3*numNodes;
        
        ////// Get node coords //////
        
        std::vector<PetscInt> pos(numCoords);
        for (PetscInt i = 0; i < numCoords; i++)
            pos[i] = i;
        
        std::vector<PetscScalar> values(numCoords);
        if (DCCtrl::IsParallel())
            VecGetValues(coordsSeq_, numCoords, pos.data(), values.data());
        else
            VecGetValues(coords_, numCoords, pos.data(), values.data());
        
        ////// Get boundary conditions //////
        
        bool *bc = new bool[numCoords];
        adapter_->GetNodesComponentsBoundaryConditionsGlobal(numCoords, pos.data(), bc);
        
        ////// Write node file //////
        
        nodeFile_.open(nodeFilename_.c_str());
        if (!nodeFile_.good())
            throw std::runtime_error(
                                     "CBApplyPressureFromFunctionNodeExport::ExportNodeFile: Couldn't create " + nodeFilename_ + ".");
        
        nodeFile_ << numNodes << " 3 1 0" << std::endl;
        
        for (PetscInt i = 0; i < numNodes; i++) {
            TFloat x = 1e3*values[3*i];
            TFloat y = 1e3*values[3*i+1];
            TFloat z = 1e3*values[3*i+2];
            TInt   b = bc[3*i] + (bc[3*i+1] << 1) + (bc[3*i+2] << 2); // 001 (1): x fixed; 010 (2): y fixed; 100 (4): z fixed; 111 (7): x,y,z fixed
            nodeFile_ << i+1 << " " << x << " " << y << " " << z << " " << b << std::endl;
        }
        
        delete[] bc;
        nodeFile_.close();
    }
} // CBApplyPressureFromFunctionNodeExport::ExportNodeFile
