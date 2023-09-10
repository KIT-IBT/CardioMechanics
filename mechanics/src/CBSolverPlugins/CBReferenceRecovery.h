/* -------------------------------------------------------
 
 CBReferenceRecovery.h
 
 Ver. 1.0.0
 
 Created:       Lukas Baron and Steffen Schuler      (01.2017)
 Last modified: Tobias Gerach  (01.02.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#pragma once

#ifndef CB_REFERENCE_RECOVERY
# define CB_REFERENCE_RECOVERY

# include <map>
# include <vector>
# include "Matrix3.h"
# include "Matrix4.h"
# include <list>

# include <petscsnes.h>

# include <unordered_map>
# include <unordered_set>

# include "CBSolverPlugin.h"

# include "CBElementSurface.h"
# include "CBCirculationCavity.h"
# include "CBModel.h"

# include "CBSolverPluginCavities.h"


using namespace math_pack;

class CBReferenceRecovery : public CBSolverPluginCavities {
public:
    // contructor and desctructor
    CBReferenceRecovery() {}
    
    virtual ~CBReferenceRecovery() {}
    
    // mandatory plugin functions -- basic structure
    void Init() override;
    void Apply(TFloat time) override;
    void ApplyToNodalForces() override;
    void ApplyToNodalForcesJacobian() override;
    void Export(TFloat time) override;
    
    // some other, mandatory functions
    void StepBack() override;
    
    bool WantsToAnalyzeResults() override {return true; }
    
    void AnalyzeResults() override;
    void WriteToFile(TFloat time) override;
    
    CBStatus GetStatus() override {return status_; }
    
    std::string GetName() override {return "ReferenceRecovery"; }
    
    // optional functions
    bool ExitCheck(TFloat time) override;
    
private:
    // functions that are needed by the public / main / top-level functions above
    // Initialization
    void   InitNodesAssociation();
    void   InitParamsFromXML();
    void   InitTargetCoordsFromCurrentCoords();
    void   InitExportPressureVolumeInfo();
    void   InitExportCycleInfo();
    void   InitPetscVec();
    void   CalcPressures(TFloat time);
    void   CalcResidual();
    TFloat CalcResidualNorm(NormType normType = NORM_INFINITY);
    void   ApplyBackwardDisplacement();
    void   ExportPressureVolumeInfo(TFloat time);
    void   ExportCycleInfo();
    void   ExportCoordsAsNodeFile(Vec coords, bool isFinalCycle = false);
    void   ExportFibersAsBasesFile(bool isFinalCycle = false);
    void   Finalize();
    double UpdateAugmentationParameter(TFloat augmentationParameter);
    void   SaveStateVariablesAsNeeded(TFloat time);
    
    std::string exportDir_;
    std::string pressureVolumeInfoFilename_;
    std::string cycleInfoFilename_;
    std::ofstream pressureVolumeInfoFile_;
    std::ofstream cycleInfoFile_;
    std::ofstream nodeFile_;
    std::ofstream basesFile_;
    
    TInt numIncrements_;
    
    TFloat startTime_;
    TFloat inflationStartTime_;
    TFloat inflationDuration_;
    TFloat inflationEndTime_;
    TInt   cycleMax_;
    
    /// Target precision between inflated and target mesh, that we want the Bols algorithm to reach
    TFloat tolerance_;
    
    /// factor to scale the negative displacement, 1.0 corresponds to original Bols, with highly nonlinear materials,
    /// factors below 1 can be beneficial
    TFloat displacementFactor_;
    
    /// node coordinates of target mesh that the Bols method ideally should reach
    Vec targetCoords_;
    
    /// switch to turn on augmentation
    bool augmentation_;
    
    /// state variables, organized as struct/class to facilitate StepBack
    class StateVariables {
    public:
        bool inflationFinished_   = false;
        TInt cycleCount_          = 0;
        TInt incrementCount_      = 1;
        TFloat time_              = 0;
        TFloat cycleFinishedTime_ = 0.0;
        TFloat residualNorm_ = INFINITY;
        TFloat residualNormL2_ = INFINITY;
        std::map<TInt, TFloat> unloadedVolumes_;
        std::map<TInt, TFloat> currentPressures_;              // < All pressures by surface that should be applied in the
                                                               // current time step.
        std::map<TInt, TFloat> currentVolumes_;
        Vec loadedCoords_;
        Vec unloadedCoords_;
        Vec prevLoadedCoords_;  // < loaded coordinates from previous cycle
        Vec prevUnloadedCoords_;  // < unloaded coordinates from previous cycle
        
        Vec residual_;
        
        StateVariables() {}
        
        void Init(const Vec &coords) {
            VecDuplicate(coords, &loadedCoords_);
            VecZeroEntries(loadedCoords_);
            VecDuplicate(coords, &unloadedCoords_);
            VecZeroEntries(unloadedCoords_);
            VecDuplicate(coords, &prevLoadedCoords_);
            VecZeroEntries(prevLoadedCoords_);
            VecDuplicate(coords, &prevUnloadedCoords_);
            VecZeroEntries(prevUnloadedCoords_);
            VecDuplicate(coords, &residual_);
            VecZeroEntries(residual_);
        }
        
        StateVariables operator=(const StateVariables &other) {
            inflationFinished_      = other.inflationFinished_;
            cycleCount_             = other.cycleCount_;
            incrementCount_         = other.incrementCount_;
            time_                   = other.time_;
            cycleFinishedTime_      = other.cycleFinishedTime_;
            residualNorm_           = other.residualNorm_;
            residualNormL2_         = other.residualNormL2_;
            unloadedVolumes_        = other.unloadedVolumes_;
            currentPressures_       = other.currentPressures_;
            currentVolumes_         = other.currentVolumes_;
            VecCopy(other.loadedCoords_, loadedCoords_);
            VecCopy(other.unloadedCoords_, unloadedCoords_);
            VecCopy(other.prevLoadedCoords_, prevLoadedCoords_);
            VecCopy(other.prevUnloadedCoords_, prevUnloadedCoords_);
            VecCopy(other.residual_, residual_);
            return *this;
        }
    }; // class StateVariables
    
    /// state variables from last successful guess
    StateVariables S_prev_;
    
    /// last guess (might have been non-successful)
    StateVariables S_curr_;
    
    /// temporary variable actually used to store the freshly computed guess
    StateVariables S_;
    
    
    /// needed to gather node coords to process zero for exporting
    Vec coordsSeq_;
    VecScatter coordsScatter_;
    
    /// All surfaces of cavities, where pressure should be applied to
    std::vector<TInt> surfaces_;
    
    std::map<TInt, TFloat> initialVolumes_;
    
    /// All pressures by surface that should be reached in the target mesh in its fully inflated configuration (loaded
    /// state in the Bols method).
    std::map<TInt, TFloat> targetPressures_;
    
    bool shallExit_ = false;
}; // class CBReferenceRecovery

#endif  // ifndef CB_BOLS_UNLOADING
