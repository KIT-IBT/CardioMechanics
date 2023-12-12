/*
 * File: CBCircModel.h
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


#ifndef CB_CIRC_MODEL_H
#define CB_CIRC_MODEL_H

#include "Matrix3.h"

#include "CBElementCavity.h"
#include "ParameterMap.h"
#include "DCType.h"

using namespace math_pack;

// Base class needed as common type to store pointers to instances of class template CBCircModel<nC,nS,nR> in a std::vector<CBCircModelBase*>
class CBCircModelBase
{
public:
    virtual ~CBCircModelBase() { };
    
    virtual void InitSteadyStateCheck(ParameterMap* parameters, std::string parameterPrefix, TInt iSVD) = 0;
    virtual bool SteadyStateCheck(TFloat time) = 0;
    
    virtual TInt GetNumberOfCavities() = 0;
    virtual void GetCavitySurfaceIndices(TInt cavitySurfaceIndices[]) = 0;
    virtual void GetPreloadingPressures(TFloat preloadingPressures[]) = 0;
    
    virtual void SetInitialCavityVolumes(TFloat cavityVolumes[]) = 0;
    virtual void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]) = 0;
    
    virtual void AcceptEstState() = 0;
    virtual void RevertPrevState() = 0;
    
    virtual void SetFilename(std::string filename) = 0;
    virtual void WriteHeaderToFile() = 0;
    virtual void WriteResultsToFile(TFloat time) = 0;
};


template <TInt nC, TInt nS, TInt nR>
class CBCircModel : public CBCircModelBase
{
public:
    CBCircModel(ParameterMap* parameters, std::string parameterPrefix);
    virtual ~CBCircModel() { };
    
    void InitSteadyStateCheck(ParameterMap* parameters, std::string parameterPrefix, TInt iSVD);
    bool SteadyStateCheck(TFloat time);
    
    TInt GetNumberOfCavities() { return nC_; }
    void GetCavitySurfaceIndices(TInt cavitySurfaceIndices[]);
    void GetPreloadingPressures(TFloat preloadingPressures[]);
    
    virtual void SetInitialCavityVolumes(TFloat cavityVolumes[]) = 0;
    virtual void GetEstCavityVolumes(TFloat cavityPressures[], TFloat timeStep, TFloat cavityVolumes[]) = 0;
    
    void AcceptEstState();
    void RevertPrevState();
    
    void SetFilename(std::string filename) { filename_ = filename; }
    void WriteHeaderToFile();
    void WriteResultsToFile(TFloat time);
    
protected:
    std::string name_;
    
    static const TInt nC_ = nC;
    static const TInt nS_ = nS;
    static const TInt nR_ = nR;
    
    //----- State Variables -----
    
    TFloat stateVars_[nS_];
    TFloat prevStateVars_[nS_];
    TFloat estStateVars_[nS_];
    std::string stateVarsNames_[nS_];
    
    //----- Results -----
    
    TFloat* results_[nR_];
    std::string resultsNames_[nR_];
    
    //----- Cavity Pressures and Indices -----
    
    TFloat cavitySurfaceIndices_[nC_];
    TFloat cavityPreloadingPressures_[2*nC_];
    
    TFloat cavityPressures_[nC_] = {0}; // needs to be initialized
    TFloat prevCavityPressures_[nC_] = {0};
    TFloat* estCavityPressures_[nC_];
    
    //----- Steady State Check -----
    
    TInt iSVD_; // index of the stroke volume difference in stateVars_
    bool steadyChkActive_, steadyChkTemporal_ = false;
    TFloat steadyChkStart_, steadyChkPeriod_, steadyChkThresh_;
    TFloat *steadyChkParam_, prevSteadyChkParam_ = 0;
    TFloat prevSteadyChkTime_ = 0, prevSteadyChkTimeActive_ = 0;
    
    //----------
    
    std::ofstream file_;
    std::string filename_;
    TFloat maxTimeStep_;
    
    //----------
    
    void Integrate(TFloat cavityPressures[], TFloat timeStep);
    
private:
    virtual void Algebraics(TFloat cavityPressures[], TFloat stateVars[]) = 0;
    virtual void Rates(TFloat cavityPressures[], TFloat stateVars[], TFloat stateVarsDot[]) = 0;
};

// Methods of class templates must be implemented in the header
#include "CBCircModel.hpp"

#endif
