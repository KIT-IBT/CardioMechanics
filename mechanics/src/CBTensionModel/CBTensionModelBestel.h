/*
 * File: CBTensionModelBestel.h
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


#pragma once

#include "CBTensionModel.h"
#include "ParameterMap.h"

/// Tension model based on three ordinary differential equations. Detail can be found in the publication "Three-Wall Segment (TriSeg) Model Describing Mechanics and Hemodynamics of Ventricular Interaction" by Joost Lumens 2009.
class CBTensionModelBestel : public CBTensionModel {
protected:
    /// flags to control the exporting of the state variables
    bool doExport_ = false;
    bool headerWritten_ = false;
    std::string filename_ = "";
    TInt ei_ = -1;
    
    // cycle length
    TFloat cycleLength_ = 1.0;
    
    // start time of activation
    TFloat startTime_ = 0.0;
    
    /// CONSTANTS / PARAMETER VALUES
    TFloat contractility_;
    TFloat maxActivationRate_;
    TFloat minActivationRate_;
    TFloat steepness_;
    TFloat onsetSystole_;
    TFloat onsetDiastole_;
    
    /// initialize all constant / parameter values from a given xml file, fallback  is used to evaluate "Mat_Default" cases
    void InitParamsFromXml(ParameterMap *params, std::string parameterKey, std::string parameterKeyFallback);
    
    /// state variable container, these exist multiple times and are going to be tracked for step_back
    struct StateVariables {
        TFloat t;
        TFloat a;
        TFloat f;
        TFloat stress;
    };
    
    /// state variables from last successful guess
    StateVariables S_prev_;
    
    /// last guess (might have been non-successful)
    StateVariables S_curr_;
    
    /// temporary variable actually used to store the freshly computed guess
    StateVariables S_;
    
    TFloat kff_ = 1.0;
    TFloat kss_ = 0.0;
    TFloat knn_ = 0.0;
    TFloat ksn_ = 0.0;
    
    // FUNCTIONS needed to compute the state variables
    /// time-dependent stress function
    TFloat dTaudt(TFloat tau, TFloat a, TFloat contractility);
    TFloat a(TFloat f);
    TFloat f(TFloat Sp, TFloat Sm);
    
public:
    CBTensionModelBestel(CBElement *e, ParameterMap *params);
    
    ~CBTensionModelBestel() {}
    
    /// Care for the interrelation between S_prev, S_curr and a temporary S
    /// do not save step-back results, discard timesteps that are too small
    void SaveStateVariablesAsNeeded(TFloat time);
    
    /// write state variables for specific indices to a file after each successful step
    void WriteToFile(const StateVariables &S);
    
    /// active tension computation function of the model.
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override;
};
