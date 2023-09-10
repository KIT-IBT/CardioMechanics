/*
 *  CBTensionModelTanH.h
 *  CardioMechanics
 *
 *  Created by Tobias Gerach, Fr. 08. May 2021
 *
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#pragma once

#include "CBTensionModel.h"
#include "ParameterMap.h"

class CBTensionModelTanH : public CBTensionModel {
protected:
    /// flags to control the exporting of the state variables
    bool doExport_        = false;
    bool headerWritten_   = false;
    std::string filename_ = "";
    TInt ei_              = -1; // element index
    
    // cycle length
    TFloat cycleLength_ = 1.0;
    
    // start time of activation
    TFloat startTime_ = 0.0;
    
    /// CONSTANTS / PARAMETER VALUES
    // minimum fiber stretch
    TFloat l_0_ = 0.7;
    
    // electromechanical delay [s]
    TFloat t_emd_ = 0.015;
    
    // peak isometric tension [kPa]
    TFloat S_peak_ = 120.0;
    
    // duration of active contraction [s]
    TFloat t_dur_ = 0.3;
    
    // baseline time constant of active contraction [s]
    TFloat tau_c0_ = 0.1;
    
    // degree of length dependence
    TFloat ld_ = 5.0;
    
    // length dependence of upstroke time [s]
    TFloat ld_up_ = 0.5;
    
    // time constant of relaxation [s]
    TFloat tau_r_ = 0.1;
    
    /// initialize all constant / parameter values from a given xml file, fallback  is used to evaluate "Mat_Default"
    /// cases
    void InitParamsFromXml(ParameterMap *params, std::string parameterKey, std::string parameterKeyFallback);
    
    /// state variable container, these exist multiple times and are going to be tracked for step_back
    struct StateVariables {
        TFloat t;
        TFloat lambda;
        TFloat Ta;
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
    
public:
    CBTensionModelTanH(CBElementSolid *e, ParameterMap *params);
    
    ~CBTensionModelTanH() {}
    
    /// Care for the interrelation between S_prev, S_curr and a temporary S
    /// do not save step-back results, discard timesteps that are too small
    void SaveStateVariablesAsNeeded(TFloat time);
    
    /// write state variables for specific indices to a file after each successful step
    void WriteToFile(const StateVariables &S);
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override;
};
