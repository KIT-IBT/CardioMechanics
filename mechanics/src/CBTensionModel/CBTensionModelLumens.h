/*
 * File: CBTensionModelLumens.h
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
class CBLumensTension : public CBTensionModel {
protected:
    
    /// flags to control the exporting of the state variables
    bool doExport_ = false;
    bool headerWritten_ = false;
    std::string filename_ = "";
    TInt ei_ = -1;
    
    // CONSTANTS / PARAMETER VALUES
    
    /// length of one heartbeat in s, 0.8s ^= 75bpm
    TFloat cycleLength_ = 0.8;
    /// start time of the active contraction in s
    //TFloat startTime_ = 0.0;
    
    /// diastolic resting level of activation, unit-free
    TFloat Crest_ = 0.0; // Lumens2009: 0.02, but we want a zero tension at rest
                         /// contractile element length with zero active stress, in meter
    TFloat Lsc0_ = 1.51 * 1e-6;
    /// reference sarcomere length at zero strain, in meter
    TFloat Lsref_ = 2.0 * 1e-6;
    /// length of isometrically stressed series elastic element, in meter
    TFloat Lseiso_ = 0.04 * 1e-6;
    /// sarcomere shortening velocity with zero load, in meter per second
    TFloat vmax_ = 7. * 1e-6;
    /// factor scaling contraction decay time, in s
    TFloat tauD_ = 32.0 * 1e-3;
    /// factor scaling contraction rise time, in s
    TFloat tauR_ = 48.0 * 1e-3;
    /// factor scaling duration of contraction, in s
    TFloat tausc_ = 425.0 * 1e-3;
    /// factor scaling active myofiber stress, in Pa
    TFloat sigma_act_ = 120.0 * 1e3;
    /// factor scaling passive myofiber stress, in Pa
    TFloat sigma_pas_ = 7.0 * 1e3;
    /// initialize all constant / parameter values from a given xml file, fallback  is used to evaluate "Mat_Default" cases
    void InitParamsFromXml(ParameterMap* params, std::string parameterKey, std::string parameterKeyFallback);
    
    /// state variable container, these exist multiple times and are going to be tracked for step_back
    struct StateVariables {
        TFloat t;
        TFloat eps_f;
        TFloat f_rise;
        TFloat ls;
        TFloat lsc;
        TFloat cl;
        TFloat T;
        TFloat c;
        TFloat sigmafact;
    };
    
    /// state variables from last successful guess
    StateVariables S_prev_;
    /// last guess (might have been non-successful)
    StateVariables S_curr_;
    /// temporary variable actually used to store the freshly computed guess
    StateVariables S_;
    
    // FUNCTIONS needed to compute the state variables
    
    /// sarcomere length from natural myofiber strain, Lumens: Lsref*exp(eps_f), Gurev: Lsref*(1+eps_f)
    TFloat Ls(TFloat eps_f);
    /// rise of mechanical activation (driving function for contraction)
    TFloat Frise(TFloat t);
    /// contractile element length (ode)
    TFloat dLscdt(TFloat Lsc, TFloat eps_f);
    /// increase of activation with sarcomere length
    TFloat CL(TFloat Lsc);
    /// decrease of activation duration with decrease of sarcomere length (as function, since a variable/quantity T is defined as well)
    TFloat T_fun(TFloat Lsc);
    /// mechanical activation (ode)
    TFloat dCdt(TFloat C, TFloat Lsc, TFloat eps_f, TFloat t);
    /// active myofiber stress
    TFloat sigma_fact(TFloat Lsc, TFloat C, TFloat eps_f);
    
public:
    CBLumensTension(CBElement* e, ParameterMap* params);
    
    ~CBLumensTension() { }
    
    /// Care for the interrelation between S_prev, S_curr and a temporary S
    /// do not save step-back results, discard timesteps that are too small
    void SaveStateVariablesAsNeeded(TFloat time);
    
    /// write state variables for specific indices to a file after each successful step
    void WriteToFile(const StateVariables& S);
    
    /// active tension computation function of the Lumens model. Performs an explicit euler step for the ode parts and tracks active step backs.
    virtual double CalcActiveTension(const math_pack::Matrix3<double>& deformation, const double time) override;
};
