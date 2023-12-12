/*
 * File: CBTensionModelLand17.h
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

class CBTensionModelLand17 : public CBTensionModel {
protected:
    /// solid element
    CBElement *e_;
    
    /// flags to control the exporting of the state variables
    bool doExport_ = false;
    bool headerWritten_ = false;
    std::string filename_ = "";
    TInt ei_ = -1; // element index
    
    // fixed cycle length for CaCoppini and CaBrixius
    TFloat cycleLength_ = 1.0;
    
    // rate dependance on/off
    std::string rateDependancy_ = "";
    
    // start time of activation
    TFloat startTime_ = 0.0;
    
    // decision flag for atrial or ventricular Ca
    std::string calciumTransientType_ = "";
    
    // external Cai
    std::vector<double> humanCai_;
    std::string calciumFile_ = "";
    std::vector<TFloat> caiQP_;
    
    // CONSTANTS / PARAMETER VALUES
    // cooperativity of the calcium-troponin C binding rate
    TFloat TRPNn_ = 2.0;
    
    // troponin unbinding rate in 1/ms
    TFloat TRPNk_ = 0.1;
    
    // half activation point in uM
    TFloat Ca50_ = 0.805;
    
    // tropomyosin rate constant in 1/ms
    TFloat ku_ = 1.0;
    
    // Hill coefficient
    TFloat nTM_ = 5.0;
    
    // Value of CaTRPN where the blocked tropomyosin binding sites on actin TmBlocked = 0.5 in steady state
    TFloat TRPN50_ = 0.35;
    
    // unbound-to-weak transition rate kuw = 0.026 * nu in 1/ms
    TFloat kuw_ = 0.182;
    
    // weak-to-strong crossbridge transition rate kws = 0.004 * mu in 1/ms
    TFloat kws_ = 0.012;
    
    // steady-state duty ratio
    TFloat rs_ = 0.25;
    
    // steady-state ratio between pre-powerstroke and non-strongly bound
    TFloat rw_ = 0.5;
    
    // distortion-dependent unbinding rate factor of S gate
    TFloat gammas_ = 0.0085;
    
    // distortion-dependent unbinding rate factor of W gate
    TFloat gammaw_ = 0.615;
    
    // combination of (cw,cs). relate to the decay rate of distortion
    TFloat phi_ = 2.23;
    
    // rescaling of the distortion induced by relative movement of the filaments
    TFloat Aeff_ = 25.0;
    
    // change in maximal tension based on changes in filament overlap
    TFloat beta0_ = 2.3;
    
    // change in calcium sensitivity
    TFloat beta1_ = -2.4;
    
    // maximal active tension at resting length in kPa
    TFloat Tref_ = 80.0;
    
    // rescaling parameter for total passive force in kPa
    TFloat a_ = 2.1;
    
    // serial spring stiffness
    TFloat k_ = 7.0;
    
    // viscocity of linear dashpot element in 1/ms
    TFloat eta_l_ = 200.0;
    TFloat eta_s_ = 20.0;
    
    // scaling factor for all crossbridge cycling rates
    TFloat xi_ = 1.0;
    
    // constants that will be calculated once during initialize
    TFloat cw_, cs_, kwu_, ksu_, A_, XSSS_, XWSS_, ktm_block_;
    
    /// initialize all constant / parameter values from a given xml file, fallback  is used to evaluate "Mat_Default" cases
    void InitParamsFromXml(ParameterMap *params, std::string parameterKey, std::string parameterKeyFallback);
    
    /// state variable container, these exist multiple times and are going to be tracked for step_back
    struct StateVariables {
        TFloat t;
        TFloat delta_t;
        TFloat XS;
        TFloat XW;
        TFloat TRPN;
        TFloat TmBlocked;
        TFloat ZetaS;
        TFloat ZetaW;
        TFloat Cd;
        TFloat lambda;
        TFloat dlambdadt;
        TFloat Ta;
        TFloat Tension;
    };
    
    /// functions needed for state variable calculation
    TFloat dXSdt(TFloat xb_ws, TFloat xb_su, TFloat xb_su_gamma);
    TFloat dXWdt(TFloat xb_uw, TFloat xb_wu, TFloat xb_ws, TFloat xb_wu_gamma);
    TFloat dTRPNdt(TFloat TRPNk, TFloat Cai, TFloat Ca50, TFloat TRPNn, TFloat TRPN);
    TFloat dTmBlockeddt(TFloat ktm_block, TFloat TRPN_NP, TFloat XU, TFloat ku, TFloat TRPN, TFloat TRPNn,
                        TFloat TmBlocked);
    TFloat dZetaSdt(TFloat A, TFloat dlambdadt, TFloat cs, TFloat ZetaS);
    TFloat dZetaWdt(TFloat A, TFloat dlambdadt, TFloat cw, TFloat ZetaW);
    TFloat dC_ddt(TFloat k, TFloat eta, TFloat C_s);
    TFloat CaBrixius(TFloat t);
    TFloat CaCoppini(TFloat t);
    TFloat CaExternal(TFloat t);
    TFloat CaElphy(TFloat t);
    TFloat Cai(std::string CaDecisionFlag, TFloat t);
    TFloat Overlap(TFloat lambda);
    
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
    CBTensionModelLand17(CBElementSolid *e, ParameterMap *params);
    
    ~CBTensionModelLand17() {}
    
    /// Care for the interrelation between S_prev, S_curr and a temporary S
    /// do not save step-back results, discard timesteps that are too small
    void SaveStateVariablesAsNeeded(TFloat time);
    
    /// write state variables for specific indices to a file after each successful step
    void WriteToFile(const StateVariables &S);
    void ReadExternalCai(std::string filename);
    
    /// Even though the name suggest this function sets active tension, it sets calcium at a QP
    /// This way, we can use both staggered coupling schemes without additional code
    CBStatus SetActiveTensionAtQuadraturePoint(TInt, TFloat) override;
    
    /// active tension computation function of the  model. Performs an explicit euler step for the ode parts and tracks active step backs.
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override;
    
    /// active stiffness required for stabilization
    TFloat CalcActiveStiffness();
};
