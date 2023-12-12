/*
 * File: CBTensionModelBestel.cpp
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


// Reference: Bestel et al., 2001: https://link.springer.com/chapter/10.1007/3-540-45468-3_143

#include "CBTensionModelBestel.h"
#include <algorithm>
#include <cassert>

TFloat CBTensionModelBestel::dTaudt(TFloat tau, TFloat a, TFloat contractility) {
    /// stress
    TFloat x = -abs(a) * tau + contractility * std::max(a, 0.0);
    
    return x;
}

TFloat CBTensionModelBestel::a(TFloat f) {
    /// activation function
    return maxActivationRate_ * f + minActivationRate_ * (1 - f);
}

TFloat CBTensionModelBestel::f(TFloat S_p, TFloat S_m) {
    /// function f(t) in ]0,1[ indicates systole
    return S_p * S_m;
}

/// Read a value from xml with fallback. Priorities: key+tail > keyFallback+tail > defaultValue
static TFloat InitKey(ParameterMap *params, std::string key, std::string keyFallback, std::string tail,
                      TFloat defaultValue) {
    if ((params->IsAvailable(key + tail) == false) && (params->IsAvailable(keyFallback + tail) == true))
        return params->Get<TFloat>(keyFallback);
    else
        return params->Get<TFloat>(key+tail, defaultValue);
}

CBTensionModelBestel::CBTensionModelBestel(CBElement *e, ParameterMap *params) {
    // set model parameters
    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    if (Tmax_ > 1) {
        throw std::runtime_error(
                                 "CBTensionModelBestel: Parameter TensionMax > 1 is not compatible with this model to ensure model output in Pa!");
    }
    int mi = e->GetMaterialIndex();
    ei_ = e->GetIndex();
    std::string key = "Materials.Mat_" + std::to_string(mi);
    std::string keyFallback = "Materials.Mat_Default";
    InitParamsFromXml(params, key+".Bestel", keyFallback+".Bestel");
    
    // set initial values for the state variables
    S_.t = 0.;
    S_.a = 0.;
    S_.f = 0.;
    S_.stress = 0.;
    
    S_curr_ = S_;
    S_prev_ = S_;
}

void CBTensionModelBestel::SaveStateVariablesAsNeeded(TFloat time) {
    // time is bigger than in the last call? -> the last call was successful, store S_curr to S_prev
    // time is smaller than in last call? -> previous run was non-successful, discard last computed values S_curr
    // time is roughly the same? -> update the current guess S_curr (with what?) ...
    if (time - S_curr_.t > 0) {
        S_prev_ = S_curr_;
        
        // export state variables on a successful export
        if (doExport_)
            this->WriteToFile(S_prev_);
    } else if (time - S_prev_.t < 0) {
        S_curr_ = S_prev_;
    } else {
        S_curr_ = S_prev_;
    }
    S_ = S_prev_;
}

double CBTensionModelBestel::CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) {
    SaveStateVariablesAsNeeded(time);
    S_.t = time;
    TFloat dt = S_.t - S_prev_.t;
    
    TFloat t = S_.t;
    while (t > cycleLength_) {
        t -=  cycleLength_;
    }
    
    // compute / integrate variables
    TFloat S_p = 0.5 * (1 + tanh((t - onsetSystole_)/steepness_) );
    TFloat S_m = 0.5 * (1 - tanh((t - onsetDiastole_)/steepness_) );
    S_.f       = f(S_p, S_m);
    S_.a       = a(S_.f);
    S_.stress  = RungeKutta4Step(S_prev_.stress, dt, [this](double u) {return dTaudt(u, S_.a, contractility_);});
    
    S_curr_ = S_;
    
    return Tmax_ * S_.stress;
}

void CBTensionModelBestel::InitParamsFromXml(ParameterMap *parameters, std::string parameterKey,
                                             std::string parameterKeyFallback) {
    startTime_   = InitKey(parameters, parameterKey, parameterKeyFallback, ".startTime", 0.0);
    cycleLength_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".cycleLength", 1.0);
    kff_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".kff", 1.0);
    kss_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".kss", 0.0);
    knn_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".knn", 0.0);
    ksn_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".ksn", 0.0);
    
    SetStressCoefficients(Matrix3<TFloat> {kff_, 0, 0,  0, kss_, ksn_,  0, ksn_, knn_});
    AddToActivationTime(startTime_);
    
    // Bestel2001:
    contractility_ = parameters->Get<TFloat>(parameterKey + ".Contractility", 1e5); // in Pa
    maxActivationRate_ = parameters->Get<TFloat>(parameterKey + ".MaxActivationRate", 5);
    minActivationRate_ = parameters->Get<TFloat>(parameterKey + ".MinActivationRate", -30);
    steepness_ = parameters->Get<TFloat>(parameterKey + ".Steepness", 5e-3); // in s
    onsetSystole_ = parameters->Get<TFloat>(parameterKey + ".OnsetSystole", 0.17); // in s
    onsetDiastole_ = parameters->Get<TFloat>(parameterKey + ".OnsetDiastole", 0.484); // in s
    
    // init variables for exporting
    std::vector<TInt> elementIndicesVec = parameters->GetArray<TInt>(parameterKey+".ExportIndices", std::vector<TInt>());
    std::set<TInt> elementIndicesSet(elementIndicesVec.begin(), elementIndicesVec.end());
    auto search = elementIndicesSet.find(ei_);
    if (search != elementIndicesSet.end())
        doExport_ = true;
    if ((parameters->IsAvailable(parameterKey + ".Filename") == false) &&
        (parameters->IsAvailable(parameterKeyFallback + ".Filename") == true) )
        filename_ = parameters->Get<std::string>(parameterKeyFallback);
    else
        filename_ = parameters->Get<std::string>(parameterKey+".Filename", "Bestel.dat");
} // CBTensionModelBestel::InitParamsFromXml

inline void CBTensionModelBestel::WriteToFile(const StateVariables &S) {
    std::ofstream file;
    
    file.open(filename_.c_str(), std::ios::app);
    if (!file.good())
        throw std::runtime_error("CBTensionModelBestel::WriteToFile: Couldn't create " + filename_ + ".");
    
    // write header
    if (!headerWritten_) {
        file << "index";
        file << "\t" << "time";
        file << "\t" << "activation";
        file << "\t" << "indicator";
        file << "\t" << "stress";
        file << "\t" << std::endl;
        headerWritten_ = true;
    }
    
    // write content to file
    file << ei_;
    file << "\t" << S.t;
    file << "\t" << S.a;
    file << "\t" << S.f;
    file << "\t" << S.stress;
    file << std::endl;
    file.close();
} // CBTensionModelBestel::WriteToFile
