/*
 *  CBTensionModelLumens.cpp
 *  CardioMechanics
 *
 * Created by Lukas Baron on Wed May 31 2017.
 *
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBTensionModelLumens.h"
#include <cassert>

TFloat CBLumensTension::Ls(TFloat eps_f) {
    TFloat x = Lsref_ * exp(eps_f);
    // TFloat x = Lsref_ * (1 + eps_f);
    return x;
}

TFloat CBLumensTension::Frise(TFloat t) {
    while (t > cycleLength_)
        t = t - cycleLength_;
    TFloat x = fmin(8., fmax(0., t / tauR_));
    return 0.02 * pow(x, 3) * pow(8 - x, 2) * exp(-x);
}

TFloat CBLumensTension::dLscdt(TFloat Lsc, TFloat eps_f) {
    TFloat x = ((Ls(eps_f) - Lsc) / Lseiso_ - 1.) * vmax_;
    // if (doExport_)
    //std::cout << " eps_f: " << eps_f << " dLscdt: " << x << std::endl;
    return x;
}

TFloat CBLumensTension::CL(TFloat Lsc) {
    return tanh(4. * 1e12 * pow(Lsc - Lsc0_, 2.));	// 1e12 is due to SI units
}

TFloat CBLumensTension::T_fun(TFloat Lsc) {
    return tausc_ * (0.29 + 0.3 * Lsc * 1e6); // 1e6: SI units
}

TFloat CBLumensTension::dCdt(TFloat C, TFloat Lsc, TFloat eps_f, TFloat t) {
    while (t > cycleLength_)
        t = t - cycleLength_;
    
    return 1 / tauR_ * CL(Lsc) * Frise(t)
    + 1 / tauD_ * (Crest_ - C) / (1 + exp((T_fun(Lsc) - t) / tauD_));
}

TFloat CBLumensTension::sigma_fact(TFloat Lsc, TFloat C, TFloat eps_f) {
    if (Lsc < Lsc0_)
        return 0;
    
    return sigma_act_ * C * (Lsc - Lsc0_) * 1e6 * (Ls(eps_f) - Lsc) / Lseiso_; // 1e6: SI units
}

/// Read a value from xml with fallback. Priorities: key+tail > keyFallback+tail > defaultValue
static TFloat InitKey(ParameterMap* params, std::string key, std::string keyFallback, std::string tail, TFloat defaultValue) {
    if(params->IsAvailable(key + tail) == false && params->IsAvailable(keyFallback + tail) == true)
        return params->Get<TFloat>(keyFallback);
    else
        return params->Get<TFloat>(key+tail, defaultValue);
}

CBLumensTension::CBLumensTension(CBElement* e, ParameterMap* params) {
    
    // set model parameters
    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    int mi = e->GetMaterialIndex();
    ei_ = e->GetIndex();
    std::string key = "Materials.Mat_" + std::to_string(mi);
    std::string keyFallback = "Materials.Mat_Default";
    InitParamsFromXml(params, key+".Lumens", keyFallback+".Lumens");
    
    // set initial values for the state variables
    S_.eps_f = 0.;
    S_.t = 0.;
    S_.f_rise = Frise(S_.t);
    S_.ls = Ls(S_.eps_f);
    S_.lsc = Ls(S_.eps_f) - Lseiso_; // chosen s.t. dLscdt(S.lsc) = 0
    S_.cl = CL(S_.lsc);
    S_.T = T_fun(S_.lsc);
    S_.c = Crest_; // chosen s.t. dCdt(S.c) = 0
    S_.sigmafact = sigma_fact(S_.lsc, S_.c, S_.eps_f);
    
    S_curr_ = S_;
    S_prev_ = S_;
}

void CBLumensTension::SaveStateVariablesAsNeeded(TFloat time) {
    // time is bigger than in the last call? -> the last call was successful, store S_curr to S_prev
    // time is smaller than in last call? -> previous run was non-successful, discard last computed values S_curr
    // time is roughly the same? -> update the current guess S_curr (with what?) ...
    if (time - S_curr_.t > 0)
    {
        S_prev_ = S_curr_;
        // export state variables on a successful export
        if (doExport_)
            this->WriteToFile(S_prev_);
    }
    else if (time - S_prev_.t < 0) {
        S_curr_ = S_prev_;
    } else {
        S_curr_ = S_prev_;
    }
    S_ = S_prev_;
    
}

double CBLumensTension::CalcActiveTension(const math_pack::Matrix3<double>& deformation, const double time) {
    SaveStateVariablesAsNeeded(time);
    S_.t = time;
    TFloat delta_t = S_.t - S_prev_.t;
    
    // compute / integrate variables
    // The order of the formula was rearranged such that the new computed data is used whenever available.
    S_.eps_f = sqrt(deformation.GetCol(0)*deformation.GetCol(0))-1;
    S_.f_rise = Frise(S_.t - GetActivationTime());
    S_.ls = Ls(S_.eps_f);
    S_.lsc = RungeKutta4Step(S_prev_.lsc, delta_t, [this](double u){return dLscdt(u, S_.eps_f);} );  // ode, solved with expl rk4
    if (S_.t < GetActivationTime())
        S_.lsc = Ls(S_.eps_f) - Lseiso_; // chosen s.t. dLscdt(S.lsc) = 0
    S_.cl = CL(S_.lsc);
    S_.T = T_fun(S_.lsc);
    S_.c = RungeKutta4Step(S_prev_.c, delta_t, [this](double u){return dCdt(u, S_.lsc, S_.eps_f, S_.t - GetActivationTime());} );  // ode, solved with expl rk4
    if (S_.t < GetActivationTime())
        S_.c = Crest_; // chosen s.t. dCdt(S.c) = 0
    S_.sigmafact = sigma_fact(S_.lsc, S_.c, S_.eps_f);
    S_curr_ = S_;
    
    // if (doExport_)
    //   std::cout << Tmax_*S_.sigmafact << " ";
    return Tmax_*S_.sigmafact;
}

void CBLumensTension::InitParamsFromXml(ParameterMap* parameters, std::string parameterKey, std::string parameterKeyFallback) {
    cycleLength_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".CycleLength", 0.8);
    TFloat startTime = InitKey(parameters, parameterKey, parameterKeyFallback, ".StartTime", 0.0);
    AddToActivationTime(startTime);
    
    // Lumens2009: Crest = 0.02, but we want a zero tension at rest
    Crest_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".Crest", 0.0);
    Lsc0_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".Lsc0" , 1.51 * 1e-6);
    Lsref_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".Lsref", 2.0 * 1e-6);
    Lseiso_      = InitKey(parameters, parameterKey, parameterKeyFallback, ".Lseiso", 0.04 * 1e-6);
    vmax_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".vmax", 7. * 1e-6);
    tauD_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".tauD", 32.0 * 1e-3);
    tauR_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".tauR", 48.0 * 1e-3);
    tausc_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".tausc", 425.0 * 1e-3);
    
    //sigma_act_   = InitKey(parameters, parameterKey, parameterKeyFallback, ".sigma_act", 120.0 * 1e3);
    sigma_act_   = 1; // we want lumens just to compute a 'normalized tension', scaling is done with Tmax
    sigma_pas_   = InitKey(parameters, parameterKey, parameterKeyFallback, ".sigma_pas", 7.0 * 1e3);
    
    if (!(cycleLength_ >= (8 * tauR_))) {  // needed for Frise() to go back to zero
        throw std::runtime_error("CBLumensTension::Frise() : CycleLength was < 8 * tauR_.");
    }
    
    // init variables for exporting
    std::vector<TInt> elementIndicesVec = parameters->GetArray<TInt>(parameterKey+".ExportIndices", std::vector<TInt>());
    std::set<TInt> elementIndicesSet(elementIndicesVec.begin(), elementIndicesVec.end());
    auto search = elementIndicesSet.find(ei_);
    if(search != elementIndicesSet.end())
        doExport_ = true;
    if(parameters->IsAvailable(parameterKey + ".Filename") == false && parameters->IsAvailable(parameterKeyFallback + ".Filename") == true)
        filename_ = parameters->Get<std::string>(parameterKeyFallback);
    else
        filename_ = parameters->Get<std::string>(parameterKey+".Filename", "lumens.dat");
    
}


inline void CBLumensTension::WriteToFile(const StateVariables& S) {
    std::ofstream file;
    file.open(filename_.c_str(), std::ios::app);
    if (!file.good())
        throw std::runtime_error("CBTensionModelLumens::WriteToFile: Couldn't create " + filename_ + ".");
    
    // write header
    if (!headerWritten_) {
        file << "\t" << "index";
        file << "\t" << "t";
        file << "\t" << "eps_f";
        file << "\t" << "f_rise";
        file << "\t" << "ls";
        file << "\t" << "lsc";
        file << "\t" << "cl";
        file << "\t" << "T";
        file << "\t" << "c";
        file << "\t" << "sigmafact";
        file << "\t" << std::endl;
        headerWritten_ = true;
    }
    
    // write content to file
    file << "\t" << ei_;
    file << "\t" << S.t;
    file << "\t" << S.eps_f;
    file << "\t" << S.f_rise;
    file << "\t" << S.ls;
    file << "\t" << S.lsc;
    file << "\t" << S.cl;
    file << "\t" << S.T;
    file << "\t" << S.c;
    file << "\t" << S.sigmafact;
    file << std::endl;
    file.close();
    
    // show content on terminal, too
    std::cout << "  " << ei_;
    std::cout << "  " << S.t;
    // std::cout << "  " << S.f_rise;
    std::cout << "  " << S.eps_f;
    // std::cout << "\t" << S.ls;
    // std::cout << "\t" << S.lsc;
    // std::cout << "\t" << S.cl;
    // std::cout << "\t" << S.T;
    // std::cout << "\t" << S.c;
    std::cout << "  " << S.sigmafact;
    std::cout << std::endl;
}
