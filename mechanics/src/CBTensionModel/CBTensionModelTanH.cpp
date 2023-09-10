/*
 *  CBTensionModelTanH.cpp
 *  CardioMechanics
 *
 *  Created by Tobias Gerach, Fr. 08. May 2021
 *
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBTensionModelTanH.h"
#include <algorithm>
#include <cassert>

/// Read a value from xml with fallback. Priorities: key+tail > keyFallback+tail > defaultValue
static TFloat InitKey(ParameterMap *params, std::string key, std::string keyFallback, std::string tail,
                      TFloat defaultValue) {
    if ((params->IsAvailable(key + tail) == false) && (params->IsAvailable(keyFallback + tail) == true))
        return params->Get<TFloat>(keyFallback);
    else
        return params->Get<TFloat>(key+tail, defaultValue);
}

CBTensionModelTanH::CBTensionModelTanH(CBElementSolid *e, ParameterMap *params) {
    // set model parameters
    Matrix3<TFloat> deformationTensor;
    
    e->GetDeformationTensor(deformationTensor);
    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    if (Tmax_ != 1000) {
        throw std::runtime_error(
                                 "CBTensionModelTanH: Parameter TensionMax != 1000 is not compatible with this model to ensure model output in Pa!");
    }
    int mi = e->GetMaterialIndex();
    ei_ = e->GetIndex();
    std::string key         = "Materials.Mat_" + std::to_string(mi);
    std::string keyFallback = "Materials.Mat_Default";
    InitParamsFromXml(params, key+".TanH", keyFallback+".TanH");
    
    // set initial values for the state variables
    S_.t      = 0.0;
    S_.lambda = sqrt(deformationTensor.GetCol(0)*deformationTensor.GetCol(0));
    S_.Ta     = 0.0;
    
    S_curr_ = S_;
    S_prev_ = S_;
}

void CBTensionModelTanH::SaveStateVariablesAsNeeded(TFloat time) {
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

double CBTensionModelTanH::CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) {
    SaveStateVariablesAsNeeded(time);
    S_.t = time;  // in s
    
    /// update extention ratio (lambda)
    S_.lambda = sqrt(deformation.GetCol(0)*deformation.GetCol(0));  // in m/m
                                                                    /// non-linear length dependence
    TFloat phi = tanh(ld_ * (S_.lambda - l_0_) );
    
    /// contraction time constant
    TFloat tau_c = tau_c0_ + ld_up_ * (1 - phi);
    
    /// onset of contraction
    TFloat t = S_.t;
    while (t > cycleLength_) {
        t -=  cycleLength_;
    }
    TFloat t_s = std::max(std::min(t - GetActivationTime() - t_emd_, t_dur_), 0.0);
    
    /// ActiveTension in kPa
    /// Tmax_ = 1000 should be used to do kPa -> Pa
    /// don't allow negative values
    S_.Ta = std::max(S_peak_ * phi * pow(tanh(t_s / tau_c), 2.0) * pow(tanh((t_dur_ - t_s)/tau_r_), 2.0), 0.0);
    
    S_curr_ = S_;
    return Tmax_ * S_.Ta;
}  // CBTensionModelTanH::CalcActiveTension

void CBTensionModelTanH::InitParamsFromXml(ParameterMap *parameters, std::string parameterKey,
                                           std::string parameterKeyFallback) {
    startTime_   = InitKey(parameters, parameterKey, parameterKeyFallback, ".startTime", 0.0);
    cycleLength_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".cycleLength", 1.0);
    kff_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".kff", 1.0);
    kss_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".kss", 0.0);
    knn_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".knn", 0.0);
    ksn_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".ksn", 0.0);
    
    SetStressCoefficients(Matrix3<TFloat> {kff_, 0, 0,  0, kss_, ksn_,  0, ksn_, knn_});
    AddToActivationTime(startTime_);
    
    l_0_    = InitKey(parameters, parameterKey, parameterKeyFallback, ".MinFiberStretch", 0.7);
    t_emd_  = InitKey(parameters, parameterKey, parameterKeyFallback, ".EMD", 0.015);
    S_peak_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".PeakTension", 120.0);
    t_dur_  = InitKey(parameters, parameterKey, parameterKeyFallback, ".Duration", 0.3);
    tau_c0_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".ContractionTime", 0.1);
    ld_     = InitKey(parameters, parameterKey, parameterKeyFallback, ".DegreeOfLengthDependence", 5.0);
    ld_up_  = InitKey(parameters, parameterKey, parameterKeyFallback, ".LengthDependenceUpstrokeTime", 0.5);
    tau_r_  = InitKey(parameters, parameterKey, parameterKeyFallback, ".RelaxationTime", 0.1);
    
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
        filename_ = parameters->Get<std::string>(parameterKey+".Filename", "TensionModel.dat");
}  // CBTensionModelTanH::InitParamsFromXml

inline void CBTensionModelTanH::WriteToFile(const StateVariables &S) {
    std::ofstream file;
    
    file.open(filename_.c_str(), std::ios::app);
    if (!file.good())
        throw std::runtime_error("CBTensionModelTanH::WriteToFile: Couldn't create " + filename_ + ".");
    
    // write header
    if (!headerWritten_) {
        file << "index";
        file << "\t" << "t";
        file << "\t" << "lambda";
        file << "\t" << "Ta";
        file << "\t" << std::endl;
        headerWritten_ = true;
    }
    
    // write content to file
    file << std::setprecision(14);  // needed to increase precision of export to compare with matlab
    file << ei_;
    file << "\t" << S.t;
    file << "\t" << S.lambda;
    file << "\t" << S.Ta;
    file << std::endl;
    file.close();
}  // CBTensionModelTanH::WriteToFile
