/*
 *  CBCircModel.hpp
 *  CardioMechanics
 *
 *  Created by Steffen Schuler on 07.01.16.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include <algorithm>
#include "filesystem.h"

#include "DCCtrl.h"

template<TInt nC, TInt nS, TInt nR>
CBCircModel<nC, nS, nR>::CBCircModel(ParameterMap *parameters, std::string parameterPrefix) {
  filename_    = parameters->Get<std::string>(parameterPrefix + ".ExportFile", "");
  maxTimeStep_ = parameters->Get<TFloat>(parameterPrefix + ".MaxIntegrationTimeStep", 1e-4);
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::InitSteadyStateCheck(ParameterMap *parameters, std::string parameterPrefix, TInt iSVD) {
  steadyChkActive_ = parameters->Get<bool>(parameterPrefix + ".SteadyStateCheck.Active", false);
  if (steadyChkActive_) {
    steadyChkStart_  = parameters->Get<TFloat>(parameterPrefix + ".SteadyStateCheck.StartTime");
    steadyChkPeriod_ = parameters->Get<TFloat>(parameterPrefix + ".SteadyStateCheck.Period");
    steadyChkThresh_ = parameters->Get<TFloat>(parameterPrefix + ".SteadyStateCheck.Threshold");

    iSVD_ = iSVD;
    std::string steadyChkMode;
    if ((iSVD_ < 0) || (iSVD_ >= nS_) ) {
      steadyChkMode = parameters->Get<std::string>(parameterPrefix + ".SteadyStateCheck.Mode", "TemporalDifference");
      if (steadyChkMode != "TemporalDifference") {
        throw std::runtime_error(
                "CBCircModel::CBCircModel(): Only TemporalDifference can be used as steady state check mode with " + name_ +
                ".");
      }
    } else {
      steadyChkMode =
        parameters->Get<std::string>(parameterPrefix + ".SteadyStateCheck.Mode", "StrokeVolumeDifference");
    }

    if (steadyChkMode == "TemporalDifference") {
      steadyChkTemporal_ = true;
      std::string paramName = parameters->Get<std::string>(parameterPrefix + ".SteadyStateCheck.Parameter");
      std::string *paramPos = std::find(stateVarsNames_, stateVarsNames_+nS_, paramName);
      if (paramPos != stateVarsNames_+nS_) {
        steadyChkParam_ = stateVars_ + (paramPos-stateVarsNames_);
      } else {
        paramPos = std::find(resultsNames_, resultsNames_+nR_, paramName);
        if (paramPos != resultsNames_+nR_) {
          steadyChkParam_ = *(results_ + (paramPos-resultsNames_));
        } else {
          throw std::runtime_error(
                  "CBCircWholeHeart::CBCircWholeHeart(): " + paramName +
                  " is not a valid parameter and thus cannot be used for steady state check.");
        }
      }
    }
  }
} // >::InitSteadyStateCheck

template<TInt nC, TInt nS, TInt nR>
bool CBCircModel<nC, nS, nR>::SteadyStateCheck(TFloat time) {
  bool inSteadyState = false;

  if (steadyChkActive_) {
    if (std::abs(time-steadyChkStart_) < 0.5*(time-prevSteadyChkTime_)) {
      if (steadyChkTemporal_)
        prevSteadyChkParam_ = *steadyChkParam_;
      else
        estStateVars_[iSVD_] = 0;

      prevSteadyChkTimeActive_ = steadyChkStart_;
    } else if ((time > steadyChkStart_) &&
               (std::abs(time-prevSteadyChkTimeActive_-steadyChkPeriod_) < 0.5*(time-prevSteadyChkTime_))) {
      TFloat steadyChkVal;

      if (steadyChkTemporal_) {
        // Check using the difference of the chosen parameter between two succeeding cycles
        steadyChkVal = std::abs(*steadyChkParam_-prevSteadyChkParam_);
        prevSteadyChkParam_ = *steadyChkParam_;
      } else {
        // Check using stroke volume difference (between LV and RV)
        steadyChkVal = std::abs(estStateVars_[iSVD_]);
        estStateVars_[iSVD_] = 0;
      }

      if (steadyChkVal < steadyChkThresh_) {
        inSteadyState = true;
        DCCtrl::print << name_ << " is in steady state!" << std::endl;
      } else {
        DCCtrl::print << name_ << " is not yet in steady state! Difference is " << steadyChkVal << " ml." << std::endl;
      }

      prevSteadyChkTimeActive_ = time;
    }
  }
  prevSteadyChkTime_ = time;

  return inSteadyState;
} // >::SteadyStateCheck

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::GetCavitySurfaceIndices(TInt cavitySurfaceIndices[]) {
  for (TInt i = 0; i < nC_; i++)
    cavitySurfaceIndices[i] = cavitySurfaceIndices_[i];
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::GetPreloadingPressures(TFloat preloadingPressures[]) {
  for (TInt i = 0; i < 2*nC_; i++)
    preloadingPressures[i] = cavityPreloadingPressures_[i];
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::AcceptEstState() {
  for (TInt i = 0; i < nS_; i++) {
    prevStateVars_[i] = stateVars_[i];
    stateVars_[i] = estStateVars_[i];
  }

  for (TInt i = 0; i < nC_; i++) {
    prevCavityPressures_[i] = cavityPressures_[i];
    cavityPressures_[i] = *(estCavityPressures_[i]);
  }
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::RevertPrevState() {
  // Not really needed

  for (TInt i = 0; i < nS_; i++)
    stateVars_[i] = prevStateVars_[i];

  for (TInt i = 0; i < nC_; i++)
    cavityPressures_[i] = prevCavityPressures_[i];
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::WriteHeaderToFile() {
  if (filename_ != "") {
    file_.open(filename_.c_str());
    if (!file_.good())
      throw std::runtime_error("CBCircModel::WriteHeaderToFile: Couldn't create " + filename_ + ".");

    file_ << std::setw(18) << "time";

    for (TInt i = 0; i < nS_; i++)
      file_ << std::setw(18) << stateVarsNames_[i];
    for (TInt i = 0; i < nR_; i++)
      file_ << std::setw(18) << resultsNames_[i];

    file_.close();
  } else {
    throw std::runtime_error("CBCircModel::WriteHeaderToFile: Filename is empty.");
  }
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::WriteResultsToFile(TFloat time) {
  if (filename_ != "") {
    file_.open(filename_.c_str(), std::ios::app);
    if (!file_.good())
      throw std::runtime_error("CBCircModel::WriteResultsToFile: Couldn't create " + filename_ + ".");

    file_ << std::endl << std::setprecision(5) << std::fixed << std::setw(18) << time;

    file_ << std::setprecision(9) << std::fixed;
    for (TInt i = 0; i < nS_; i++)
      file_ << std::setw(18) << estStateVars_[i];
    for (TInt i = 0; i < nR_; i++)
      file_ << std::setw(18) << *(results_[i]);

    file_.close();
  } else {
    throw std::runtime_error("CBCircModel::WriteResultsToFile: Filename is empty.");
  }
}

template<TInt nC, TInt nS, TInt nR>
void CBCircModel<nC, nS, nR>::Integrate(TFloat cavityPressures[], TFloat timeStep) {
  TInt nSteps = ceil(timeStep / maxTimeStep_);
  TFloat dt = timeStep / nSteps;

  // TFloat dp[nC_];
  for (TInt i = 0; i < nC_; i++) {
    *(estCavityPressures_[i]) = cavityPressures[i];

    // dp[i] = (cavityPressures[i]-cavityPressures_[i]) / nSteps;
  }

  for (TInt i = 0; i < nS_; i++)
    estStateVars_[i] = stateVars_[i];

  for (TInt s = 0; s < nSteps; s++) {
    TFloat udot0[nS_], udot1[nS_], udot2[nS_], udot3[nS_];
    TFloat u1[nS_], u2[nS_], u3[nS_];
    TFloat p1[nC_], p2[nC_], p3[nC_];

    for (TInt i = 0; i < nC_; i++) {
      //// interpolate pressures linearly at points needed for RK4
      // p1[i] = cavityPressures_[i] + s * dp[i];
      // p2[i] = p1[i] + dp[i] / 2.0;
      // p3[i] = p1[i] + dp[i];

      p1[i] = cavityPressures[i];
      p2[i] = cavityPressures[i];
      p3[i] = cavityPressures[i];
    }

    Rates(p1, estStateVars_, udot0);

    for (TInt i = 0; i < nS_; i++)
      u1[i] = estStateVars_[i] + dt * udot0[i] / 2.0;
    Rates(p2, u1, udot1);

    for (TInt i = 0; i < nS_; i++)
      u2[i] = estStateVars_[i] + dt * udot1[i] / 2.0;
    Rates(p2, u2, udot2);

    for (TInt i = 0; i < nS_; i++)
      u3[i] = estStateVars_[i] + dt * udot2[i];
    Rates(p3, u3, udot3);

    for (TInt i = 0; i < nS_; i++)
      estStateVars_[i] += dt * (udot0[i] + 2.0*udot1[i] + 2.0*udot2[i] + udot3[i]) / 6.0;
  }

  // calculate results_
  Algebraics(cavityPressures, estStateVars_);
} // >::Integrate
