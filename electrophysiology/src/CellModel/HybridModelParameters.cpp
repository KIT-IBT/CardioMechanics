/*
 * File: HybridModelParameters.cpp
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


#include <HybridModelParameters.h>

HybridModelParameters::HybridModelParameters(const char *initFile) {
  P = new Parameter[vtLast];
  Init(initFile);
}

void HybridModelParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"HybridModelParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}

void HybridModelParameters::Init(const char *initFile) {
  P[VT_AMATPa].name      = "AMATP";
  P[VT_MATPa].name       = "MATP";
  P[VT_MADPPa].name      = "MADPP";
  P[VT_AwMADPPa].name    = "AwMADPP";
  P[VT_AsMADPPa].name    = "AsMADPP";
  P[VT_AsMADPa].name     = "AsMADP";
  P[VT_AMADPa].name      = "AMADP";
  P[VT_TCaa].name        = "TCa";
  P[VT_TMona].name       = "TMon";
  P[VT_MADPa].name       = "MADP";
  P[VT_Ma].name          = "M";
  P[VT_k_1].name         = "k_1";
  P[VT_k_m1].name        = "k_m1";
  P[VT_k_2].name         = "k_2";
  P[VT_k_3].name         = "k_3";
  P[VT_k_m3].name        = "k_m3";
  P[VT_k_4].name         = "k_4";
  P[VT_k_m4].name        = "k_m4";
  P[VT_k_5].name         = "k_5";
  P[VT_k_m5].name        = "k_m5";
  P[VT_k_6].name         = "k_6";
  P[VT_k_m6].name        = "k_m6";
  P[VT_k_7].name         = "k_7";
  P[VT_k_8].name         = "k_8";
  P[VT_k_m8].name        = "k_m8";
  P[VT_k_on].name        = "k_on";
  P[VT_k_off].name       = "k_off";
  P[VT_k_9].name         = "k_9";
  P[VT_k_10].name        = "k_10";
  P[VT_k_11].name        = "k_11";
  P[VT_k_12].name        = "k_12";
  P[VT_k_13].name        = "k_13";
  P[VT_k_14].name        = "k_14";
  P[VT_tm_on].name       = "tm_on";
  P[VT_tm_off].name      = "tm_off";
  P[VT_TCaMax].name      = "TCaMax";
  P[VT_TCaMin].name      = "TCaMin";
  P[VT_Fmax].name        = "Fmax";
  P[VT_dFmax].name       = "dFmax";
  P[VT_N_v].name         = "N_v";
  P[VT_v50].name         = "v50";
  P[VT_TCa_stretch].name = "TCa_stretch";
  P[VT_TMon_coop].name   = "TMon_coop";
  P[VT_TMon_pow].name    = "TMon_pow";
  P[VT_ATP].name         = "ATP";
  P[VT_F_physiol].name   = "F_physiol";
  P[VT_k7_base].name     = "k7_base";
  P[VT_k7_stretch].name  = "k7_stretch";
  P[VT_k7_force].name    = "k7_force";
  P[VT_k5_stretch].name  = "k5_stretch";
  P[VT_k5_xb].name       = "k5_xb";
  P[VT_detach_vel].name  = "detach_vel";

  P[VT_dFmax].readFromFile = false;

  ParameterLoader FPL(initFile, FMT_Hybrid);
  this->setOverlapParameters(FPL.getOverlapString());
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = FPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  Calculate();
}  // HybridModelParameters::Init

void HybridModelParameters::Calculate() {
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();

  P[VT_dFmax].value = 1.0/P[VT_Fmax].value;
}
