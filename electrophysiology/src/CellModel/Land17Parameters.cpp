/* -------------------------------------------------------

   Land17Parameters.cpp

   Ver. 1.1.0

   Created:       Robin Moss (07.06.2017)
   Last modified: Tobias Gerach (08.02.2021)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

// based on the available Matlab code from cemrg.co.uk

#include <Land17Parameters.h>
Land17Parameters::Land17Parameters(const char *initFile) {
  P = new Parameter[vtLast];
  Init(initFile);
}

void Land17Parameters::PrintParameters() {
  cout << "Land17Parameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value <<endl;
  }
}

void Land17Parameters::Init(const char *initFile) {
#if KADEBUG
  cerr << "Loading the Land17 parameter from " << initFile << " ...\n";
#endif  // if KADEBUG
  P[VT_perm50].name      = "perm50";
  P[VT_TRPN_n].name      = "TRPN_n";
  P[VT_koff].name        = "koff";
  P[VT_dr].name          = "dr";
  P[VT_wfrac].name       = "wfrac";
  P[VT_A_eff].name       = "A_eff";
  P[VT_ktm_unblock].name = "ktm_unblock";
  P[VT_beta_1].name      = "beta_1";
  P[VT_beta_0].name      = "beta_0";
  P[VT_gamma].name       = "gamma";
  P[VT_gamma_wu].name    = "gamma_wu";
  P[VT_phi].name         = "phi";
  P[VT_nperm].name       = "nperm";
  P[VT_ca50].name        = "ca50";
  P[VT_Tref].name        = "Tref";
  P[VT_nu].name          = "nu";
  P[VT_mu].name          = "mu";
  P[VT_xi].name          = "xi";
  P[VT_a].name           = "a";
  P[VT_k].name           = "k";
  P[VT_eta_l].name       = "eta_l";
  P[VT_eta_s].name       = "eta_s";

  /// parameters for HFM by Bollen et al.
  #ifdef HFM
  P[VT_HFM_Multiplier].name = "HFM_Multiplier";
  #endif // ifdef HFM

  P[VT_k_ws].name           = "VT_k_ws";
  P[VT_k_uw].name           = "VT_k_uw";
  P[VT_cdw].name            = "VT_cdw";
  P[VT_cds].name            = "VT_cds";
  P[VT_k_wu].name           = "VT_k_wu";
  P[VT_k_su].name           = "VT_k_su";
  P[VT_A].name              = "VT_A";
  P[VT_XS_init].name        = "XS_init";
  P[VT_XW_init].name        = "XW_init";
  P[VT_TRPN_init].name      = "TRPN_init";
  P[VT_TmBlocked_init].name = "TmBlocked_init";
  P[VT_ZETAS_init].name     = "ZETAS_init";
  P[VT_ZETAW_init].name     = "ZETAW_init";
  P[VT_Cd_init].name        = "Cd_init";

  P[VT_k_ws].readFromFile      = false;
  P[VT_k_uw].readFromFile      = false;
  P[VT_cdw].readFromFile       = false;
  P[VT_cds].readFromFile       = false;
  P[VT_k_wu].readFromFile      = false;
  P[VT_k_su].readFromFile      = false;
  P[VT_A].readFromFile         = false;
  P[VT_XSSS].readFromFile      = false;
  P[VT_XWSS].readFromFile      = false;
  P[VT_ktm_block].readFromFile = false;


  ParameterLoader FPL(initFile, FMT_Land17);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = FPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();

#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // Land17Parameters::Init

void Land17Parameters::Calculate() {
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
#if KADEBUG
  cerr << "#Land17Parameters - Calculate ..." << endl;
#endif  // if KADEBUG

  P[VT_k_ws].value = 0.004 * P[VT_mu].value * P[VT_xi].value;
  P[VT_k_uw].value = 0.026  * P[VT_nu].value;

  P[VT_cdw].value = P[VT_phi].value * P[VT_k_uw].value * (1.0-P[VT_dr].value)*(1.0-P[VT_wfrac].value) /
    ((1.0-P[VT_dr].value)*P[VT_wfrac].value);
  P[VT_cds].value = P[VT_phi].value * P[VT_k_ws].value * ((1.0-P[VT_dr].value)*P[VT_wfrac].value) / P[VT_dr].value;

  P[VT_k_wu].value = P[VT_k_uw].value * (1.0/P[VT_wfrac].value - 1.0) - P[VT_k_ws].value;
  P[VT_k_su].value = P[VT_k_ws].value * (1.0/P[VT_dr].value - 1.0) * P[VT_wfrac].value;

  P[VT_A].value = (P[VT_A_eff].value) *
    (P[VT_dr].value) / ((1.0-P[VT_dr].value) * P[VT_wfrac].value + P[VT_dr].value);

  P[VT_XSSS].value      = P[VT_dr].value * 0.5;
  P[VT_XWSS].value      = (1.0 - P[VT_dr].value) * P[VT_wfrac].value * 0.5;
  P[VT_ktm_block].value = P[VT_ktm_unblock].value * pow(P[VT_perm50].value, P[VT_nperm].value) * 0.5 /
    (0.5 - P[VT_XSSS].value - P[VT_XWSS].value);
}  // Land17Parameters::Calculate
