/* -------------------------------------------------------

   Land17.cpp

   Ver. 1.1.0

   Created:       Robin Moss (07.06.2017)
   Last modified: Tobias Gerach (08.02.2021)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

// based on the available Matlab code from cemrg.co.uk

#include <Land17.h>


Land17::Land17(Land17Parameters *pp) {
  L17p = pp;
  Init();
}

void Land17::Init() {
  #if KADEBUG
  std::cerr << "#initializing Class: Land17 ... " << endl;
  #endif  // if KADEBUG
  XS        = (v(VT_XS_init));
  XW        = (v(VT_XW_init));
  TRPN      = (v(VT_TRPN_init));
  TmBlocked = (v(VT_TmBlocked_init));
  ZETAS     = (v(VT_ZETAS_init));
  ZETAW     = (v(VT_ZETAW_init));
  Cd        = (v(VT_Cd_init));
}

ML_CalcType Land17::Calc(double tinc, ML_CalcType stretch, ML_CalcType velocity, ML_CalcType &Ca_i, int euler = 1) {
  double HT  = tinc*1000.0;    // sec to msec
  double VEL = velocity/1000.0;  // 1/s to 1/ms


  // XB model
  // -------------------------------------------------------------------------------
  double lambda_m = std::min(1.2, stretch);
  double overlap  = 1.0 + v(VT_beta_0) * (lambda_m + std::min(0.87, lambda_m) - 1.87);
  double h        = std::max(0.0, overlap);

  // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
  double XU    = (1.0 - TmBlocked) - XW - XS;
  double xb_ws = v(VT_k_ws) * XW;
  double xb_uw = v(VT_k_uw) * XU;
  double xb_wu = v(VT_k_wu) * XW;
  double xb_su = v(VT_k_su) * XS;

  double zs_pos      = (ZETAS > 0.0) * ZETAS;
  double zs_neg      = (ZETAS < -1.0) * (-ZETAS - 1.0);
  double zs_         = std::max(zs_pos, zs_neg);  // should be zero if ZETAS + 1 is in [0,1] interval
  double gamma_rate  = (v(VT_gamma)) * zs_;
  double xb_su_gamma = gamma_rate * XS;
  double diff_XS     = xb_ws - xb_su - xb_su_gamma;

  XS += HT * diff_XS;

  double gr_w_        = ((ZETAW < 0.0) ? -ZETAW : ZETAW); // use absolute value of ZetaW
  double gamma_rate_w = (v(VT_gamma_wu)) * gr_w_;  // weak xbs don't like being strained
  double xb_wu_gamma  = gamma_rate_w * XW;
  double diff_XW      = xb_uw - xb_wu - xb_ws - xb_wu_gamma;

  XW += HT * diff_XW;

  double ca50_     = (v(VT_ca50) + v(VT_beta_1) * (lambda_m - 1.0))
    #ifdef HFM
    * v(VT_HFM_Multiplier)  // Change values for HF 0.6
    #endif  // ifdef HF
  ;
  double diff_TRPN = v(VT_koff) * (pow((Ca_i/ca50_), v(VT_TRPN_n)) * (1.0 - TRPN) - TRPN);

  TRPN += HT * diff_TRPN;

  double trpn_np_       = pow(TRPN, -v(VT_nperm)/2.0);
  double trpn_np        = std::min(100.0, trpn_np_);
  double diff_TmBlocked = v(VT_ktm_block) * trpn_np * XU - v(VT_ktm_unblock) * pow(TRPN, v(VT_nperm)/2.0) * TmBlocked;

  TmBlocked += HT * diff_TmBlocked;

  // Velocity dependence -- assumes distortion resets on W->S
  double diff_ZETAS = v(VT_A) * VEL - v(VT_cds) * ZETAS;  // - gamma_rate * ZETAS;
  double diff_ZETAW = v(VT_A) * VEL - v(VT_cdw) * ZETAW;  // - gamma_rate_w * ZETAW;

  ZETAS += HT * diff_ZETAS;
  ZETAW += HT * diff_ZETAW;

  // Active and Total Force
  Ta = h * (v(VT_Tref) / v(VT_dr)) * ((ZETAS + 1.0) * XS + ZETAW * XW);

  // Minimal implementation of the passive cell model
  // Similar to a standard linear solid model. It is used for the viscoelastic response.
  double C_s     = (stretch - 1.0) - Cd;
  double eta     = (C_s > 0.0) ? v(VT_eta_l) : v(VT_eta_s);
  double diff_Cd = v(VT_k) * C_s / eta;

  Cd += HT * diff_Cd;

  double F_d = v(VT_a) * v(VT_k) * C_s;


  // Total Tension
  T = Ta + F_d;

  return T;
}  // Land17::Calc

void Land17::Print(ostream &tempstr) {
  tempstr << XS << ' ' << XW << ' ' << TRPN
          << ' ' << TmBlocked << ' ' << ZETAS << ' ' << ZETAW << ' ' << Ta
          << ' ' << T << ' ' << Cd;
}

void Land17::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] = {"XS", "XW", "TRPN", "TmBlocked", "ZETAS", "ZETAW", "Ta", "T", "Cd"};

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
