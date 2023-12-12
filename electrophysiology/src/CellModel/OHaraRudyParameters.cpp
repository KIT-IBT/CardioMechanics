/*
 * File: OHaraRudyParameters.cpp
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


#include <OHaraRudyParameters.h>
OHaraRudyParameters::OHaraRudyParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

OHaraRudyParameters::~OHaraRudyParameters() {}

void OHaraRudyParameters::PrintParameters() {
  cout << "OHaraRudyParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void OHaraRudyParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the OHaraRudy parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  P[VT_amplitude].name        = "amplitude";
  P[VT_duration].name         = "duration";
  P[VT_celltype].name         = "celltype";
  P[VT_Na_o].name             = "Na_o";
  P[VT_Ca_o].name             = "Ca_o";
  P[VT_K_o].name              = "K_o";
  P[VT_INaF_Multiplier].name  = "INaF_Multiplier";
  P[VT_INaL_Multiplier].name  = "INaL_Multiplier";
  P[VT_Ito_Multiplier].name   = "Ito_Multiplier";
  P[VT_PCa_Multiplier].name   = "PCa_Multiplier";
  P[VT_IKr_Multiplier].name   = "IKr_Multiplier";
  P[VT_IKs_Multiplier].name   = "IKs_Multiplier";
  P[VT_IK1_Multiplier].name   = "IK1_Multiplier";
  P[VT_INaCa_Multiplier].name = "INaCa_Multiplier";
  P[VT_INaK_Multiplier].name  = "INaK_Multiplier";
  P[VT_IKb_Multiplier].name   = "IKb_Multiplier";
  P[VT_INab_Multiplier].name  = "INab_Multiplier";
  P[VT_ICab_Multiplier].name  = "ICab_Multiplier";
  P[VT_IpCa_Multiplier].name  = "IpCa_Multiplier";

  /// state variables
  P[VT_Init_Vm].name             = "Init_Vm";
  P[VT_Init_Na_i].name           = "Init_Na_i";
  P[VT_Init_Na_ss].name          = "Init_Na_ss";
  P[VT_Init_K_i].name            = "Init_K_i";
  P[VT_Init_K_ss].name           = "Init_K_ss";
  P[VT_Init_Ca_i].name           = "Init_Ca_i";
  P[VT_Init_Ca_ss].name          = "Init_Ca_ss";
  P[VT_Init_Ca_nsr].name         = "Init_Ca_nsr";
  P[VT_Init_Ca_jsr].name         = "Init_Ca_jsr";
  P[VT_Init_m].name              = "Init_m";
  P[VT_Init_h_fast].name         = "Init_h_fast";
  P[VT_Init_h_slow].name         = "Init_h_slow";
  P[VT_Init_j].name              = "Init_j";
  P[VT_Init_h_CaMK_slow].name    = "Init_h_CaMK_slow";
  P[VT_Init_j_CaMK].name         = "Init_j_CaMK";
  P[VT_Init_m_L].name            = "Init_m_L";
  P[VT_Init_h_L].name            = "Init_h_L";
  P[VT_Init_h_L_CaMK].name       = "Init_h_L_CaMK";
  P[VT_Init_a].name              = "Init_a";
  P[VT_Init_i_fast].name         = "Init_i_fast";
  P[VT_Init_i_slow].name         = "Init_i_slow";
  P[VT_Init_a_CaMK].name         = "Init_a_CaMK";
  P[VT_Init_i_CaMK_fast].name    = "Init_i_CaMK_fast";
  P[VT_Init_i_CaMK_slow].name    = "Init_i_CaMK_slow";
  P[VT_Init_d].name              = "Init_d";
  P[VT_Init_f_fast].name         = "Init_f_fast";
  P[VT_Init_f_slow].name         = "Init_f_slow";
  P[VT_Init_f_Ca_fast].name      = "Init_f_Ca_fast";
  P[VT_Init_f_Ca_slow].name      = "Init_f_Ca_slow";
  P[VT_Init_j_Ca].name           = "Init_j_Ca";
  P[VT_Init_n].name              = "Init_n";
  P[VT_Init_f_CaMK_fast].name    = "Init_f_CaMK_fast";
  P[VT_Init_f_Ca_CaMK_fast].name = "Init_f_Ca_CaMK_fast";
  P[VT_Init_x_r_fast].name       = "Init_x_r_fast";
  P[VT_Init_x_r_slow].name       = "Init_x_r_slow";
  P[VT_Init_x_s1].name           = "Init_x_s1";
  P[VT_Init_x_s2].name           = "Init_x_s2";
  P[VT_Init_x_K1].name           = "Init_x_K1";
  P[VT_Init_J_rel_NP].name       = "Init_J_rel_NP";
  P[VT_Init_J_rel_CaMK].name     = "Init_J_rel_CaMK";
  P[VT_Init_CaMK_trap].name      = "Init_CaMK_trap";

  #ifdef ISAC
  P[VT_ISAC_SWITCH].name   =     "ISAC_SWITCH";
  P[VT_alpha].name         =     "alpha";
  P[VT_G_sac].name             =     "G_sac";
  P[VT_K_sac].name             =     "K_sac";
  P[VT_E_sac].name         =     "E_sac";
  #endif // ifdef ISAC

  /// parameters for HF by Gomez et al.
  #ifdef HF
  P[VT_HF_Gomez_1].name = "HF_Gomez_1";
  P[VT_HF_Gomez_2].name = "HF_Gomez_2";
  P[VT_HF_Gomez_3].name = "HF_Gomez_3";
  P[VT_HF_Gomez_4].name = "HF_Gomez_4";
  P[VT_HF_Gomez_5].name = "HF_Gomez_5";
  P[VT_HF_Gomez_6].name = "HF_Gomez_6";
  P[VT_HF_Gomez_7].name = "HF_Gomez_7";
  P[VT_HF_Gomez_8].name = "HF_Gomez_8";
  P[VT_HF_Gomez_9].name = "HF_Gomez_9";
  P[VT_HF_Gomez_0].name = "HF_Gomez_0";
  #endif // ifdef HF

  /// parameters for HFM by Bollen et al.
  #ifdef HFM
  P[VT_HFM_Multiplier].name = "HFM_Multiplier";
  #endif // ifdef HFM

  /// parameters for strong coupling with Land17
  #ifdef TRPN
  P[VT_Init_Ca_TRPN].name = "Init_Ca_TRPN";
  P[VT_beta1].name        = "beta1";
  P[VT_Ca_T50].name       = "Ca_T50";
  P[VT_k_TRPN].name       = "k_TRPN";
  P[VT_n_TRPN].name       = "n_TRPN";
  #endif  // if TRPN

  /// Parameters that are not read from .ev file
  P[VT_R].name              = "R";
  P[VT_T].name              = "T";
  P[VT_F].name              = "F";
  P[VT_L].name              = "L";
  P[VT_r].name              = "r";
  P[VT_tau_h_L].name        = "tau_h_L";
  P[VT_A_h_fast].name       = "A_h_fast";
  P[VT_A_h_slow].name       = "A_h_slow";
  P[VT_v_cell].name         = "v_cell";
  P[VT_A_geo].name          = "A_geo";
  P[VT_A_cap].name          = "A_cap";
  P[VT_v_myo].name          = "v_myo";
  P[VT_v_nsr].name          = "v_nsr";
  P[VT_v_jsr].name          = "v_jsr";
  P[VT_v_ss].name           = "v_ss";
  P[VT_P_CaNa].name         = "P_CaNa";
  P[VT_P_CaK].name          = "P_CaK";
  P[VT_P_Ca_CaMK].name      = "P_Ca_CaMK";
  P[VT_P_CaNa_CaMK].name    = "P_CaNa_CaMK";
  P[VT_P_CaK_CaMK].name     = "P_CaK_CaMK";
  P[VT_A_f_fast].name       = "A_f_fast";
  P[VT_A_f_slow].name       = "A_f_slow";
  P[VT_k_Na1].name          = "k_Na1";
  P[VT_k_Na2].name          = "k_Na2";
  P[VT_k_asymm].name        = "k_asymm";
  P[VT_k_Ca_on].name        = "k_Ca_on";
  P[VT_k_Ca_off].name       = "k_Ca_off";
  P[VT_k_m_1].name          = "k_m_1";
  P[VT_k_p_2].name          = "k_p_2";
  P[VT_k_p_4].name          = "k_p_4";
  P[VT_MgADP].name          = "MgADP";
  P[VT_MgATP].name          = "MgATP";
  P[VT_K_MgATP].name        = "K_MgATP";
  P[VT_beta_tau].name       = "beta_tau";
  P[VT_P_Ca].name           = "P_Ca";
  P[VT_h_10].name           = "h_10";
  P[VT_h_11].name           = "h_11";
  P[VT_h_12].name           = "h_12";
  P[VT_k_1].name            = "k_1";
  P[VT_k_2].name            = "k_2";
  P[VT_k_5].name            = "k_5";
  P[VT_beta_1].name         = "beta_1";
  P[VT_alpha_2].name        = "alpha_2";
  P[VT_alpha_4].name        = "alpha_4";
  P[VT_alpha_rel].name      = "alpha_rel";
  P[VT_beta_tau_CaMK].name  = "beta_tau_CaMK";
  P[VT_alpha_rel_CaMK].name = "alpha_rel_CaMK";
  P[VT_tau_h_LCaMK].name    = "tau_h_LCaMK";
  P[VT_RToverF].name        = "RToverF";

  P[VT_R].readFromFile              = false;
  P[VT_T].readFromFile              = false;
  P[VT_F].readFromFile              = false;
  P[VT_L].readFromFile              = false;
  P[VT_r].readFromFile              = false;
  P[VT_tau_h_L].readFromFile        = false;
  P[VT_A_h_fast].readFromFile       = false;
  P[VT_A_h_slow].readFromFile       = false;
  P[VT_v_cell].readFromFile         = false;
  P[VT_A_geo].readFromFile          = false;
  P[VT_A_cap].readFromFile          = false;
  P[VT_v_myo].readFromFile          = false;
  P[VT_v_nsr].readFromFile          = false;
  P[VT_v_jsr].readFromFile          = false;
  P[VT_v_ss].readFromFile           = false;
  P[VT_P_CaNa].readFromFile         = false;
  P[VT_P_CaK].readFromFile          = false;
  P[VT_P_Ca_CaMK].readFromFile      = false;
  P[VT_P_CaNa_CaMK].readFromFile    = false;
  P[VT_P_CaK_CaMK].readFromFile     = false;
  P[VT_A_f_fast].readFromFile       = false;
  P[VT_A_f_slow].readFromFile       = false;
  P[VT_k_Na1].readFromFile          = false;
  P[VT_k_Na2].readFromFile          = false;
  P[VT_k_asymm].readFromFile        = false;
  P[VT_k_Ca_on].readFromFile        = false;
  P[VT_k_Ca_off].readFromFile       = false;
  P[VT_k_m_1].readFromFile          = false;
  P[VT_k_p_2].readFromFile          = false;
  P[VT_k_p_4].readFromFile          = false;
  P[VT_MgADP].readFromFile          = false;
  P[VT_MgATP].readFromFile          = false;
  P[VT_K_MgATP].readFromFile        = false;
  P[VT_beta_tau].readFromFile       = false;
  P[VT_P_Ca].readFromFile           = false;
  P[VT_h_10].readFromFile           = false;
  P[VT_h_11].readFromFile           = false;
  P[VT_h_12].readFromFile           = false;
  P[VT_k_1].readFromFile            = false;
  P[VT_k_2].readFromFile            = false;
  P[VT_k_5].readFromFile            = false;
  P[VT_beta_1].readFromFile         = false;
  P[VT_alpha_2].readFromFile        = false;
  P[VT_alpha_4].readFromFile        = false;
  P[VT_alpha_rel].readFromFile      = false;
  P[VT_beta_tau_CaMK].readFromFile  = false;
  P[VT_alpha_rel_CaMK].readFromFile = false;
  P[VT_tau_h_LCaMK].readFromFile    = false;
  P[VT_RToverF].readFromFile        = false;

  ParameterLoader EPL(initFile, EMT_OHaraRudy);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // OHaraRudyParameters::Init

void OHaraRudyParameters::Calculate() {
#if KADEBUG
  cerr << "#OHaraRudyParameters - Calculate ..." << endl;
#endif  // if KADEBUG

  /// constants
  P[VT_R].value        = 8314.0;
  P[VT_T].value        = 310.0;
  P[VT_F].value        = 96485.0;
  P[VT_L].value        = 0.01;
  P[VT_r].value        = 0.0011;
  P[VT_tau_h_L].value  = 200.0;
  P[VT_A_h_fast].value = 0.99;
  P[VT_A_h_slow].value = 0.01;
  P[VT_A_f_fast].value = 0.6;
  P[VT_k_Na1].value    = 15.0;
  P[VT_k_Na2].value    = 5.0;
  P[VT_k_asymm].value  = 12.5;
  P[VT_k_Ca_on].value  = 1.5e6;
  P[VT_k_Ca_off].value = 5.0e3;
  P[VT_k_m_1].value    = 182.4;
  P[VT_k_p_2].value    = 687.2;
  P[VT_k_p_4].value    = 639.0;
  P[VT_MgADP].value    = 0.05;
  P[VT_MgATP].value    = 9.8;
  P[VT_K_MgATP].value  = 1.698e-7;
  P[VT_beta_tau].value = 4.75;
  P[VT_P_Ca].value     = P[VT_PCa_Multiplier].value * 0.0001;

  /// cell geometry
  P[VT_v_cell].value =  1000.0 * 3.14 * P[VT_r].value * P[VT_r].value * P[VT_L].value;
  P[VT_A_geo].value  =  2.0 * 3.14 * P[VT_r].value * P[VT_r].value + 2.0 * 3.14 * P[VT_r].value * P[VT_L].value;
  P[VT_A_cap].value  =  2.0 * P[VT_A_geo].value;
  P[VT_v_myo].value  =  0.68 * P[VT_v_cell].value;
  P[VT_v_nsr].value  =  0.0552 * P[VT_v_cell].value;
  P[VT_v_jsr].value  =  0.0048 * P[VT_v_cell].value;
  P[VT_v_ss].value   =  0.02 * P[VT_v_cell].value;

  /// permeability to ion X
  P[VT_P_CaNa].value      = 0.00125 * P[VT_P_Ca].value;
  P[VT_P_CaK].value       = 3.574e-04 * P[VT_P_Ca].value;
  P[VT_P_Ca_CaMK].value   = 1.1 * P[VT_P_Ca].value;
  P[VT_P_CaNa_CaMK].value = 0.00125 * P[VT_P_Ca_CaMK].value;
  P[VT_P_CaK_CaMK].value  = 3.5743e-04 * P[VT_P_Ca_CaMK].value;

  /// fraction of channels with gate X undergoing fast/slow process
  P[VT_A_f_slow].value = 1.0 - P[VT_A_f_fast].value;

  /// ??
  P[VT_h_10].value = P[VT_k_asymm].value + 1.0 + (P[VT_Na_o].value * (1.0 + P[VT_Na_o].value / P[VT_k_Na2].value)) /
    P[VT_k_Na1].value;
  P[VT_h_11].value = (P[VT_Na_o].value * P[VT_Na_o].value) / (P[VT_h_10].value * P[VT_k_Na1].value * P[VT_k_Na2].value);
  P[VT_h_12].value = 1.0 / P[VT_h_10].value;
  P[VT_k_1].value  = P[VT_h_12].value * P[VT_Ca_o].value * P[VT_k_Ca_on].value;
  P[VT_k_2].value  = P[VT_k_Ca_off].value;
  P[VT_k_5].value  = P[VT_k_Ca_off].value;

  // (de)phosphorylation rates
  P[VT_beta_1].value  = P[VT_k_m_1].value * P[VT_MgADP].value;
  P[VT_alpha_2].value = P[VT_k_p_2].value;
  P[VT_alpha_4].value = (P[VT_k_p_4].value * (P[VT_MgATP].value / P[VT_K_MgATP].value)) /
    (1.0 + (P[VT_MgATP].value / P[VT_K_MgATP].value));
  P[VT_alpha_rel].value      = 0.5 * P[VT_beta_tau].value;
  P[VT_beta_tau_CaMK].value  = 1.25 * P[VT_beta_tau].value;
  P[VT_alpha_rel_CaMK].value = 0.5 * P[VT_beta_tau_CaMK].value;

  // time constant of gate X
  P[VT_tau_h_LCaMK].value = 3.0 * P[VT_tau_h_L].value;

  P[VT_RToverF].value = (P[VT_R].value * P[VT_T].value) / P[VT_F].value;
}  // OHaraRudyParameters::Calculate

void OHaraRudyParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#OHaraRudyParameters - InitTable ..." << endl;
#endif  // if KADEBUG

  tinc *= 1000.0;  // second to millisecond conversion

  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    // V in mV
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);

    //// time dependent gates
    //// I_Na
    // m_inf[Vi] = 1.0 / (1.0 + exp(-(V + 39.57)/9.871));
    // const ML_CalcType tau_m = 1.0 / (6.765 * exp((V + 11.64)/34.77) + 8.552 * exp(-(V + 77.42)/5.955));
    // exptau_m[Vi] = exp(- tinc / tau_m);
    // h_inf[Vi] = 1.0 / (1.0 + exp((V + 82.9)/6.086));
    // const ML_CalcType tau_h_fast = 1.0 / (1.432e-05 * exp(-(V + 1.196)/6.285) + 6.149 * exp((V + 0.5096)/20.27));
    // const ML_CalcType tau_h_slow = 1.0 / (0.009746 * exp(-(V + 17.95)/28.05) + 0.3343 * exp((V + 5.730)/56.66));
    // exptau_h_fast[Vi] = exp(- tinc / tau_h_fast);
    // exptau_h_slow[Vi] = exp(- tinc / tau_h_slow);
    // j_inf[Vi] = h_inf[Vi];
    // const ML_CalcType tau_j = 2.038 + 1.0 / (0.02136 * exp(-(V + 100.6)/8.281) + 0.3052 * exp((V + 0.9941)/38.45));
    // exptau_j[Vi] = exp(- tinc / tau_j);
    // h_CaMK_inf[Vi] = 1.0 / (1.0 + exp((V + 89.1) / 6.086));
    // const ML_CalcType tau_h_CaMK_slow = 3.0 * tau_h_slow;
    // exptau_h_CaMK_slow[Vi] = exp(-tinc / tau_h_CaMK_slow);
    // j_CaMK_inf[Vi] = j_inf[Vi];
    // const ML_CalcType tau_j_CaMK = 1.46 * tau_j;
    // exptau_j_CaMK[Vi] = exp(-tinc / tau_j_CaMK);
    // m_L_inf[Vi] = 1.0 / (1.0 + exp(-(V + 42.85) / 5.264));
    // const ML_CalcType tau_m_L = tau_m;
    // exptau_m_L[Vi] = exp(-tinc / tau_m_L);
    // h_L_inf[Vi] = 1.0 / (1.0 + exp((V + 87.61) / 7.488));
    // exptau_h_L[Vi] = exp(-tinc / P[VT_tau_h_L].value);
    // h_L_CaMK_inf[Vi] = 1.0 / (1.0 + exp((V + 93.81) / 7.488));
    // const ML_CalcType tau_h_L_CaMK = 3.0 * P[VT_tau_h_L].value;
    // exptau_h_L_CaMK[Vi] = exp(-tinc / tau_h_L_CaMK);
    //
    //// I_to
    // a_inf[Vi] = 1.0 / (1.0 + exp(-(V - 14.34) / 14.82));
    // const ML_CalcType tau_a = 1.0515 / (1.0/(1.2089*(1.0 + exp(-(V - 18.4099) / 29.3814))) + 3.5/(1.0 + exp((V +
    // 100.0) / 29.3814)));
    // exptau_a[Vi] = exp(-tinc / tau_a);
    // i_inf[Vi] = 1.0 / (1.0 + exp((V + 43.94) / 5.711));
    // const ML_CalcType delta_epi = P[VT_celltype].value == 1.0 ? (1.0 - 0.95 / (1.0 + exp((V + 70.0) / 5.0))) : 1.0;
    // const ML_CalcType tau_i_fast_base =  4.562 + 1.0 / (0.3933 * exp(-(V + 100.0) / 100.0) + 0.08004 * exp((V + 50.0)
    // / 16.59));
    // const ML_CalcType tau_i_slow_base = 23.62 + 1.0 / (0.001416 * exp(-(V + 96.52) / 59.05) + 1.7808e-08 * exp((V +
    // 114.1) / 8.079));
    // const ML_CalcType tau_i_fast = tau_i_fast_base * delta_epi;
    // const ML_CalcType tau_i_slow = tau_i_slow_base * delta_epi;
    // exptau_i_fast[Vi] = exp(-tinc / tau_i_fast);
    // exptau_i_slow[Vi] = exp(-tinc / tau_i_slow);
    // A_i_fast[Vi] = 1.0 / (1.0 + exp((V - 213.6) / 151.2));
    // A_i_slow[Vi] = 1.0 - A_i_fast[Vi];
    // a_CaMK_inf[Vi] = 1.0 / (1.0 + exp(-(V - 24.34) / 14.82));
    // const ML_CalcType tau_a_CaMK = tau_a;
    // exptau_a_CaMK[Vi] = exp(-tinc / tau_a_CaMK);
    // i_CaMK_inf[Vi] = i_inf[Vi];
    // const ML_CalcType delta_CaMK_develop = 1.354 + 0.0001 / (exp((V - 167.4) / 15.89) + exp(-(V - 12.23) / 0.2154));
    // const ML_CalcType delta_CaMK_recover = 1.0 - 0.5 / (1.0 + exp((V + 70.0) / 20.0));
    // const ML_CalcType tau_i_CaMK_fast = tau_i_fast * delta_CaMK_develop * delta_CaMK_recover;
    // const ML_CalcType tau_i_CaMK_slow = tau_i_slow * delta_CaMK_develop * delta_CaMK_recover;
    // exptau_i_CaMK_fast[Vi] = exp(-tinc / tau_i_CaMK_fast);
    // exptau_i_CaMK_slow[Vi] = exp(-tinc / tau_i_CaMK_slow);
    //
    //
    //// I_CaL
    // d_inf[Vi] = 1.0 / (1.0 + exp(-(V + 3.940) / 4.230));
    // const ML_CalcType tau_d = 0.6 + 1.0 / (exp(-0.05 * (V + 6.0)) + exp(0.09 * (V + 14.0)));
    // exptau_d[Vi] = exp(-tinc / tau_d);
    // f_inf[Vi] = 1.0 / (1.0 + exp((V + 19.58) / 3.696));
    // const ML_CalcType tau_f_fast = 7.0 + 1.0 / (0.0045 * exp(-(V + 20.0) / 10.0) + 0.0045 * exp((V + 20.0) / 10.0));
    // const ML_CalcType tau_f_slow = 1000.0 + 1.0 / (3.5e-05 * exp(-(V + 5.0) / 4.0) + 3.5e-05 * exp((V + 5.0) / 6.0));
    // exptau_f_fast[Vi] = exp(-tinc / tau_f_fast);
    // exptau_f_slow[Vi] = exp(-tinc / tau_f_slow);
    // f_Ca_inf[Vi] = f_inf[Vi];
    // const ML_CalcType tau_f_Ca_fast = 7.0 + 1.0 / (0.04 * exp(-(V - 4.0) / 7.0) + 0.04 * exp((V - 4.0) / 7.0));
    // const ML_CalcType tau_f_Ca_slow = 100.0 + 1.0 / (0.00012 * exp(- V / 3.0) + 0.00012 * exp(V/7.0));
    // exptau_f_Ca_fast[Vi] = exp(-tinc / tau_f_Ca_fast);
    // exptau_f_Ca_slow[Vi] = exp(-tinc / tau_f_Ca_slow);
    // A_f_Ca_fast[Vi] = 0.3 + 0.6 / (1.0 + exp((V - 10.0) / 10.0));
    // A_f_Ca_slow[Vi] = 1.0 - A_f_Ca_fast[Vi];
    // j_Ca_inf[Vi] = f_Ca_inf[Vi];
    // exptau_j_Ca[Vi] = exp(-tinc / P[VT_tau_j_Ca].value);
    // f_CaMK_inf[Vi] = f_inf[Vi];
    // const ML_CalcType tau_f_CaMK_fast = 2.5 * tau_f_fast;
    // exptau_f_CaMK_fast[Vi] = exp(-tinc / tau_f_CaMK_fast);
    // f_Ca_CaMK_inf[Vi] = f_inf[Vi];
    // const ML_CalcType tau_f_Ca_CaMK_fast = 2.5 * tau_f_Ca_fast;
    // exptau_f_Ca_CaMK_fast[Vi] = exp(-tinc / tau_f_Ca_CaMK_fast);
    //
    //// I_Kr
    // x_r_inf[Vi] = 1.0 / (1.0 + exp(-(V + 8.337) / 6.789));
    // const ML_CalcType tau_x_r_fast = 12.98 + 1.0 / (0.3652 * exp((V - 31.66) / 3.869) + 4.123e-05 * exp(-(V - 47.78)
    // / 20.38));
    // const ML_CalcType tau_x_r_slow = 1.865 + 1.0 / (0.06629 * exp((V - 34.7) / 7.355) + 1.128e-05 * exp(- (V - 29.74)
    // / 25.94));
    // exptau_x_r_fast[Vi] = exp(-tinc / tau_x_r_fast);
    // exptau_x_r_slow[Vi] = exp(-tinc / tau_x_r_slow);
    // A_x_r_fast[Vi] = 1.0 / (1.0 + exp((V + 54.81) / 38.21));
    // A_x_r_slow[Vi] = 1.0 - A_x_r_fast[Vi];
    // R_Kr[Vi] = 1.0 / ((1.0 + exp((V + 55.0) / 75.0)) * (1.0 + exp((V - 10.0) / 30.0)));
    //
    //// I_Ks
    // x_s1_inf[Vi] = 1.0/(1.0 + exp(- (V + 11.6) / 8.932));
    // const ML_CalcType tau_x_s1 = 817.3 + 1.0 / (0.0002326 * exp((V + 48.28) / 17.80)+ 0.001292 * exp(- (V + 210.0) /
    // 230.0));
    // exptau_x_s1[Vi] = exp(-tinc / tau_x_s1);
    // x_s2_inf[Vi] = x_s1_inf[Vi];
    // const ML_CalcType tau_x_s2 = 1.0 / (0.01 * exp((V - 50.0) / 20.0) + 0.0193 * exp(- (V + 66.54) / 31.0));
    // exptau_x_s2[Vi] = exp(-tinc / tau_x_s2);
    //
    //// I_K1
    // x_K1_inf[Vi] = 1.0 / (1.0 + exp(-(V + 2.5538 * P[VT_K_o].value + 144.59) / (1.5692 * P[VT_K_o].value + 3.8115)));
    // const ML_CalcType tau_x_K1 = 122.2 / (exp(-(V + 127.2) / 20.36) + exp((V + 236.8) / 69.33));
    // exptau_x_K1[Vi] = exp(-tinc / tau_x_K1);
    // R_K1[Vi] = 1.0 / (1.0 + exp((V + 105.8 - 2.6 * P[VT_K_o].value) / 9.493));
    //
    //// utility parameters
    // x_Kb[Vi] = 1.0 / (1.0 + exp(-(V - 14.48) / (18.34)));
  }
}  // OHaraRudyParameters::InitTable
