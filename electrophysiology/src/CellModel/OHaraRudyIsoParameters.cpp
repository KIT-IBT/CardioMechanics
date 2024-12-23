/*
 * File: OHaraRudyIsoParameters.cpp
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


#include <OHaraRudyIsoParameters.h>
OHaraRudyIsoParameters::OHaraRudyIsoParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

OHaraRudyIsoParameters::~OHaraRudyIsoParameters() {}

void OHaraRudyIsoParameters::PrintParameters() {
  cout << "OHaraRudyIsoParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void OHaraRudyIsoParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the OHaraRudyIso parameter from " << initFile << " ...\n";
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
  P[VT_Init_h_PKA_fast].name     = "Init_h_PKA_fast";
  P[VT_Init_h_PKA_slow].name     = "Init_h_PKA_slow";
  P[VT_Init_j_PKA].name          = "Init_j_PKA";
  P[VT_Init_h_both_fast].name    = "Init_h_both_fast";
  P[VT_Init_h_both_slow].name    = "Init_h_both_slow";
  P[VT_Init_j_both].name         = "Init_j_both";
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
  P[VT_Init_d_PKA].name          = "Init_d_PKA";
  P[VT_Init_f_PKA_fast].name     = "Init_f_PKA_fast";
  P[VT_Init_f_PKA_slow].name     = "Init_f_PKA_slow";
  P[VT_Init_f_Ca_PKA_fast].name  = "Init_f_Ca_PKA_fast";
  P[VT_Init_f_Ca_PKA_slow].name  = "Init_f_Ca_PKA_slow";
  P[VT_Init_f_both_fast].name    = "Init_f_both_fast";
  P[VT_Init_f_Ca_both_fast].name = "Init_f_Ca_both_fast";
  P[VT_Init_x_r_fast].name       = "Init_x_r_fast";
  P[VT_Init_x_r_slow].name       = "Init_x_r_slow";
  P[VT_Init_x_s1].name           = "Init_x_s1";
  P[VT_Init_x_s2].name           = "Init_x_s2";
  P[VT_Init_x_s1_PKA].name       = "Init_x_s1_PKA";
  P[VT_Init_x_s2_PKA].name       = "Init_x_s2_PKA";
  P[VT_Init_x_K1].name           = "Init_x_K1";
  P[VT_Init_J_rel_NP].name       = "Init_J_rel_NP";
  P[VT_Init_J_rel_CaMK].name     = "Init_J_rel_CaMK";
  P[VT_Init_J_rel_PKA].name     = "Init_J_rel_PKA";
  P[VT_Init_J_rel_both].name     = "Init_J_rel_both";
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
    
  /// parameters for isoprenaline implementation
  P[VT_ICaL_P_frac].name    = "ICaL_P_frac";
  P[VT_IKs_P_frac].name     = "IKs_P_frac";
  P[VT_IKb_P_frac].name     = "IKb_P_frac";
  P[VT_INaK_P_frac].name    = "INaK_P_frac";
  P[VT_INa_P_frac].name     = "INa_P_frac";
  P[VT_RyR_P_frac].name     = "RyR_P_frac";
  P[VT_SERCA_P_frac].name   = "SERCA_P_frac";
  P[VT_TnI_P_frac].name     = "TnI_P_frac";


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
  P[VT_P_Ca_PKA].name       = "P_Ca_PKA";
  P[VT_P_CaNa_PKA].name     = "P_CaNa_PKA";
  P[VT_P_CaK_PKA].name      = "P_CaK_PKA";
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
  P[VT_alpha_rel_PKA].name  = "alpha_rel_PKA";
  P[VT_alpha_rel_both].name = "alpha_rel_both";
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
  P[VT_P_Ca_PKA].readFromFile       = false;
  P[VT_P_CaNa_PKA].readFromFile       = false;
  P[VT_P_CaK_PKA].readFromFile       = false;
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
  P[VT_alpha_rel_PKA].readFromFile  = false;
  P[VT_alpha_rel_both].readFromFile = false;
  P[VT_tau_h_LCaMK].readFromFile    = false;
  P[VT_RToverF].readFromFile        = false;

  ParameterLoader EPL(initFile, EMT_OHaraRudyIso);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // OHaraRudyIsoParameters::Init

void OHaraRudyIsoParameters::Calculate() {
#if KADEBUG
  cerr << "#OHaraRudyIsoParameters - Calculate ..." << endl;
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
  P[VT_P_Ca_PKA].value    = 1.9 * P[VT_P_Ca].value;
  P[VT_P_CaNa_PKA].value  = 0.00125 * P[VT_P_Ca_PKA].value;
  P[VT_P_CaK_PKA].value   = 3.5743e-04 * P[VT_P_Ca_PKA].value;

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
  P[VT_alpha_rel_PKA].value = 2.5 * P[VT_alpha_rel].value;
  P[VT_alpha_rel_both].value = 2.5 * P[VT_alpha_rel_CaMK].value;

  // time constant of gate X
  P[VT_tau_h_LCaMK].value = 3.0 * P[VT_tau_h_L].value;

  P[VT_RToverF].value = (P[VT_R].value * P[VT_T].value) / P[VT_F].value;
  
}  // OHaraRudyIsoParameters::Calculate

void OHaraRudyIsoParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#OHaraRudyIsoParameters - InitTable ..." << endl;
#endif  // if KADEBUG

  tinc *= 1000.0;  // second to millisecond conversion

  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    // V in mV
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);

  }
}  // OHaraRudyIsoParameters::InitTable
