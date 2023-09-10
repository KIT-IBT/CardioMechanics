/*      File: HimenoParameters.cpp */

#include <HimenoParameters.h>

HimenoParameters::HimenoParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

HimenoParameters::~HimenoParameters() {}

void HimenoParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"HimenoParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}

void HimenoParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "#Loading the HimenoEtAl parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  // Initialization of the Parameters ...
  P[VT_R].name                  = "R";
  P[VT_Tx].name                 = "Tx";
  P[VT_F].name                  = "F";
  P[VT_stim_duration].name      = "stim_duration";
  P[VT_Amp].name                = "Amp";
  P[VT_k_on_CaM].name           = "k_on_CaM";
  P[VT_B_tot_CaM].name          = "B_tot_CaM";
  P[VT_k_off_CaM].name          = "k_off_CaM";
  P[VT_k_on_TnCh].name          = "k_on_TnCh";
  P[VT_B_tot_TnCh].name         = "B_tot_TnCh";
  P[VT_k_off_TnCh].name         = "k_off_TnCh";
  P[VT_k_on_SR].name            = "k_on_SR";
  P[VT_B_tot_SR].name           = "B_tot_SR";
  P[VT_k_off_SR].name           = "k_off_SR";
  P[VT_B_tot_L_iz].name         = "B_tot_L_iz";
  P[VT_k_off_L_iz].name         = "k_off_L_iz";
  P[VT_k_on_L_iz].name          = "k_on_L_iz";
  P[VT_B_tot_H_iz].name         = "B_tot_H_iz";
  P[VT_k_off_H_iz].name         = "k_off_H_iz";
  P[VT_k_on_H_iz].name          = "k_on_H_iz";
  P[VT_B_tot_L_jnc].name        = "B_tot_L_jnc";
  P[VT_Ca_2_jnc].name           = "Ca_2_jnc";
  P[VT_k_off_L_jnc].name        = "k_off_L_jnc";
  P[VT_k_on_L_jnc].name         = "k_on_L_jnc";
  P[VT_B_tot_H_jnc].name        = "B_tot_H_jnc";
  P[VT_k_off_H_jnc].name        = "k_off_H_jnc";
  P[VT_k_on_H_jnc].name         = "k_on_H_jnc";
  P[VT_B_tot_CSQN].name         = "B_tot_CSQN";
  P[VT_k_off_CSQN].name         = "k_off_CSQN";
  P[VT_k_on_CSQN].name          = "k_on_CSQN";
  P[VT_G_dCa_jnciz].name        = "G_dCa_jnciz";
  P[VT_Sc_Cell].name            = "Sc_Cell";
  P[VT_G_dCa_izblk].name        = "G_dCa_izblk";
  P[VT_P_trans].name            = "P_trans";
  P[VT_T_L].name                = "T_L";
  P[VT_K_L].name                = "K_L";
  P[VT_ATP].name                = "ATP";
  P[VT_P_CaL_Ca].name           = "P_CaL_Ca";
  P[VT_Cao].name                = "Cao";
  P[VT_Nao].name                = "Nao";
  P[VT_Ko].name                 = "Ko";
  P[VT_f_CaL_blk].name          = "f_CaL_blk";
  P[VT_f_CaL_iz].name           = "f_CaL_iz";
  P[VT_f_CaL_jnc].name          = "f_CaL_jnc";
  P[VT_k_I2O].name              = "k_I2O";
  P[VT_k_I1I2].name             = "k_I1I2";
  P[VT_k_I1O].name              = "k_I1O";
  P[VT_f_LSM].name              = "f_LSM";
  P[VT_P_Na].name               = "P_Na";
  P[VT_Mg_2_cyt].name           = "Mg_2_cyt";
  P[VT_SPM].name                = "SPM";
  P[VT_G_K1].name               = "G_K1";
  P[VT_f_mode1].name            = "f_mode1";
  P[VT_G_Kr].name               = "G_Kr";
  P[VT_f_Ks_iz].name            = "f_Ks_iz";
  P[VT_f_Ks_blk].name           = "f_Ks_blk";
  P[VT_P_Ks_K].name             = "P_Ks_K";
  P[VT_G_Kto].name              = "G_Kto";
  P[VT_P_Kpl].name              = "P_Kpl";
  P[VT_f_Cab_blk].name          = "f_Cab_blk";
  P[VT_P_Cab].name              = "P_Cab";
  P[VT_f_Cab_iz].name           = "f_Cab_iz";
  P[VT_P_bNSC_K].name           = "P_bNSC_K";
  P[VT_P_bNSC_Na].name          = "P_bNSC_Na";
  P[VT_f_l_Ca_blk].name         = "f_l_Ca_blk";
  P[VT_P_l_Ca_Na].name          = "P_l_Ca_Na";
  P[VT_f_l_Ca_iz].name          = "f_l_Ca_iz";
  P[VT_ATP_cyt].name            = "ATP_cyt";
  P[VT_G_KATP].name             = "G_KATP";
  P[VT_delta_Nai].name          = "delta_Nai";
  P[VT_K_d_Nai_0].name          = "K_d_Nai_0";
  P[VT_delta_Nao].name          = "delta_Nao";
  P[VT_K_d_Nao_0].name          = "K_d_Nao_0";
  P[VT_delta_Ki].name           = "delta_Ki";
  P[VT_K_d_Ki_0].name           = "K_d_Ki_0";
  P[VT_delta_Ko].name           = "delta_Ko";
  P[VT_K_d_Ko_0].name           = "K_d_Ko_0";
  P[VT_K_d_MgATP].name          = "K_d_MgATP";
  P[VT_MgATP_cyt].name          = "MgATP_cyt";
  P[VT_k_1_plus].name           = "k_1_plus";
  P[VT_k_1_minus].name          = "k_1_minus";
  P[VT_k_2_plus].name           = "k_2_plus";
  P[VT_k_2_minus].name          = "k_2_minus";
  P[VT_k_3_plus].name           = "k_3_plus";
  P[VT_k_3_minus].name          = "k_3_minus";
  P[VT_k_4_plus].name           = "k_4_plus";
  P[VT_k_4_minus].name          = "k_4_minus";
  P[VT_Pi].name                 = "Pi";
  P[VT_H].name                  = "H";
  P[VT_Amp_NaK].name            = "Amp_NaK";
  P[VT_K_m_act].name            = "K_m_act";
  P[VT_alpha_1_on].name         = "alpha_1_on";
  P[VT_alpha_1_off].name        = "alpha_1_off";
  P[VT_alpha_2_on].name         = "alpha_2_on";
  P[VT_alpha_2_off].name        = "alpha_2_off";
  P[VT_beta_1_on].name          = "beta_1_on";
  P[VT_beta_1_off].name         = "beta_1_off";
  P[VT_beta_2_on].name          = "beta_2_on";
  P[VT_beta_2_off].name         = "beta_2_off";
  P[VT_K_m_Nai].name            = "K_m_Nai";
  P[VT_K_m_Nao].name            = "K_m_Nao";
  P[VT_K_m_Cai].name            = "K_m_Cai";
  P[VT_K_m_Cao].name            = "K_m_Cao";
  P[VT_k_3].name                = "k_3";
  P[VT_k_4].name                = "k_4";
  P[VT_Amp_NCX].name            = "Amp_NCX";
  P[VT_f_NCX_blk].name          = "f_NCX_blk";
  P[VT_f_NCX_iz].name           = "f_NCX_iz";
  P[VT_K_m].name                = "K_m";
  P[VT_Amp_PMCA].name           = "Amp_PMCA";
  P[VT_f_PMCA_blk].name         = "f_PMCA_blk";
  P[VT_f_PMCA_iz].name          = "f_PMCA_iz";
  P[VT_J_L].name                = "J_L";
  P[VT_J_R].name                = "J_R";
  P[VT_g_D].name                = "g_D";
  P[VT_Q_10].name               = "Q_10";
  P[VT_sloc0].name              = "sloc0";
  P[VT_f_n].name                = "f_n";
  P[VT_N_RyR].name              = "N_RyR";
  P[VT_p_O_RyR_base].name       = "p_O_RyR_base";
  P[VT_P_RyR].name              = "P_RyR";
  P[VT_K_dCai].name             = "K_dCai";
  P[VT_K_dCasr].name            = "K_dCasr";
  P[VT_MgADP_cyt].name          = "MgADP_cyt";
  P[VT_Amp_SERCA].name          = "Amp_SERCA";
  P[VT_I_Kto_Na].name           = "I_Kto_Na";
  P[VT_I_app].name              = "I_app";
  P[VT_halfSL].name             = "halfSL";
  P[VT_TS_tot].name             = "TS_tot";
  P[VT_propFh].name             = "propFh";
  P[VT_Za].name                 = "Za";
  P[VT_Yv].name                 = "Yv";
  P[VT_Yd].name                 = "Yd";
  P[VT_Yc].name                 = "Yc";
  P[VT_Lc].name                 = "Lc";
  P[VT_Zb].name                 = "Zb";
  P[VT_Yb].name                 = "Yb";
  P[VT_rate_f].name             = "rate_f";
  P[VT_convertF].name           = "convertF";
  P[VT_eqvhalfSL].name          = "eqvhalfSL";
  P[VT_Zp].name                 = "Zp";
  P[VT_Yp].name                 = "Yp";
  P[VT_Zr].name                 = "Zr";
  P[VT_Yr].name                 = "Yr";
  P[VT_Zq].name                 = "Zq";
  P[VT_Yq].name                 = "Yq";
  P[VT_hwr].name                = "hwr";
  P[VT_rate_B].name             = "rate_B";
  P[VT_hpr].name                = "hpr";
  P[VT_tau_i_f].name            = "tau_i_f";
  P[VT_V_init].name             = "init_V_m";
  P[VT_init_CaMCa].name         = "init_CaMCa";
  P[VT_init_SRCa].name          = "init_SRCa";
  P[VT_init_TnChCa].name        = "init_TnChCa";
  P[VT_init_p_O_NaT].name       = "init_p_O_NaT";
  P[VT_init_p_I_2_NaT].name     = "init_p_I_2_NaT";
  P[VT_init_p_I_s_NaT].name     = "init_p_I_s_NaT";
  P[VT_init_p_O_NaL].name       = "init_p_O_NaL";
  P[VT_init_p_I_1_NaL].name     = "init_p_I_1_NaL";
  P[VT_init_p_I_2_NaL].name     = "init_p_I_2_NaL";
  P[VT_init_p_I_s_NaL].name     = "init_p_I_s_NaL";
  P[VT_init_chi_r_fast].name    = "init_chi_r_fast";
  P[VT_init_chi_r_slow].name    = "init_chi_r_slow";
  P[VT_init_para_Xs1].name      = "init_para_Xs1";
  P[VT_init_para_Xs2].name      = "init_para_Xs2";
  P[VT_init_i_fast].name        = "init_i_fast";
  P[VT_init_i_slow].name        = "init_i_slow";
  P[VT_init_P_7].name           = "init_P_7";
  P[VT_init_P_8_13].name        = "init_P_8_13";
  P[VT_init_P_1_6].name         = "init_P_1_6";
  P[VT_init_p_E_1_NCX_blk].name = "init_p_E_1_NCX_blk";
  P[VT_init_p_I_2_NCX_blk].name = "init_p_I_2_NCX_blk";
  P[VT_init_p_I_1_NCX_blk].name = "init_p_I_1_NCX_blk";
  P[VT_init_p_E_1_NCX_iz].name  = "init_p_E_1_NCX_iz";
  P[VT_init_p_I_1_NCX_iz].name  = "init_p_I_1_NCX_iz";
  P[VT_init_p_I_2_NCX_iz].name  = "init_p_I_2_NCX_iz";
  P[VT_init_Y_ooo].name         = "init_Y_ooo";
  P[VT_init_Y_ooc].name         = "init_Y_ooc";
  P[VT_init_Y_occ].name         = "init_Y_occ";
  P[VT_init_Y_coc].name         = "init_Y_coc";
  P[VT_init_Y_coo].name         = "init_Y_coo";
  P[VT_init_Y_cco].name         = "init_Y_cco";
  P[VT_init_Y_oco].name         = "init_Y_oco";
  P[VT_init_Y_co_iz].name       = "init_Y_co_iz";
  P[VT_init_Y_oo_iz].name       = "init_Y_oo_iz";
  P[VT_init_Y_oc_iz].name       = "init_Y_oc_iz";
  P[VT_init_Y_co_blk].name      = "init_Y_co_blk";
  P[VT_init_Y_oo_blk].name      = "init_Y_oo_blk";
  P[VT_init_Y_oc_blk].name      = "init_Y_oc_blk";
  P[VT_init_Ca_2_tot_jnc].name  = "init_Ca_2_tot_jnc";
  P[VT_init_Ca_2_tot_iz].name   = "init_Ca_2_tot_iz";
  P[VT_init_Ca_2_tot_blk].name  = "init_Ca_2_tot_blk";
  P[VT_init_Ca_2_SRup].name     = "init_Ca_2_SRup";
  P[VT_init_Ca_2_tot_SRrl].name = "init_Ca_2_tot_SRrl";
  P[VT_init_Nai].name           = "init_Nai";
  P[VT_init_Ki].name            = "init_Ki";
  P[VT_init_TSCa_3].name        = "init_TSCa_3";
  P[VT_init_TSCa_3W].name       = "init_TSCa_3W";
  P[VT_init_TSCa_3S].name       = "init_TSCa_3S";
  P[VT_init_hw].name            = "init_hw";
  P[VT_init_hp].name            = "init_hp";
  P[VT_init_Pb_spm].name        = "init_Pb_spm";
  P[VT_init_a].name             = "init_a";
  P[VT_init_L_bound_iz].name    = "init_L_bound_iz";
  P[VT_init_H_bound_iz].name    = "init_H_bound_iz";
  P[VT_init_TS_S].name          = "init_TS_S";
  P[VT_init_TS_W].name          = "init_TS_W";
  P[VT_vol].name                = "vol";


  P[VT_pi].readFromFile            = false;
  P[VT_RTONF].readFromFile         = false;
  P[VT_inverseRTONF].readFromFile  = false;
  P[VT_vol].readFromFile           = false;
  P[VT_C].readFromFile             = false;
  P[VT_K_dL_iz].readFromFile       = false;
  P[VT_K_dH_iz].readFromFile       = false;
  P[VT_K_dL_jnc].readFromFile      = false;
  P[VT_K_dH_jnc].readFromFile      = false;
  P[VT_K_d_CSQN_Ca].readFromFile   = false;
  P[VT_P_CaL_Na].readFromFile      = false;
  P[VT_P_CaL_K].readFromFile       = false;
  P[VT_P_Ks_Na].readFromFile       = false;
  P[VT_P_l_Ca_K].readFromFile      = false;
  P[VT_alpha_2_plus].readFromFile  = false;
  P[VT_f_L].readFromFile           = false;
  P[VT_f_R].readFromFile           = false;
  P[VT_k_oc].readFromFile          = false;
  P[VT_V_jnc].readFromFile         = false;
  P[VT_V_iz].readFromFile          = false;
  P[VT_V_blk].readFromFile         = false;
  P[VT_V_SRt].readFromFile         = false;
  P[VT_V_SRrl].readFromFile        = false;
  P[VT_V_SRup].readFromFile        = false;
  P[VT_V_cyt].readFromFile         = false;
  P[VT_Yvd].readFromFile           = false;
  P[VT_T_L_K_L].readFromFile       = false;
  P[VT_ATPfactor].readFromFile     = false;
  P[VT_chi_Kr].readFromFile        = false;
  P[VT_chi_K1].readFromFile        = false;
  P[VT_chi_Kpl].readFromFile       = false;
  P[VT_p_O_KATP].readFromFile      = false;
  P[VT_chi_KATP].readFromFile      = false;
  P[VT_p_O_KATP].readFromFile      = false;
  P[VT_MgATP_bar].readFromFile     = false;
  P[VT_alpha_1_minus].readFromFile = false;
  P[VT_alpha_3_minus].readFromFile = false;
  P[VT_alpha_4_plus].readFromFile  = false;
  P[VT_q_E_2_Na].readFromFile      = false;
  P[VT_q_E_2_Ca].readFromFile      = false;
  P[VT_alpha_1].readFromFile       = false;

  ParameterLoader EPL(initFile, EMT_HimenoEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  // End Initialization of the Parameters ...

  Calculate();
  InitTable(tinc);
}  // HimenoParameters::Init

void HimenoParameters::Calculate() {
#if KADEBUG
  cerr << "#HimenoEtAlParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_pi].value            = M_PI;
  P[VT_RTONF].value         = (P[VT_R].value)*(P[VT_Tx].value)/(P[VT_F].value);
  P[VT_inverseRTONF].value  = 1/(P[VT_RTONF].value);
  P[VT_vol].value           = 120 * (37.62 * P[VT_Sc_Cell].value) * 8.4;
  P[VT_C].value             = 192.46 * P[VT_Sc_Cell].value;
  P[VT_K_dL_iz].value       = P[VT_k_off_L_iz].value / P[VT_k_on_L_iz].value;
  P[VT_K_dH_iz].value       = P[VT_k_off_H_iz].value / P[VT_k_on_H_iz].value;
  P[VT_K_dL_jnc].value      = P[VT_k_off_L_jnc].value / P[VT_k_on_L_jnc].value;
  P[VT_K_dH_jnc].value      = P[VT_k_off_H_jnc].value / P[VT_k_on_H_jnc].value;
  P[VT_K_d_CSQN_Ca].value   = P[VT_k_off_CSQN].value / P[VT_k_on_CSQN].value;
  P[VT_P_CaL_Na].value      = 0.0000185 * P[VT_P_CaL_Ca].value;
  P[VT_P_CaL_K].value       = 0.000367 * P[VT_P_CaL_Ca].value;
  P[VT_P_Ks_Na].value       = 0.04 * P[VT_P_Ks_K].value;
  P[VT_P_l_Ca_K].value      = P[VT_P_l_Ca_Na].value;
  P[VT_alpha_2_plus].value  = P[VT_k_2_plus].value;
  P[VT_f_L].value           = P[VT_J_L].value / P[VT_g_D].value;
  P[VT_f_R].value           = P[VT_J_R].value / P[VT_g_D].value;
  P[VT_k_oc].value          =  P[VT_Q_10].value * 0.5664;
  P[VT_V_jnc].value         = 0.008 * P[VT_vol].value;
  P[VT_V_iz].value          = 0.035 * P[VT_vol].value;
  P[VT_V_blk].value         = 0.68 * P[VT_vol].value;
  P[VT_V_SRt].value         = 0.06 * P[VT_vol].value;
  P[VT_V_SRrl].value        = 0.2 * P[VT_V_SRt].value;
  P[VT_V_SRup].value        = 0.8 * P[VT_V_SRt].value;
  P[VT_V_cyt].value         = P[VT_V_jnc].value + P[VT_V_iz].value + P[VT_V_blk].value;
  P[VT_Yvd].value           = P[VT_Yv].value;
  P[VT_T_L_K_L].value       = P[VT_T_L].value * P[VT_K_L].value;
  P[VT_ATPfactor].value     = 1 / (1 + pow(1.4 / P[VT_ATP].value, 3));
  P[VT_chi_Kr].value        = sqrt(P[VT_Ko].value / 4.5);
  P[VT_chi_K1].value        = pow(P[VT_Ko].value / 4.5, 0.4) / (1.0 + exp(-(P[VT_Ko].value - 2.2) / 0.6));
  P[VT_chi_Kpl].value       = pow(P[VT_Ko].value / 5.4, 0.16);
  P[VT_chi_KATP].value      = 0.0236 * pow(P[VT_Ko].value, 0.24);
  P[VT_p_O_KATP].value      = 0.8 / (1.0 + pow(P[VT_ATP_cyt].value / 0.1, 2));
  P[VT_MgATP_bar].value     = P[VT_MgATP_cyt].value / P[VT_K_d_MgATP].value;
  P[VT_alpha_1_minus].value = P[VT_k_1_minus].value * P[VT_MgADP_cyt].value;
  P[VT_alpha_3_minus].value = P[VT_k_3_minus].value * P[VT_Pi].value * P[VT_H].value / (1 + P[VT_MgATP_bar].value);
  P[VT_alpha_4_plus].value  = P[VT_k_4_plus].value * P[VT_MgATP_bar].value/ (1 + P[VT_MgATP_bar].value);
  P[VT_q_E_2_Na].value      = 1.0 /
    (1.0 + pow(P[VT_K_m_Nao].value / P[VT_Nao].value, 3) * (1.0 + P[VT_Cao].value / P[VT_K_m_Cao].value));
  P[VT_q_E_2_Ca].value = 1.0 /
    (1.0 + P[VT_K_m_Cao].value / P[VT_Cao].value * (1.0 + pow(P[VT_Nao].value / P[VT_K_m_Nao].value, 3)));
  P[VT_alpha_1].value = 25900 * P[VT_MgATP_cyt].value;
}  // HimenoParameters::Calculate

void HimenoParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#HimenoEtAlParameters - InitTable ..." << endl;
#endif  // if KADEBUG
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    alpha_plus[Vi]    = 1 / (3.734 * exp(-V / 8.5) + 0.35 * exp(-V / 3500));
    alpha_minus[Vi]   = 1 / (4.65 * exp(V / 15) + 1.363 * exp(V / 100));
    exp_dRTFVm[Vi]    = exp(-2* P[VT_inverseRTONF].value * V);
    epsilon_minus[Vi] = 1 / (8084 * exp(V / 10) + 158 * exp(V / 1000)) + 1 /
      (134736 * exp(-V / 5) + 337 * exp(-V / 2000));
    f_C_Na[Vi]            = 1 / (1 + exp(-(V + 48) / 7));
    k_C2O[Vi]             = 1 / (0.0025 * exp(V / -8.0) + 0.15 * exp(V / -100.0));
    k_OC[Vi]              = 1 / (30 * exp(V / 12.0) + 0.53 * exp(V / 50.0));
    k_OI2[Vi]             = 1 / (0.0433 * exp(V / -27.0) + 0.34 * exp(V / -2000.0));
    k_C2I2[Vi]            = 0.5 / (1.0 + (P[VT_k_I2O].value * k_OC[Vi]) / (k_OI2[Vi] * k_C2O[Vi]));
    k_I2C[Vi]             = 0.5 - k_C2I2[Vi];
    k_Isb[Vi]             = 1 / (300000.0 * exp(V / 10.0) + 50000 * exp(V / 16.0));
    k_Isf[Vi]             = 1 / (0.016 * exp(V / -9.9) + 8.0 * exp(V / -45.0));
    chi_r_infinity[Vi]    = 1 / (1 + exp(-(V + 8.337) / 6.789));
    tau_chi_r_fast[Vi]    = 12.98 + 1 / (0.3652 * exp((V - 31.66) / 3.869) + 0.00004123 * exp(-(V - 47.78) / 20.38));
    tau_chi_r_slow[Vi]    = 1.865 + 1 / (0.06629 * exp((V - 34.70) / 7.355) + 0.00001128 * exp(-(V - 29.74) / 25.94));
    A_chi_r_fast[Vi]      = 1 / (1 + exp((V + 4.81) / 38.21));
    A_chi_r_slow[Vi]      = 1 - A_chi_r_fast[Vi];
    R_Kr[Vi]              = 1 / ((1 + exp((V + 55) / 75)) * (1 + exp((V - 10) / 30)));
    para_Xs1_infinity[Vi] = 1 / (1 + exp(-(V + 11.60) / 8.932));
    tau_Xs1[Vi]           = 817.3 + 1 / (0.0002326 * exp((V + 48.28) / 17.80) + 0.001292 * exp(-(V + 210.0) / 230.0));
    tau_Xs2[Vi]           = 1 / (0.01 * exp((V - 50) / 20) + 0.0193 * exp(-(V + 66.54) / 31));
    a_infinity[Vi]        = 1 / (1 + exp(-(V - 14.34) / 14.82));
    tau_a[Vi]             = 1.0515 /
      (1 / (1.2089 * (1 + exp(-(V - 18.41) / 29.38))) + 3.5 / (1 + exp((V + 100) / 29.38)));
    i_infinity[Vi] = 1 / (1 + exp((V + 43.94) / 5.711));

    // tau_i_fast[Vi] = 4.562 + 1 / (0.3933 * exp(-(V + 100) / 100) + 0.08004 * exp((V + 50) / 16.59));
    // tau_i_slow[Vi] = 23.62 + 1 / (0.001416 * exp(-(V + 96.52) / 59.05) + 0.000000017808 * exp((V + 114.1) / 8.079));
    tau_i_fast[Vi] = (4.562 + 1 / (0.3933 * exp(-(V + 100) / 100) + 0.08004 * exp((V + 50) / 16.59)))*
      (1.0-(P[VT_tau_i_f].value/(1.0+exp((V +70.0)/5.0))));
    tau_i_slow[Vi] = (23.62 + 1 / (0.001416 * exp(-(V + 96.52) / 59.05) + 0.000000017808 * exp((V + 114.1) / 8.079)))*
      (1.0-(P[VT_tau_i_f].value/(1.0+exp((V +70.0)/5.0))));


    A_i_fast[Vi]     = 1 / (1 + exp((V - 213.6) / 151.2));
    A_i_slow[Vi]     = 1 - A_i_fast[Vi];
    p_O_Kpl[Vi]      = V / (1 - exp(-V / 13.0));
    K_d_Nai[Vi]      = P[VT_K_d_Nai_0].value * exp((P[VT_delta_Nai].value * V) / P[VT_RTONF].value);
    K_d_Nao[Vi]      = P[VT_K_d_Nao_0].value * exp((P[VT_delta_Nao].value * V) / P[VT_RTONF].value);
    K_d_Ki[Vi]       = P[VT_K_d_Ki_0].value * exp((P[VT_delta_Ki].value * V) / P[VT_RTONF].value);
    K_d_Ko[Vi]       = P[VT_K_d_Ko_0].value * exp((P[VT_delta_Ko].value * V) / P[VT_RTONF].value);
    Nao_bar[Vi]      = P[VT_Nao].value / K_d_Nao[Vi];
    Ko_bar[Vi]       = P[VT_Ko].value / K_d_Ko[Vi];
    alpha_3_plus[Vi] = P[VT_k_3_plus].value *
      pow(Ko_bar[Vi], 2) / (pow(1 + Nao_bar[Vi], 3) + pow(1 + Ko_bar[Vi], 2) - 1);
    alpha_2_minus[Vi] = P[VT_k_2_minus].value *
      pow(Nao_bar[Vi], 3) / (pow(1 + Nao_bar[Vi], 3) + pow(1 + Ko_bar[Vi], 2) - 1);
    k_1[Vi]          = exp((0.32 * V) / P[VT_RTONF].value);
    k_2[Vi]          = exp(((0.32 - 1) * V) / P[VT_RTONF].value);
    alpha_E[Vi]      = k_2[Vi] * P[VT_q_E_2_Na].value + P[VT_k_4].value * P[VT_q_E_2_Ca].value;
    Ca_2_nd_L02s[Vi] =  P[VT_f_L].value * 2* P[VT_inverseRTONF].value * V * exp_dRTFVm[Vi] / (1 - exp_dRTFVm[Vi]) *
      P[VT_Cao].value;
    Ca_2_nd_L0d[Vi] = (1 + P[VT_f_L].value * 2* P[VT_inverseRTONF].value * V / (1 - exp_dRTFVm[Vi]));
  }
}  // HimenoParameters::InitTable
