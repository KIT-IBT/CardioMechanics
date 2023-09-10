/*      File: HimenoParameters.h        */

#ifndef HimenoPARAMETERS_H
#define HimenoPARAMETERS_H


#include <ParameterLoader.h>

namespace NS_HimenoParameters {
enum varType {
  VT_R = vtFirst,
  VT_Tx,
  VT_F,
  VT_stim_duration,
  VT_Amp,
  VT_k_on_CaM,
  VT_B_tot_CaM,
  VT_k_off_CaM,
  VT_k_on_TnCh,
  VT_B_tot_TnCh,
  VT_k_off_TnCh,
  VT_k_on_SR,
  VT_B_tot_SR,
  VT_k_off_SR,
  VT_B_tot_L_iz,
  VT_k_off_L_iz,
  VT_k_on_L_iz,
  VT_B_tot_H_iz,
  VT_k_off_H_iz,
  VT_k_on_H_iz,
  VT_B_tot_L_jnc,
  VT_Ca_2_jnc,
  VT_k_off_L_jnc,
  VT_k_on_L_jnc,
  VT_B_tot_H_jnc,
  VT_k_off_H_jnc,
  VT_k_on_H_jnc,
  VT_B_tot_CSQN,
  VT_k_off_CSQN,
  VT_k_on_CSQN,
  VT_G_dCa_jnciz,
  VT_Sc_Cell,
  VT_G_dCa_izblk,
  VT_P_trans,
  VT_T_L,
  VT_K_L,
  VT_ATP,
  VT_P_CaL_Ca,
  VT_Cao,
  VT_Nao,
  VT_Ko,
  VT_f_CaL_blk,
  VT_f_CaL_iz,
  VT_f_CaL_jnc,
  VT_k_I2O,
  VT_k_I1I2,
  VT_k_I1O,
  VT_f_LSM,
  VT_P_Na,
  VT_Mg_2_cyt,
  VT_SPM,
  VT_G_K1,
  VT_f_mode1,
  VT_G_Kr,
  VT_f_Ks_iz,
  VT_f_Ks_blk,
  VT_P_Ks_K,
  VT_G_Kto,
  VT_P_Kpl,
  VT_f_Cab_blk,
  VT_P_Cab,
  VT_f_Cab_iz,
  VT_P_bNSC_K,
  VT_P_bNSC_Na,
  VT_f_l_Ca_blk,
  VT_P_l_Ca_Na,
  VT_f_l_Ca_iz,
  VT_ATP_cyt,
  VT_G_KATP,
  VT_delta_Nai,
  VT_K_d_Nai_0,
  VT_delta_Nao,
  VT_K_d_Nao_0,
  VT_delta_Ki,
  VT_K_d_Ki_0,
  VT_delta_Ko,
  VT_K_d_Ko_0,
  VT_K_d_MgATP,
  VT_MgATP_cyt,
  VT_k_1_plus,
  VT_k_1_minus,
  VT_k_2_plus,
  VT_k_2_minus,
  VT_k_3_plus,
  VT_k_3_minus,
  VT_k_4_plus,
  VT_k_4_minus,
  VT_Pi,
  VT_H,
  VT_Amp_NaK,
  VT_K_m_act,
  VT_alpha_1_on,
  VT_alpha_1_off,
  VT_alpha_2_on,
  VT_alpha_2_off,
  VT_beta_1_on,
  VT_beta_1_off,
  VT_beta_2_on,
  VT_beta_2_off,
  VT_K_m_Nai,
  VT_K_m_Nao,
  VT_K_m_Cai,
  VT_K_m_Cao,
  VT_k_3,
  VT_k_4,
  VT_Amp_NCX,
  VT_f_NCX_blk,
  VT_f_NCX_iz,
  VT_K_m,
  VT_Amp_PMCA,
  VT_f_PMCA_blk,
  VT_f_PMCA_iz,
  VT_J_L,
  VT_J_R,
  VT_g_D,
  VT_Q_10,
  VT_sloc0,
  VT_f_n,
  VT_N_RyR,
  VT_p_O_RyR_base,
  VT_P_RyR,
  VT_K_dCai,
  VT_K_dCasr,
  VT_MgADP_cyt,
  VT_Amp_SERCA,
  VT_I_Kto_Na,
  VT_I_app,
  VT_halfSL,
  VT_TS_tot,
  VT_propFh,
  VT_Za,
  VT_Yv,
  VT_Yd,
  VT_Yc,
  VT_Lc,
  VT_Zb,
  VT_Yb,
  VT_rate_f,
  VT_convertF,
  VT_eqvhalfSL,
  VT_Zp,
  VT_Yp,
  VT_Zr,
  VT_Yr,
  VT_Zq,
  VT_Yq,
  VT_hwr,
  VT_rate_B,
  VT_hpr,
  VT_tau_i_f,
  VT_V_init,
  VT_init_CaMCa,
  VT_init_SRCa,
  VT_init_TnChCa,
  VT_init_p_O_NaT,
  VT_init_p_I_2_NaT,
  VT_init_p_I_s_NaT,
  VT_init_p_O_NaL,
  VT_init_p_I_1_NaL,
  VT_init_p_I_2_NaL,
  VT_init_p_I_s_NaL,
  VT_init_chi_r_fast,
  VT_init_chi_r_slow,
  VT_init_para_Xs1,
  VT_init_para_Xs2,
  VT_init_i_fast,
  VT_init_i_slow,
  VT_init_P_7,
  VT_init_P_8_13,
  VT_init_P_1_6,
  VT_init_p_E_1_NCX_blk,
  VT_init_p_I_2_NCX_blk,
  VT_init_p_I_1_NCX_blk,
  VT_init_p_E_1_NCX_iz,
  VT_init_p_I_1_NCX_iz,
  VT_init_p_I_2_NCX_iz,
  VT_init_Y_ooo,
  VT_init_Y_ooc,
  VT_init_Y_occ,
  VT_init_Y_coc,
  VT_init_Y_coo,
  VT_init_Y_cco,
  VT_init_Y_oco,
  VT_init_Y_co_iz,
  VT_init_Y_oo_iz,
  VT_init_Y_oc_iz,
  VT_init_Y_co_blk,
  VT_init_Y_oo_blk,
  VT_init_Y_oc_blk,
  VT_init_Ca_2_tot_jnc,
  VT_init_Ca_2_tot_iz,
  VT_init_Ca_2_tot_blk,
  VT_init_Ca_2_SRup,
  VT_init_Ca_2_tot_SRrl,
  VT_init_Nai,
  VT_init_Ki,
  VT_init_TSCa_3,
  VT_init_TSCa_3W,
  VT_init_TSCa_3S,
  VT_init_hw,
  VT_init_hp,
  VT_init_Pb_spm,
  VT_init_a,
  VT_init_L_bound_iz,
  VT_init_H_bound_iz,
  VT_init_TS_W,
  VT_init_TS_S,
  VT_inverseRTONF,
  VT_RTONF,
  VT_vol,
  VT_pi,
  VT_C,
  VT_K_dL_iz,
  VT_K_dH_iz,
  VT_K_dL_jnc,
  VT_K_dH_jnc,
  VT_K_d_CSQN_Ca,
  VT_P_CaL_Na,
  VT_P_CaL_K,
  VT_P_Ks_Na,
  VT_P_l_Ca_K,
  VT_alpha_2_plus,
  VT_f_L,
  VT_f_R,
  VT_k_oc,
  VT_V_jnc,
  VT_V_iz,
  VT_V_blk,
  VT_V_SRt,
  VT_V_SRrl,
  VT_V_SRup,
  VT_V_cyt,
  VT_Yvd,
  VT_T_L_K_L,
  VT_ATPfactor,
  VT_chi_Kr,
  VT_chi_K1,
  VT_chi_Kpl,
  VT_p_O_KATP,
  VT_chi_KATP,
  VT_MgATP_bar,
  VT_alpha_1_minus,
  VT_alpha_3_minus,
  VT_alpha_4_plus,
  VT_q_E_2_Na,
  VT_q_E_2_Ca,
  VT_alpha_1,
  vtLast
};
}  // namespace NS_HimenoParameters

using namespace NS_HimenoParameters;

class HimenoParameters : public vbNewElphyParameters {
 public:
  HimenoParameters(const char *, ML_CalcType);
  ~HimenoParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  ML_CalcType alpha_plus[RTDT];
  ML_CalcType alpha_minus[RTDT];
  ML_CalcType exp_dRTFVm[RTDT];
  ML_CalcType epsilon_minus[RTDT];
  ML_CalcType exp_VdRTF[RTDT];
  ML_CalcType exp_2VdRTF[RTDT];
  ML_CalcType f_C_Na[RTDT];
  ML_CalcType k_C2O[RTDT];
  ML_CalcType k_OC[RTDT];
  ML_CalcType k_OI2[RTDT];
  ML_CalcType k_C2I2[RTDT];
  ML_CalcType k_I2C[RTDT];
  ML_CalcType k_Isb[RTDT];
  ML_CalcType k_Isf[RTDT];
  ML_CalcType chi_r_infinity[RTDT];
  ML_CalcType tau_chi_r_fast[RTDT];
  ML_CalcType tau_chi_r_slow[RTDT];
  ML_CalcType A_chi_r_fast[RTDT];
  ML_CalcType A_chi_r_slow[RTDT];
  ML_CalcType R_Kr[RTDT];
  ML_CalcType para_Xs1_infinity[RTDT];
  ML_CalcType tau_Xs1[RTDT];
  ML_CalcType tau_Xs2[RTDT];
  ML_CalcType a_infinity[RTDT];
  ML_CalcType tau_a[RTDT];
  ML_CalcType i_infinity[RTDT];
  ML_CalcType tau_i_fast[RTDT];
  ML_CalcType tau_i_slow[RTDT];
  ML_CalcType A_i_fast[RTDT];
  ML_CalcType A_i_slow[RTDT];
  ML_CalcType p_O_Kpl[RTDT];
  ML_CalcType K_d_Nai[RTDT];
  ML_CalcType K_d_Nao[RTDT];
  ML_CalcType K_d_Ki[RTDT];
  ML_CalcType K_d_Ko[RTDT];
  ML_CalcType Nao_bar[RTDT];
  ML_CalcType Ko_bar[RTDT];
  ML_CalcType alpha_3_plus[RTDT];
  ML_CalcType alpha_2_minus[RTDT];
  ML_CalcType k_1[RTDT];
  ML_CalcType k_2[RTDT];
  ML_CalcType alpha_E[RTDT];
  ML_CalcType Ca_2_nd_L02s[RTDT];
  ML_CalcType Ca_2_nd_L0d[RTDT];
};  // class HimenoParameters
#endif  // ifndef HimenoPARAMETERS_H
