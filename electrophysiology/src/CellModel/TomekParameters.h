/*
 * File: TomekParameters.h
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


#ifndef TOMEKPARAMETERS_H
#define TOMEKPARAMETERS_H

/// Tomek model variations.
/// Stretch avtivated channel by Pueyo et al. 2016 Implemented by Albert Dasi
#define ISAC

/// Troponin dynamics from Land et. al tension model to use mechano-electric feedback. Implemented by Jan-Erik Duhme.
#define TRPN

#include <ParameterLoader.h>

namespace NS_TomekParameters {
enum varType {
  VT_amplitude = vtFirst,
  VT_duration,
  VT_celltype,
  VT_Na_o,
  VT_Ca_o,
  VT_K_o,
  VT_Cl_o,
  VT_INaF_Multiplier,
  VT_INaL_Multiplier,
  VT_Ito_Multiplier,
  VT_PCa_Multiplier,
  VT_IKr_Multiplier,
  VT_IKs_Multiplier,
  VT_IK1_Multiplier,
  VT_INaCa_Multiplier,
  VT_INaK_Multiplier,
  VT_IKb_Multiplier,
  VT_INab_Multiplier,
  VT_ICab_Multiplier,
  VT_IpCa_Multiplier,

  /// state variables
  VT_Init_Vm,
  VT_Init_Na_i,
  VT_Init_Na_ss,
  VT_Init_K_i,
  VT_Init_K_ss,
  VT_Init_Ca_i,
  VT_Init_Ca_ss,
  VT_Init_Ca_nsr,
  VT_Init_Ca_jsr,
  VT_Init_Cl_i,
  VT_Init_m,
  VT_Init_A_h,
  VT_Init_B_h,
  VT_Init_h,
  VT_Init_A_j,
  VT_Init_B_j,  
  VT_Init_j,
  VT_Init_h_p,
  VT_Init_j_p,
  VT_Init_m_L,
  VT_Init_h_L,
  VT_Init_h_L_CaMK,
  VT_Init_a,
  VT_Init_i_fast,
  VT_Init_i_slow,
  VT_Init_a_CaMK,
  VT_Init_i_CaMK_fast,
  VT_Init_i_CaMK_slow,
  VT_Init_d,
  VT_Init_f_fast,
  VT_Init_f_slow,
  VT_Init_f_Ca_fast,
  VT_Init_f_Ca_slow,
  VT_Init_j_Ca,
  VT_Init_n_ss,
  VT_Init_n_i,
  VT_Init_f_CaMK_fast,
  VT_Init_f_Ca_CaMK_fast,
  VT_Init_C_0,
  VT_Init_C_1,
  VT_Init_C_2,
  VT_Init_O,
  VT_Init_I,
  VT_Init_x_s1,
  VT_Init_x_s2,
  VT_Init_J_rel_NP,
  VT_Init_J_rel_CaMK,
  VT_Init_CaMK_trap,

  /// parameters that are not read from .ev
  VT_R,
  VT_T,
  VT_F,
  VT_L,
  VT_r,
  VT_tau_h_L,
  VT_A_h,
  VT_v_cell,
  VT_A_geo,
  VT_A_cap,
  VT_v_myo,
  VT_v_nsr,
  VT_v_jsr,
  VT_v_ss,
  VT_P_CaNa,
  VT_P_CaK,
  VT_P_Ca_CaMK,
  VT_P_CaNa_CaMK,
  VT_P_CaK_CaMK,
  VT_A_f_fast,
  VT_A_f_slow,
  VT_k_Na1,
  VT_k_Na2,
  VT_k_asymm,
  VT_k_Ca_on,
  VT_k_Ca_off,
  VT_k_m_1,
  VT_k_p_2,
  VT_k_p_4,
  VT_MgADP,
  VT_MgATP,
  VT_K_MgATP,
  VT_beta_tau,
  VT_P_Ca,
  VT_h_10,
  VT_h_11,
  VT_h_12,
  VT_k_1,
  VT_k_2,
  VT_k_5,
  VT_beta_1,
  VT_alpha_2,
  VT_alpha_4,
  VT_alpha_rel,
  VT_beta_tau_CaMK,
  VT_alpha_rel_CaMK,
  VT_tau_h_LCaMK,
  VT_RToverF,
  
  #ifdef ISAC
  VT_ISAC_SWITCH,
  VT_alpha,
  VT_G_sac,
  VT_K_sac,
  VT_E_sac,
#endif // ifdef ISAC
  #ifdef TRPN
  VT_Init_Ca_TRPN,
  VT_beta1,
  VT_Ca_T50,
  VT_k_TRPN,
  VT_n_TRPN,
#endif  // ifdef TRPN
  vtLast
};
}  // namespace NS_TomekParameters

using namespace NS_TomekParameters;

class TomekParameters : public vbNewElphyParameters {
 public:
  TomekParameters(const char *, ML_CalcType);
  ~TomekParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);

  ML_CalcType m_inf[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType h_inf[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType j_inf[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType h_p_inf[RTDT];
  ML_CalcType exptau_h_p[RTDT];
  ML_CalcType j_p_inf[RTDT];
  ML_CalcType exptau_j_p[RTDT];
  ML_CalcType m_L_inf[RTDT];
  ML_CalcType exptau_m_L[RTDT];
  ML_CalcType h_L_inf[RTDT];
  ML_CalcType exptau_h_L[RTDT];
  ML_CalcType h_L_CaMK_inf[RTDT];
  ML_CalcType exptau_h_L_CaMK[RTDT];
  ML_CalcType a_inf[RTDT];
  ML_CalcType exptau_a[RTDT];
  ML_CalcType i_inf[RTDT];
  ML_CalcType exptau_i_fast[RTDT];
  ML_CalcType exptau_i_slow[RTDT];
  ML_CalcType A_i_fast[RTDT];
  ML_CalcType A_i_slow[RTDT];
  ML_CalcType a_CaMK_inf[RTDT];
  ML_CalcType exptau_a_CaMK[RTDT];
  ML_CalcType i_CaMK_inf[RTDT];
  ML_CalcType exptau_i_CaMK_fast[RTDT];
  ML_CalcType exptau_i_CaMK_slow[RTDT];
  ML_CalcType d_inf[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType f_inf[RTDT];
  ML_CalcType exptau_f_fast[RTDT];
  ML_CalcType exptau_f_slow[RTDT];
  ML_CalcType f_Ca_inf[RTDT];
  ML_CalcType exptau_f_Ca_fast[RTDT];
  ML_CalcType exptau_f_Ca_slow[RTDT];
  ML_CalcType A_f_Ca_fast[RTDT];
  ML_CalcType A_f_Ca_slow[RTDT];
  ML_CalcType j_Ca_inf[RTDT];
  ML_CalcType exptau_j_Ca[RTDT];
  ML_CalcType f_CaMK_inf[RTDT];
  ML_CalcType exptau_f_CaMK_fast[RTDT];
  ML_CalcType f_Ca_CaMK_inf[RTDT];
  ML_CalcType exptau_f_Ca_CaMK_fast[RTDT];
  ML_CalcType x_r_inf[RTDT];
  ML_CalcType exptau_x_r_fast[RTDT];
  ML_CalcType exptau_x_r_slow[RTDT];
  ML_CalcType A_x_r_fast[RTDT];
  ML_CalcType A_x_r_slow[RTDT];
  ML_CalcType C_0[RTDT];
  ML_CalcType C_1[RTDT];
  ML_CalcType C_2[RTDT];
  ML_CalcType O[RTDT];
  ML_CalcType I[RTDT];
  ML_CalcType R_Kr[RTDT];
  ML_CalcType x_s1_inf[RTDT];
  ML_CalcType exptau_x_s1[RTDT];
  ML_CalcType x_s2_inf[RTDT];
  ML_CalcType exptau_x_s2[RTDT];
  ML_CalcType x_K1_inf[RTDT];
  ML_CalcType exptau_x_K1[RTDT];
  ML_CalcType R_K1[RTDT];
  ML_CalcType x_Kb[RTDT];
};  // class TomekParameters

#endif  // ifndef TOMEKPARAMETERS_H
