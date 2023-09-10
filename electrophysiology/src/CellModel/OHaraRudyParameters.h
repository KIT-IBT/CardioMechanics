/* -------------------------------------------------------

   OHaraRudyParameters.h

   Ver. 1.3.0

   Created:       Tobias Gerach (01.2018)
   Last modified: Tobias Gerach (07.05.2023)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#ifndef OHARARUDYPARAMETERS_H
#define OHARARUDYPARAMETERS_H

/// OHara Rudy model variations.
/// Stretch avtivated channel by Pueyo et al. 2016 Implemented by Albert Dasi
#define ISAC

/// Troponin dynamics from Land et al. tension model to use mechano-electric feedback. Implemented by Stephanie Appel
#define TRPN

/// Heart failure specific alterations from Gomez et al. 2014 Implemented by Albert Dasi
#define HF

/// Heart failure specific alterations from Bollen et al. 2017. Implemented by Albert Dasi
#define HFM

/// Changed I_Na h and j gate time constants following Dutta et al. 2017. Should be used together with PASSINI
/// Implemented by Tobias Gerach
#define DUTTA

/// Changed I_Na h and j gate following Passini et al. 2016. Implemented by Tobias Gerach
#define PASSINI

#include <ParameterLoader.h>

namespace NS_OHaraRudyParameters {
enum varType {
  VT_amplitude = vtFirst,
  VT_duration,
  VT_celltype,
  VT_Na_o,
  VT_Ca_o,
  VT_K_o,
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
  VT_Init_m,
  VT_Init_h_fast,
  VT_Init_h_slow,
  VT_Init_j,
  VT_Init_h_CaMK_slow,
  VT_Init_j_CaMK,
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
  VT_Init_n,
  VT_Init_f_CaMK_fast,
  VT_Init_f_Ca_CaMK_fast,
  VT_Init_x_r_fast,
  VT_Init_x_r_slow,
  VT_Init_x_s1,
  VT_Init_x_s2,
  VT_Init_x_K1,
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
  VT_A_h_fast,
  VT_A_h_slow,
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
  #ifdef HFM
  VT_HFM_Multiplier,
  #endif // ifdef HFM
  #ifdef HF
  VT_HF_Gomez_1,
  VT_HF_Gomez_2,
  VT_HF_Gomez_3,
  VT_HF_Gomez_4,
  VT_HF_Gomez_5,
  VT_HF_Gomez_6,
  VT_HF_Gomez_7,
  VT_HF_Gomez_8,
  VT_HF_Gomez_9,
  VT_HF_Gomez_0,
  #endif // ifdef HF
  vtLast
};
}  // namespace NS_OHaraRudyParameters

using namespace NS_OHaraRudyParameters;

class OHaraRudyParameters : public vbNewElphyParameters {
 public:
  OHaraRudyParameters(const char *, ML_CalcType);
  ~OHaraRudyParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);

  ML_CalcType m_inf[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType h_inf[RTDT];
  ML_CalcType exptau_h_fast[RTDT];
  ML_CalcType exptau_h_slow[RTDT];
  ML_CalcType j_inf[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType h_CaMK_inf[RTDT];
  ML_CalcType exptau_h_CaMK_slow[RTDT];
  ML_CalcType j_CaMK_inf[RTDT];
  ML_CalcType exptau_j_CaMK[RTDT];
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
  ML_CalcType R_Kr[RTDT];
  ML_CalcType x_s1_inf[RTDT];
  ML_CalcType exptau_x_s1[RTDT];
  ML_CalcType x_s2_inf[RTDT];
  ML_CalcType exptau_x_s2[RTDT];
  ML_CalcType x_K1_inf[RTDT];
  ML_CalcType exptau_x_K1[RTDT];
  ML_CalcType R_K1[RTDT];
  ML_CalcType x_Kb[RTDT];
};  // class OHaraRudyParameters

#endif  // ifndef OHARARUDYPARAMETERS_H
