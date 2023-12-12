/*
 * File: FabbriParameters.h
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


#ifndef FabbriParameters_H
#define FabbriParameters_H


#include <ParameterLoader.h>

namespace NS_FabbriParameters {
enum varType {
  VT_Vari_Nai = vtFirst,
  VT_Vari_Ki,
  VT_R,
  VT_T,
  VT_F,
  VT_C_m,
  VT_ACh,
  VT_Iso_1_uM,
  VT_Nao,
  VT_Init_Ki,
  VT_Ko,
  VT_Cao,
  VT_g_fNa,
  VT_g_fK,
  VT_Km_Kp,
  VT_Km_Nap,
  VT_i_NaK_max,
  VT_K_NaCa,
  VT_Qci,
  VT_Qn,
  VT_Qco,
  VT_K3ni,
  VT_Kci,
  VT_K1ni,
  VT_K2ni,
  VT_Kcni,
  VT_K3no,
  VT_K1no,
  VT_K2no,
  VT_Kco,
  VT_g_na,
  VT_P_CaL,
  VT_k_dl,
  VT_V_dl,
  VT_k_fl,
  VT_V_fl,
  VT_alpha_fCa,
  VT_Km_fCa,
  VT_P_CaT,
  VT_ks,
  VT_MaxSR,
  VT_MinSR,
  VT_EC50_SR,
  VT_EC50_SK,
  VT_n_SK,
  VT_HSR,
  VT_koCa,
  VT_kiCa,
  VT_kim,
  VT_kom,
  VT_tau_dif_Ca,
  VT_tau_tr,
  VT_K_up,
  VT_slope_up,
  VT_TC_tot,
  VT_TMC_tot,
  VT_CM_tot,
  VT_CQ_tot,
  VT_kf_TC,
  VT_kf_TMM,
  VT_kf_TMC,
  VT_kf_CM,
  VT_kf_CQ,
  VT_kb_TC,
  VT_kb_TMC,
  VT_kb_TMM,
  VT_kb_CM,
  VT_kb_CQ,
  VT_Init_Mgi,
  VT_V_jsr_part,
  VT_V_i_part,
  VT_V_nsr_part,
  VT_R_cell,
  VT_L_cell,
  VT_L_sub,
  VT_g_Kur,
  VT_g_to,
  VT_g_Kr,
  VT_g_Ks,
  VT_g_KACh,
  VT_E_K,
  VT_k34,
  VT_alpha_a,
  VT_RTd2F,
  VT_Init_Vm,
  VT_Init_Ca_sub,
  VT_Init_Nai,
  VT_Init_y,
  VT_Init_m,
  VT_Init_h,
  VT_Init_dL,
  VT_Init_fL,
  VT_Init_fCa,
  VT_Init_dT,
  VT_Init_fT,
  VT_Init_R_Ca_rel,
  VT_Init_O_Ca_rel,
  VT_Init_I_Ca_rel,
  VT_Init_RI_Ca_rel,
  VT_Init_Ca_jsr,
  VT_Init_Ca_nsr,
  VT_Init_Cai,
  VT_Init_fTMM,
  VT_Init_fCMi,
  VT_Init_fCMs,
  VT_Init_fTC,
  VT_Init_fTMC,
  VT_Init_fCQ,
  VT_Init_r_Kur,
  VT_Init_s,
  VT_Init_q,
  VT_Init_r_to,
  VT_Init_paS,
  VT_Init_paF,
  VT_Init_piy,
  VT_Init_n,
  VT_Init_a,
  VT_Init_x,
  VT_g_na_l,
  VT_g_na_mut,
  VT_g_SK,
  VT_V_cell,
  VT_V_sub,
  VT_V_jsr,
  VT_V_i,
  VT_V_nsr,
  VT_Iso_shift_If,
  VT_Iso_INak_increase,
  VT_Iso_ICaL_increase,
  VT_Iso_shift_dL,
  VT_Iso_slope_dL,
  VT_Iso_IKs_increase,
  VT_Iso_shift_n,
  VT_Iso_b_up,
  VT_K_Iso_shift,
  VT_K_Iso_increase,
  VT_K_Iso_slope_dL,
  VT_ACh_if_shift,
  VT_ACh_iCaL_block,
  VT_y_inf_a,
  VT_y_inf_b,
  VT_V_0p5_y_inf,
  VT_tau_y_a_shift,
  VT_tau_y_b_shift,
  VT_V_0p5_m_inf,
  VT_k_m,
  VT_tau_m_shift,
  VT_V_0p5_h_inf,
  VT_k_h,
  VT_tau_h_gain,
  VT_V_0p5_n_inf,
  VT_k_n,
  VT_R231C,
  vtLast
};
} // namespace NS_FabbriParameters

using namespace NS_FabbriParameters;

class FabbriParameters : public vbNewElphyParameters {
 public:
  FabbriParameters(const char *, ML_CalcType);
  ~FabbriParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);

  ML_CalcType m_y[RTDT];
  ML_CalcType exptau_y[RTDT];
  ML_CalcType m_m[RTDT];
  ML_CalcType m_m_mut[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType exptau_m_mut[RTDT];
  ML_CalcType m_h[RTDT];
  ML_CalcType m_h_mut[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType exptau_h_mut[RTDT];
  ML_CalcType m_dL[RTDT];
  ML_CalcType exptau_dL[RTDT];
  ML_CalcType m_fL[RTDT];
  ML_CalcType exptau_fL[RTDT];
  ML_CalcType m_dT[RTDT];
  ML_CalcType exptau_dT[RTDT];
  ML_CalcType m_fT[RTDT];
  ML_CalcType exptau_fT[RTDT];
  ML_CalcType m_r_Kur[RTDT];
  ML_CalcType exptau_r_Kur[RTDT];
  ML_CalcType m_s[RTDT];
  ML_CalcType exptau_s[RTDT];
  ML_CalcType m_q[RTDT];
  ML_CalcType exptau_q[RTDT];
  ML_CalcType m_r_to[RTDT];
  ML_CalcType exptau_r_to[RTDT];
  ML_CalcType m_pa[RTDT];
  ML_CalcType exptau_paF[RTDT];
  ML_CalcType exptau_paS[RTDT];
  ML_CalcType m_pi[RTDT];
  ML_CalcType exptau_pi[RTDT];
  ML_CalcType m_n[RTDT];
  ML_CalcType exptau_n[RTDT];
  ML_CalcType m_a[RTDT];
  ML_CalcType exptau_a[RTDT];
  ML_CalcType ICs_on_Icontrol[RTDT];
  ML_CalcType doo[RTDT];
  ML_CalcType k21[RTDT];
  ML_CalcType k41[RTDT];
  ML_CalcType k23[RTDT];
  ML_CalcType k32[RTDT];
  ML_CalcType dNai[RTDT];
  ML_CalcType dKi[RTDT];
};  // class FabbriParameters
#endif  // ifndef FabbriParameters_H
