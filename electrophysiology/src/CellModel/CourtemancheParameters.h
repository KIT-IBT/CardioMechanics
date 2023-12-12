/*
 * File: CourtemancheParameters.h
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


#ifndef CourtemancheParameters_H
#define CourtemancheParameters_H

/// Model variations
/// Stretch avtivated channel by Kuijpers et al. 2007
#define ISAC

/// Troponin dynamics from Land et al. tension model to use mechano-electric feedback.
#define TRPN

#include <ParameterLoader.h>

namespace NS_CourtemancheParameters {
enum varType {
  VT_R = vtFirst,
  VT_Tx,
  VT_F,
  VT_C_m,
  VT_Vcell,
  VT_Vcell_outside,
  VT_part_i,
  VT_part_up,
  VT_part_rel,
  VT_Tfac,
  VT_K_o,
  VT_Na_o,
  VT_Ca_o,
  VT_g_Na,
  VT_g_to,
  VT_g_Kr,
  VT_g_Ks,
  VT_g_CaL,
  VT_g_K1,
  VT_g_bNa,
  VT_g_bCa,
  VT_I_NaKmax,
  VT_k_mNai,
  VT_k_mKo,
  VT_I_pCamax,
  VT_k_mpCa,
  VT_k_NaCa,
  VT_k_mNa,
  VT_k_mCa,
  VT_k_sat,
  VT_gamm,
  VT_CMDNmax,
  VT_TRPNmax,
  VT_CSQNmax,
  VT_K_mcmdn,
  VT_K_mtrpn,
  VT_K_mcsqn,
  VT_k_rel,
  VT_k_up,
  VT_I_upmax,
  VT_Ca_upmax,
  VT_t_tr,
  VT_shiftm,
  VT_shifth,
  VT_shiftj,
  VT_shiftoa,
  VT_C_oa_1,
  VT_C_oa_2,
  VT_C_oa_3,
  VT_C_oa_4,
  VT_C_oa_5,
  VT_C_oa_6,
  VT_CaL_fT_1,
  VT_Init_Na_i,
  VT_Init_K_i,
  VT_Init_m,
  VT_Init_oa,
  VT_Init_ua,
  VT_Init_d,
  VT_Init_f_Ca,
  VT_Init_Ca_i,
  VT_Init_Ca_up,
  VT_Init_Ca_rel,
  VT_Init_h,
  VT_Init_j,
  VT_Init_oi,
  VT_Init_ui,
  VT_Init_Xr,
  VT_Init_Xs,
  VT_Init_f,
  VT_Init_u,
  VT_Init_v,
  VT_Init_w,
  VT_Init_Vm,
  VT_ua_a1,
  VT_ua_a2,
  VT_ua_a3,
  VT_ua_a4,
  VT_ua_a5,
  VT_ua_b1,
  VT_ua_b2,
  VT_ua_b3,
  VT_ua_m1,
  VT_ua_m2,
  VT_ua_KQ10,
  VT_ui_a1,
  VT_ui_a2,
  VT_ui_a3,
  VT_ui_a4,
  VT_ui_b1,
  VT_ui_b2,
  VT_ui_m1,
  VT_ui_m2,
  VT_ui_KQ10,
  VT_g_Kur1,
  VT_g_Kur2,
  VT_g_Kur3,
  VT_g_Kur4,
  VT_g_K1_1,
  VT_g_K1_2,
  VT_Xr_a1,
  VT_Xr_a2,
  VT_Xr_a3,
  VT_Xr_b1,
  VT_Xr_b2,
  VT_Xr_KQ10,
  VT_Xr_m1,
  VT_Xr_m2,
  VT_g_Kr1,
  VT_g_Kr2,
  VT_Xs_a1,
  VT_Xs_a2,
  VT_Xs_a3,
  VT_Xs_b1,
  VT_Xs_b2,
  VT_Xs_KQ10,
  VT_Xs_m1,
  VT_Xs_m2,
  VT_m_a1,
  VT_m_a2,
  VT_m_b1,
  VT_m_b2,
  VT_h_a1,
  VT_h_a2,
  VT_h_b1,
  VT_h_b2,
  VT_h_b3,
  VT_h_b4,
  VT_h_b5,
  VT_h_b6,
  VT_j_a1,
  VT_j_a2,
  VT_j_a3,
  VT_j_a4,
  VT_j_a5,
  VT_j_b1,
  VT_j_b2,
  VT_j_b3,
  VT_j_b4,
  VT_j_b5,
  VT_j_b6,
  VT_tau_j_mult,
  VT_tau_f_mult,
  VT_E_rL,
  VT_RTdF,
  VT_FdRT,
  VT_RTd2F,
  VT_Vi,
  VT_Vup,
  VT_Vrel,
  VT_VupdVi,
  VT_VreldVi,
  VT_VreldVup,
  VT_CmdFvi,
  VT_Cmd2F,
  VT_csqnkm,
  VT_kmcsqnm2,
  VT_kkmcsqn,
  VT_cmdnkm,
  VT_kmcmdnm2,
  VT_kkmcmdn,
  VT_trpnkm,
  VT_kmtrpnm2,
  VT_kkmtrpn,
  VT_kupleak,
  VT_dt_tr,
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
}  // namespace NS_CourtemancheParameters

using namespace NS_CourtemancheParameters;

class CourtemancheParameters : public vbNewElphyParameters {
 public:
  CourtemancheParameters(const char *, ML_CalcType);
  ~CourtemancheParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);

  ML_CalcType exp250;
  ML_CalcType exp50;


  ML_CalcType m_m[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType m_h[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType m_j[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType m_oa[RTDT];
  ML_CalcType exptau_oa[RTDT];
  ML_CalcType m_oi[RTDT];
  ML_CalcType exptau_oi[RTDT];
  ML_CalcType m_ua[RTDT];
  ML_CalcType exptau_ua[RTDT];
  ML_CalcType m_ui[RTDT];
  ML_CalcType exptau_ui[RTDT];
  ML_CalcType m_Xr[RTDT];
  ML_CalcType exptau_Xr[RTDT];
  ML_CalcType m_Xs[RTDT];
  ML_CalcType exptau_Xs[RTDT];
  ML_CalcType m_d[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType m_f[RTDT];
  ML_CalcType exptau_f[RTDT];
  ML_CalcType m_w[RTDT];
  ML_CalcType exptau_w[RTDT];
  ML_CalcType expVm[RTDT];
  ML_CalcType exp2gamm[RTDT];
  ML_CalcType CK1[RTDT];
  ML_CalcType g_Kur[RTDT];
  ML_CalcType CKr[RTDT];
  ML_CalcType CCaL[RTDT];
  ML_CalcType CNaK[RTDT];
  ML_CalcType CNaCa[RTDT];
};  // class CourtemancheParameters
#endif  // ifndef CourtemancheParameters_H
