/* File: MaleckarEtAlParameters.h
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#ifndef MALECKARETALPARAMETERS_H
#define MALECKARETALPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_MaleckarEtAlParameters {
enum varType {
  VT_R = vtFirst,
  VT_T,
  VT_F,
  VT_Cm,
  VT_stim_duration,
  VT_Amp,
  VT_P_Na,
  VT_g_Ca_L,
  VT_E_Ca_app,
  VT_k_Ca,
  VT_g_t,
  VT_g_kur,
  VT_g_K1,
  VT_g_Ks,
  VT_g_Kr,
  VT_g_B_Na,
  VT_g_B_Ca,
  VT_K_NaK_K,
  VT_i_NaK_max,
  VT_pow_K_NaK_Na_15,
  VT_i_CaP_max,
  VT_k_CaP,
  VT_K_NaCa,
  VT_d_NaCa,
  VT_gamma_Na,
  VT_ACh,
  VT_phi_Na_en,
  VT_Vol_i,
  VT_V_corrcell,
  VT_Vol_d,
  VT_tau_di,
  VT_Mg_i,
  VT_Vol_c,
  VT_tau_Na,
  VT_tau_K,
  VT_tau_Ca,
  VT_Na_b,
  VT_Ca_b,
  VT_K_b,
  VT_I_up_max,
  VT_k_cyca,
  VT_k_srca,
  VT_k_xcs,
  VT_alpha_rel,
  VT_Vol_up,
  VT_Vol_rel,
  VT_r_recov,
  VT_tau_tr,
  VT_k_rel_i,
  VT_k_rel_d,
  VT_V_init,
  VT_Na_c_init,
  VT_Na_i_init,
  VT_m_init,
  VT_h1_init,
  VT_h2_init,
  VT_Ca_d_init,
  VT_d_L_init,
  VT_f_L1_init,
  VT_f_L2_init,
  VT_K_c_init,
  VT_K_i_init,
  VT_r_init,
  VT_s_init,
  VT_a_ur_init,
  VT_i_ur_init,
  VT_n_init,
  VT_pa_init,
  VT_Ca_c_init,
  VT_Ca_i_init,
  VT_O_C_init,
  VT_O_TC_init,
  VT_O_TMgC_init,
  VT_O_TMgMg_init,
  VT_O_init,
  VT_Ca_rel_init,
  VT_Ca_up_init,
  VT_O_Calse_init,
  VT_F1_init,
  VT_F2_init,
  VT_dT_yKur1,
  VT_dT_yKur2,
  VT_Init_yKur,
  VT_InitTableDone,
  VT_tincforTab,
  VT_RTdF,
  VT_FdRT,
  VT_FFdRT,
  VT_RTd2F,
  VT_Ci_di,
  VT_Ci_tr,
  VT_dVol_iF,
  VT_dVol_cF,
  VT_dVol_dF,
  VT_dVol_upF,
  VT_dVol_relF,
  VT_k_xcs2ds,
  vtLast
};
}  // namespace NS_MaleckarEtAlParameters

using namespace NS_MaleckarEtAlParameters;

class MaleckarEtAlParameters : public vbNewElphyParameters {
 public:
  MaleckarEtAlParameters(const char *, ML_CalcType);
  ~MaleckarEtAlParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  ML_CalcType exptau_r[RTDT];
  ML_CalcType r_infinity[RTDT];
  ML_CalcType a_ur_infinity[RTDT];
  ML_CalcType exptau_a_ur[RTDT];
  ML_CalcType i_ur_infinity[RTDT];
  ML_CalcType exptau_i_ur[RTDT];
  ML_CalcType m_infinity[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType h_infinity[RTDT];
  ML_CalcType exptau_h1[RTDT];
  ML_CalcType exptau_h2[RTDT];
  ML_CalcType d_L_infinity[RTDT];
  ML_CalcType exptau_d_L[RTDT];
  ML_CalcType f_L_infinity[RTDT];
  ML_CalcType exptau_f_L1[RTDT];
  ML_CalcType exptau_f_L2[RTDT];
  ML_CalcType exptau_s[RTDT];
  ML_CalcType s_infinity[RTDT];
  ML_CalcType exptau_n[RTDT];
  ML_CalcType n_infinity[RTDT];
  ML_CalcType exptau_pa[RTDT];
  ML_CalcType p_a_infinity[RTDT];
  ML_CalcType pip[RTDT];
  ML_CalcType Q_tot[RTDT];
  ML_CalcType KIKACH[RTDT];
  ML_CalcType CexpINa[RTDT];
  ML_CalcType CI_NaCa1[RTDT];
  ML_CalcType CI_NaCa2[RTDT];
  ML_CalcType exptau_yKur[RTDT];
  ML_CalcType yKur_infinity[RTDT];
};  // class MaleckarEtAlParameters

#endif  // ifndef MALECKARETALPARAMETERS_H
