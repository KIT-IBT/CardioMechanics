/* File: KoivumaekiEtAlParameters.h
        Institute of Biomedical Engineering, Karlsruhe Institute of Technology (KIT) */

#ifndef KOIVUMAEKIETALPARAMETERS_H
#define KOIVUMAEKIETALPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_KoivumaekiEtAlParameters {
enum varType {
  VT_R = vtFirst,
  VT_T,
  VT_F,
  VT_Cm,
  VT_stim_duration,
  VT_Amp,
  VT_Na_o,
  VT_Ca_o,
  VT_K_o,
  VT_V_ss,
  VT_V_corrcell,
  VT_r_junct,
  VT_l_cell,
  VT_cell_dillation,
  VT_dx,
  VT_totCMDN,
  VT_totSLlow,
  VT_totSLhigh,
  VT_totFura2,
  VT_totHTRPN,
  VT_totLTRPN,
  VT_KdCMDN,
  VT_KdSLlow,
  VT_KdSLhigh,
  VT_KdFura2,
  VT_kHTRPNon,
  VT_kHTRPNoff,
  VT_kLTRPNon,
  VT_kLTRPNoff,
  VT_totCSQN,
  VT_KdCSQN,
  VT_totBNa,
  VT_KdBNa,
  VT_PNa,
  VT_ECa_app,
  VT_gCaL,
  VT_kCan,
  VT_kCa,
  VT_gt,
  VT_gsus,
  VT_gKs,
  VT_gKr,
  VT_gK1,
  VT_gNab,
  VT_gCab,
  VT_INaKmax,
  VT_kNaKK,
  VT_kNaKNa,
  VT_ICaPmax,
  VT_kCaP,
  VT_kNaCa,
  VT_gam,
  VT_dNaCa,
  VT_gIf,
  VT_DCa,
  VT_DCaSR,
  VT_DCaBm,
  VT_DNa,
  VT_SERCAKmf,
  VT_SERCAKmr,
  VT_k4,
  VT_cpumps,
  VT_kSRleak,
  VT_tau_fca,
  VT_RyRtauadapt,
  VT_RyRtauactss,
  VT_RyRtauinactss,
  VT_RyRtauact,
  VT_RyRtauinact,
  VT_V_init,
  VT_m_init,
  VT_h1_init,
  VT_h2_init,
  VT_d_init,
  VT_f1_init,
  VT_fca_init,
  VT_r_init,
  VT_s_init,
  VT_sus_r_init,
  VT_sus_s_init,
  VT_n_init,
  VT_pa_init,
  VT_y_init,
  VT_RyR_oss_init,
  VT_RyR_css_init,
  VT_RyR_ass_init,
  VT_RyR_o1_init,
  VT_RyR_c1_init,
  VT_RyR_a1_init,
  VT_RyR_o2_init,
  VT_RyR_c2_init,
  VT_RyR_a2_init,
  VT_RyR_o3_init,
  VT_RyR_c3_init,
  VT_RyR_a3_init,
  VT_SERCA_Ca1_init,
  VT_SERCA_Ca2_init,
  VT_SERCA_Ca3_init,
  VT_SERCA_Cass_init,
  VT_Na_ss_init,
  VT_Na_i_init,
  VT_K_i_init,
  VT_HTRPNCa1_init,
  VT_HTRPNCa2_init,
  VT_HTRPNCa3_init,
  VT_HTRPNCa4_init,
  VT_LTRPNCa1_init,
  VT_LTRPNCa2_init,
  VT_LTRPNCa3_init,
  VT_LTRPNCa4_init,
  VT_Ca_ss_init,
  VT_Ca_i1_init,
  VT_Ca_i2_init,
  VT_Ca_i3_init,
  VT_Ca_i4_init,
  VT_Ca_SR1_init,
  VT_Ca_SR2_init,
  VT_Ca_SR3_init,
  VT_Ca_SR4_init,
  VT_RTdF,
  VT_FdRT,
  VT_RTd2F,
  VT_Dlcell,
  VT_Ddcell,
  VT_Dvcell,
  VT_rstart,
  VT_rend,
  VT_Aj_nj,
  VT_xj_nj,
  VT_xj_nj_Nai,
  VT_V_nonjunct1,
  VT_V_nonjunct2,
  VT_V_nonjunct3,
  VT_V_nonjunct4,
  VT_V_nonjunct_Nai,
  VT_V_cytosol,
  VT_V_SR1,
  VT_V_SR2,
  VT_V_SR3,
  VT_V_SR4,
  VT_exp_tau_fca,
  VT_pow_kNaKNa,
  VT_k1,
  VT_k2,
  VT_k3,
  VT_powKo,
  VT_Na_o3,
  VT_CJj_nj,
  VT_CJNa,
  vtLast
};
}  // namespace NS_KoivumaekiEtAlParameters

using namespace NS_KoivumaekiEtAlParameters;

class KoivumaekiEtAlParameters : public vbNewElphyParameters {
 public:
  KoivumaekiEtAlParameters(const char *, ML_CalcType);
  ~KoivumaekiEtAlParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  ML_CalcType exptau_r[RTDT];
  ML_CalcType r_infinity[RTDT];
  ML_CalcType sus_r_infinity[RTDT];
  ML_CalcType exptau_sus_r[RTDT];
  ML_CalcType sus_s_infinity[RTDT];
  ML_CalcType exptau_sus_s[RTDT];
  ML_CalcType m_infinity[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType h_infinity[RTDT];
  ML_CalcType exptau_h1[RTDT];
  ML_CalcType exptau_h2[RTDT];
  ML_CalcType d_infinity[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType f_infinity[RTDT];
  ML_CalcType exptau_f1[RTDT];
  ML_CalcType exptau_s[RTDT];
  ML_CalcType s_infinity[RTDT];
  ML_CalcType exptau_n[RTDT];
  ML_CalcType n_infinity[RTDT];
  ML_CalcType exptau_pa[RTDT];
  ML_CalcType p_a_infinity[RTDT];
  ML_CalcType pip[RTDT];
  ML_CalcType exptau_y[RTDT];
  ML_CalcType y_infinity[RTDT];
  ML_CalcType CexpINa[RTDT];
  ML_CalcType CI_NaCa1[RTDT];
  ML_CalcType CI_NaCa2[RTDT];
};  // class KoivumaekiEtAlParameters

#endif  // ifndef KOIVUMAEKIETALPARAMETERS_H
