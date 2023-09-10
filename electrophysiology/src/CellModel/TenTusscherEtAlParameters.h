/*      File: TenTusscherEtAlParameters.h
    automatically created by ExtractParameterClass.pl - done by dw (20.04.2007)
    Institute of Biomedical Engineering, Universitt Karlsruhe (TH)
    send comments to dw@ibt.uka.de      */

#ifndef TENTUSSCHERETALPARAMETERS_H
#define TENTUSSCHERETALPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_TenTusscherEtAlParameters {
enum varType {
  VT_R = vtFirst,
  VT_Tx,
  VT_F,
  VT_m_init,
  VT_h_init,
  VT_j_init,
  VT_xr1_init,
  VT_xr2_init,
  VT_xs_init,
  VT_r_init,
  VT_s_init,
  VT_d_init,
  VT_f_init,
  VT_fCa_init,
  VT_g_init,
  VT_K_o,
  VT_Ca_o,
  VT_Na_o,
  VT_Vc,
  VT_Vsr,
  VT_Bufc,
  VT_Kbufc,
  VT_Bufsr,
  VT_Kbufsr,
  VT_taufca,
  VT_taug,
  VT_Vmaxup,
  VT_Kup,
  VT_C,
  VT_g_Kr,
  VT_pKNa,
  VT_g_Ks,
  VT_g_K1,
  VT_g_to,
  VT_g_Na,
  VT_g_bNa,
  VT_KmK,
  VT_KmNa,
  VT_knak,
  VT_g_CaL,
  VT_g_bCa,
  VT_kNaCa,
  VT_KmNai,
  VT_KmCa,
  VT_ksat,
  VT_n,
  VT_g_pCa,
  VT_KpCa,
  VT_g_pK,
  VT_V_init,
  VT_Cai_init,
  VT_CaSR_init,
  VT_Nai_init,
  VT_Ki_init,
  VT_s_inf_vHalf,
  VT_tau_s_f1,
  VT_tau_s_slope1,
  VT_tau_s_vHalf1,
  VT_tau_s_f2,
  VT_tau_s_f3,
  VT_m_Xr1_1,
  VT_m_Xr1_2,
  VT_a_Xr1_1,
  VT_a_Xr1_2,
  VT_b_Xr1_1,
  VT_b_Xr1_2,
  VT_K_Q10Xr1,
  VT_m_Xr2_1,
  VT_m_Xr2_2,
  VT_a_Xr2_1,
  VT_a_Xr2_2,
  VT_b_Xr2_1,
  VT_b_Xr2_2,
  VT_K_Q10Xr2,
  VT_inverseVcF2,
  VT_inverseVcF2C,
  VT_inverseVcFC,
  VT_VcdVsr,
  VT_Kupsquare,
  VT_BufcPKbufc,
  VT_Kbufcsquare,
  VT_Kbufc2,
  VT_BufsrPKbufsr,
  VT_Kbufsrsquare,
  VT_Kbufsr2,
  VT_KopKNaNao,
  VT_KmNai3,
  VT_Nao3,
  VT_inverseRTONF,
  VT_RTONF,
  VT_Amp,
  vtLast
};
}  // namespace NS_TenTusscherEtAlParameters

using namespace NS_TenTusscherEtAlParameters;

class TenTusscherEtAlParameters : public vbNewElphyParameters {
 public:
  TenTusscherEtAlParameters(const char *, ML_CalcType);
  ~TenTusscherEtAlParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  ML_CalcType rec_ipK[RTDT];
  ML_CalcType d_inf[RTDT];
  ML_CalcType f_inf[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType exptau_Xr1[RTDT];
  ML_CalcType exptau_Xr2[RTDT];
  ML_CalcType exptau_Xs[RTDT];
  ML_CalcType exptau_r[RTDT];
  ML_CalcType exptau_s[RTDT];
  ML_CalcType m_inf[RTDT];
  ML_CalcType h_inf[RTDT];
  ML_CalcType j_inf[RTDT];
  ML_CalcType Xr1_inf[RTDT];
  ML_CalcType Xr2_inf[RTDT];
  ML_CalcType Xs_inf[RTDT];
  ML_CalcType r_inf[RTDT];
  ML_CalcType s_inf[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType exptau_f[RTDT];
  ML_CalcType NaCa_P1[RTDT];
  ML_CalcType NaCa_P2[RTDT];
  ML_CalcType NaK_P1[RTDT];
  ML_CalcType CaL_P1[RTDT];
  ML_CalcType CaL_P2[RTDT];
};  // class TenTusscherEtAlParameters
#endif  // ifndef TENTUSSCHERETALPARAMETERS_H
