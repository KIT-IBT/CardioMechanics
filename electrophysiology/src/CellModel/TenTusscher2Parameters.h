/*
 * File: TenTusscher2Parameters.h
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


#ifndef TENTUSSCHER2PARAMETERS_H
#define TENTUSSCHER2PARAMETERS_H
//#define ACTIVATE_ISAC_CHANNEL
//#define SAC_KUIJPERS
//#define SAC_SACHS

// uncomment if Markov I_Na by Clancy and Rudy should be used
// #define MARKOV_I_NA

#include <ParameterLoader.h>
#include <TenTusscher2IschemiaSetup.h>

namespace NS_TenTusscher2Parameters {
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
  VT_f2_init,
  VT_fCa_init,
  VT_Rq_init,
  VT_O_init,
  VT_K_o,
  VT_Ca_o,
  VT_Na_o,
  VT_Vc,
  VT_Vsr,
  VT_Vss,
  VT_Vrel,
  VT_Vcell,
  VT_ks1,
  VT_ks2,
  VT_k3,
  VT_k4,
  VT_EC,
  VT_max_sr,
  VT_min_sr,
  VT_Vxfer,
  VT_Bufc,
  VT_Kbufc,
  VT_Bufsr,
  VT_Kbufsr,
  VT_Bufss,
  VT_Kbufss,
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
  VT_CaSS_init,
  VT_Nai_init,
  VT_Ki_init,
  VT_lenght,
  VT_radius,
  VT_pi,
  VT_CaiMin,
  VT_CaiMax,
  VT_s_inf_vHalf,
  VT_tau_s_f1,
  VT_tau_s_slope1,
  VT_tau_s_vHalf1,
  VT_tau_s_enable,
  VT_vol,
  VT_vi,
  VT_inverseviF,
  VT_inverseviF2,
  VT_inversevssF2,
  VT_volforCall,
  VT_INVERSECAPACITANCE,
  VT_VcdVsr,
  VT_Kupsquare,
  VT_BufcPKbufc,
  VT_Kbufcsquare,
  VT_Kbufc2,
  VT_BufsrPKbufsr,
  VT_Kbufsrsquare,
  VT_Kbufsr2,
  VT_KmNai3,
  VT_Nao3,
  VT_StepCai,
  VT_inverseRTONF,
  VT_RTONF,
  VT_tau_s_add,
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
  VT_m_Xs_1,
  VT_m_Xs_2,
  VT_a_Xs_1,
  VT_a_Xs_2,
  VT_b_Xs_1,
  VT_b_Xs_2,
  VT_tau_Xs_add,
  VT_K_Q10Xs,
  VT_stim_duration,
  VT_Amp,
  VT_alpha3_div,
  VT_alpha5_div,
  VT_beta5_div,
  VT_alpha6,
  VT_beta6,
  VT_initMINALC3,
  VT_initMINALC2,
  VT_initMINALC1,
  VT_initMINALO,
  VT_initMINAUC3,
  VT_initMINAUC2,
  VT_initMINAUC1,
  VT_initMINAUO,
  VT_initMINAUIC3,
  VT_initMINAUIC2,
  VT_initMINAUIF,
  VT_initMINAUIM1,
  VT_initMINAUIM2,
#ifdef ACTIVATE_ISAC_CHANNEL
#ifdef SAC_SACHS
    VT_g_SAC,
    VT_alphaSAC,
    VT_ESAC,
    VT_KSAC,
#endif //ifdef SAC_SACHS
#ifdef SAC_KUIJPERS
    VT_PsacNa,
    VT_PsacK,
    VT_PsacCa,
    VT_Gsac,
    VT_Ksac,
    VT_alphasac,
    VT_lambdasac,
    VT_gsac,
#endif  // ifdef SAC_KUIJPERS
#endif  // ifdef ACTIVATE_ISAC_CHANNEL
#ifdef ACTIVATE_IKATP_CHANNEL
  VT_f_T,
  VT_nicholsarea,
  VT_gkatp,
  VT_gammaconst,
  VT_Mgi,
  VT_atpi,
  VT_adpi,
  VT_Km_factor,
# ifdef ISCHEMIA
  VT_IschemiaStart,
  VT_IschemiaStage1,
  VT_IschemiaStage2,
  VT_ZoneFactor,
  VT_DiffusionFactor,
  VT_RestoreIschemia,
  VT_Ko_ZoneFactor_Begin,
  VT_Ko_ZoneFactor_End,
  VT_fpH_ZoneFactor_Begin,
  VT_fpH_ZoneFactor_End,
  VT_pO_ZoneFactor_Begin,
  VT_pO_ZoneFactor_End,
  VT_K_o_stage1,
  VT_K_o_stage2,
  VT_gCaL_stage1,
  VT_gCaL_stage2,
  VT_gNa_stage1,
  VT_gNa_stage2,
  VT_Mgi_stage1,
  VT_Mgi_stage2,
  VT_ATP_stage1,
  VT_ATP_stage2,
  VT_ADP_stage1,
  VT_ADP_stage2,
  VT_dVmNa_stage0,
  VT_dVmNa_stage1,
  VT_dVmNa_stage2,
  VT_knak_1b,
  VT_kNaCa_1b,
  VT_Vmaxup_1b,
  VT_Vrel_1b,
# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
    
  vtLast
};
}  // namespace NS_TenTusscher2Parameters

using namespace NS_TenTusscher2Parameters;

class TenTusscher2Parameters : public vbNewElphyParameters {
 public:
  TenTusscher2Parameters(const char *, ML_CalcType);
  ~TenTusscher2Parameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  void InitTableWithtinc(ML_CalcType);
  ML_CalcType rec_ipK[RTDT];
  ML_CalcType d_inf[RTDT];
  ML_CalcType f_inf[RTDT];
  ML_CalcType f2_inf[RTDT];
  ML_CalcType m_inf[RTDT];
  ML_CalcType h_inf[RTDT];
  ML_CalcType j_inf[RTDT];
  ML_CalcType Xr1_inf[RTDT];
  ML_CalcType Xr2_inf[RTDT];
  ML_CalcType Xs_inf[RTDT];
  ML_CalcType r_inf[RTDT];
  ML_CalcType s_inf[RTDT];
  ML_CalcType NaCa_P1[RTDT];
  ML_CalcType NaCa_P2[RTDT];
  ML_CalcType NaK_P1[RTDT];
  ML_CalcType CaL_P1[RTDT];
  ML_CalcType CaL_P2[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType exptau_Xr1[RTDT];
  ML_CalcType exptau_Xr2[RTDT];
  ML_CalcType exptau_Xs[RTDT];
  ML_CalcType exptau_s[RTDT];
  ML_CalcType exptau_r[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType exptau_f[RTDT];
  ML_CalcType exptau_f2[RTDT];
  ML_CalcType ECA[RTDT];
#ifdef ACTIVATE_IKATP_CHANNEL
  ML_CalcType KhNa[RTDT];
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
#ifdef MARKOV_I_NA
  ML_CalcType a11[RTDT];
  ML_CalcType a12[RTDT];
  ML_CalcType a13[RTDT];
  ML_CalcType b11[RTDT];
  ML_CalcType b12[RTDT];
  ML_CalcType b13[RTDT];
  ML_CalcType a3[RTDT];
  ML_CalcType b3[RTDT];
  ML_CalcType a2[RTDT];
  ML_CalcType b2[RTDT];
  ML_CalcType a4[RTDT];
  ML_CalcType b4[RTDT];
  ML_CalcType a5[RTDT];
  ML_CalcType b5[RTDT];
#endif  // ifdef MARKOV_I_NA
#ifdef SAC_KUIJPERS
    ML_CalcType Csac1[RTDT];
    ML_CalcType Csac2[RTDT];
#endif  // ifdef SAC_KUIJPERS
};  // class TenTusscher2Parameters
#endif  // ifndef TENTUSSCHER2PARAMETERS_H
