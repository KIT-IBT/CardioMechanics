/*
 * File: TWorldParameters.h
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


#ifndef TWORLDPARAMETERS_H
#define TWORLDPARAMETERS_H

/// TWorld model variations.

#include <ParameterLoader.h>

namespace NS_TWorldParameters {
enum varType {
  /// Simulation Parameters
  VT_amplitude = vtFirst,
  VT_duration,
  VT_celltype,
  VT_Na_o,
  VT_Ca_o,
  VT_K_o,
  VT_Cl_o,
    
  /// Multipliers
  VT_ICaL_fractionSS,
  VT_INaCa_fractionSS,
  VT_INaF_Multiplier,
  VT_INaL_Multiplier,
  VT_Ito_Multiplier,
  VT_Ito_slow_Multiplier,
  VT_Ito_fast_Multiplier,
  VT_ICaL_Multiplier,
  VT_IKr_Multiplier,
  VT_IKs_Multiplier,
  VT_IK1_Multiplier,
  VT_IKb_Multiplier,
  VT_INaCa_Multiplier,
  VT_INaK_Multiplier,
  VT_INab_Multiplier,
  VT_ICab_Multiplier,
  VT_IpCa_Multiplier,
  VT_ICaCl_Multiplier,
  VT_IClb_Multiplier,
  VT_Jrel_Multiplier,
  VT_Jup_Multiplier,
  VT_CMDN_Multiplier,
    
  /// PKA fractions
  VT_fINa_PKA,
  VT_fICaL_PKA,
  VT_fINaK_PKA,
  VT_fIKs_PKA,
  VT_fIfPLB_PKA,
  VT_fTnI_PKA,
  VT_fMyBPC_PKA,
  VT_PP1_tot,

  /// State Variables
  // Transmenbrane Voltage and Concentrations
  VT_Init_Vm,
  VT_Init_Na_dyad,
  VT_Init_Na_sl,
  VT_Init_Na_myo,
  VT_Init_K_myo,
  VT_Init_Cl_myo,
  VT_Init_Ca_dyad,
  VT_Init_Ca_sl,
  VT_Init_Ca_i,
  VT_Init_Ca_SR,
  // CaMK and Ca Signalling
  VT_Init_CaMK_trap,
  VT_Init_CaMK_f_ICaL,
  VT_Init_CaMK_f_RyR,
  VT_Init_CaMK_f_PLB,
  VT_Init_casig_serca_trap,
  // Sodium current (INa, INaL)
  VT_Init_m,
  VT_Init_h,
  VT_Init_j,
  VT_Init_h_p,
  VT_Init_j_p,
  VT_Init_m_PKA,
  VT_Init_h_PKA,
  VT_Init_j_PKA,
  VT_Init_h_both,
  VT_Init_j_both,
  VT_Init_m_L,
  VT_Init_h_L,
  VT_Init_h_L_p,
  // L-type calcium current (I_CaL, I_CaNa, I_CaK)
  VT_Init_d,
  VT_Init_f_fast,
  VT_Init_f_slow,
  VT_Init_f_Ca_fast,
  VT_Init_f_Ca_slow,
  VT_Init_j_Ca,
  VT_Init_f_p_fast,
  VT_Init_f_Ca_p_fast,
  VT_Init_d_PKA,
  VT_Init_f_PKA_fast,
  VT_Init_f_PKA_slow,
  VT_Init_f_Ca_PKA_fast,
  VT_Init_f_Ca_PKA_slow,
  VT_Init_f_both_fast,
  VT_Init_f_Ca_both_fast,
  VT_Init_n_Ca_dyad,
  VT_Init_n_Ca_sl,
  VT_Init_I_CaL_pureCDI_dyad,
  VT_Init_I_CaL_pureCDI_sl,
  // Transient outward current (Ito)
  VT_Init_a_slow,
  VT_Init_a_fast,
  VT_Init_i_slow,
  VT_Init_i_fast,
  VT_Init_a_p_slow,
  VT_Init_a_p_fast,
  VT_Init_i_p_slow,
  VT_Init_i_p_fast,
  // Rapid delayed rectifier current (IKr)
  VT_Init_C_0,
  VT_Init_C_1,
  VT_Init_C_2,
  VT_Init_O,
  VT_Init_I,
  // Slow delayed rectifier current (IKs)
  VT_Init_xs_dyad,
  VT_Init_xs_sl,
  // Calcium release from SR (Jrel, Jleak)
  VT_Init_J_rel_ICaLdep_act,
  VT_Init_J_rel_ICaLdep_f1,
  VT_Init_J_rel_ICaLdep_f2,
  VT_Init_ryr_R,
  VT_Init_ryr_O,
  VT_Init_ryr_I,
  VT_Init_ryr_CaRI,
  VT_Init_ryr_R_p,
  VT_Init_ryr_O_p,
  VT_Init_ryr_I_p,
  VT_Init_ryr_CaRI_p,
  // Buffering
  VT_Init_Buffer_NaBj,
  VT_Init_Buffer_NaBsl,
  VT_Init_Buffer_TnClow,
  VT_Init_Buffer_TnCHc,
  VT_Init_Buffer_TnCHm,
  VT_Init_Buffer_CaM,
  VT_Init_Buffer_Myosin_ca,
  VT_Init_Buffer_Myosin_mg,
  VT_Init_Buffer_SRB,
  VT_Init_Buffer_SLLj,
  VT_Init_Buffer_SLLsl,
  VT_Init_Buffer_SLHj,
  VT_Init_Buffer_SLHsl,
  VT_Init_Buffer_Csqn,
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Land-Niederer model of contraction
  ////////////////////////////////////////////////////////////////////////////////////////
  VT_perm50,
  VT_TRPN_n,
  VT_koff,
  VT_dr,
  VT_wfrac,
  VT_A_eff,
  VT_ktm_unblock,
  VT_beta_1,
  VT_beta_0,
  VT_gamma,
  VT_gamma_wu,
  VT_phi,
  VT_nperm,
  VT_ca50,
  VT_Tref,
  VT_nu,
  VT_mu,
  VT_xi,
  VT_a,
  VT_k,
  VT_eta_l,
  VT_eta_s,
    
  VT_k_ws,
  VT_k_uw,
  VT_cdw,
  VT_cds,
  VT_k_wu,
  VT_k_su,
  VT_A,
  VT_XSSS,
  VT_XWSS,
  VT_ktm_block,
    
  VT_XS_init,
  VT_XW_init,
  VT_TRPN_init,
  VT_TmBlocked_init,
  VT_ZETAS_init,
  VT_ZETAW_init,
  VT_Cd_init,

  vtLast
};
}  // namespace NS_TWorldParameters

using namespace NS_TWorldParameters;

class TWorldParameters : public vbNewElphyParameters {
 public:
  TWorldParameters(const char *, ML_CalcType);
  ~TWorldParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
};  // class TWorldParameters

#endif  // ifndef TWORLDPARAMETERS_H
