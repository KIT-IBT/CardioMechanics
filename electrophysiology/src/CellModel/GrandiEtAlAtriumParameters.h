/*
 * File: GrandiEtAlAtriumParameters.h
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


#ifndef GRANDIETALATRIUMPARAMETERS_H
#define GRANDIETALATRIUMPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_GrandiEtAlAtriumParameters {
enum varType {
  VT_stim_amplitude = vtFirst,
  VT_stim_duration,
  VT_R,
  VT_Cmem,
  VT_cellLength,
  VT_junctionLength,
  VT_distSLcyto,
  VT_distJuncSL,
  VT_DcaJuncSL,
  VT_DcaSLcyto,
  VT_DnaJuncSL,
  VT_DnaSLcyto,
  VT_J_ca_juncsl,
  VT_J_ca_slmyo,
  VT_J_na_juncsl,
  VT_J_na_slmyo,
  VT_Fjunc,
  VT_Fjunc_CaL,
  VT_K_o,
  VT_Na_o,
  VT_Ca_o,
  VT_Mgi,
  VT_AF,
  VT_ISO,
  VT_RA,
  VT_GNaL,
  VT_GNa,
  VT_GK1,
  VT_GNaB,
  VT_IbarNaK,
  VT_KmNaip,
  VT_KmKo,
  VT_Q10NaK,
  VT_Q10KmNai,
  VT_pNaK,
  VT_gkp,
  VT_GClCa,
  VT_GClB,
  VT_KdClCa,
  VT_pNa,
  VT_pCa,
  VT_pK,
  VT_Q10CaL,
  VT_IbarNCX,
  VT_KmCai,
  VT_KmCao,
  VT_KmNai,
  VT_KmNao,
  VT_ksat,
  VT_nu,
  VT_Kdact,
  VT_Q10NCX,
  VT_IbarSLCaP,
  VT_KmPCa,
  VT_GCaB,
  VT_Q10SLCaP,
  VT_Q10SRCaP,
  VT_GtoFast,
  VT_Gkur,
  VT_Vmax_SRCaP,
  VT_Kmf,
  VT_Kmr,
  VT_hillSRCaP,
  VT_ks,
  VT_koCa,
  VT_kom,
  VT_kim,
  VT_ec50SR,
  VT_Bmax_Naj,
  VT_Bmax_Nasl,
  VT_koff_na,
  VT_kon_na,
  VT_Bmax_TnClow,
  VT_koff_tncl,
  VT_kon_tncl,
  VT_Bmax_TnChigh,
  VT_koff_tnchca,
  VT_kon_tnchca,
  VT_koff_tnchmg,
  VT_kon_tnchmg,
  VT_Bmax_CaM,
  VT_koff_cam,
  VT_kon_cam,
  VT_Bmax_myosin,
  VT_koff_myoca,
  VT_kon_myoca,
  VT_koff_myomg,
  VT_kon_myomg,
  VT_Bmax_SR,
  VT_koff_sr,
  VT_kon_sr,
  VT_koff_sll,
  VT_kon_sll,
  VT_koff_slh,
  VT_kon_slh,
  VT_koff_csqn,
  VT_kon_csqn,
  VT_gks_junc,
  VT_gks_sl,
  VT_fcaCaMSL,
  VT_fcaCaj,
  VT_MaxSR,
  VT_MinSR,
  VT_Frdy,
  VT_cellRadius,
  VT_junctionRadius,
  VT_Fsl,
  VT_Fsl_CaL,
  VT_sigma,
  VT_gkr,
  VT_gkr0,
  VT_kiCa,
  VT_Temp,
  VT_Vcell,
  VT_SAjunc,
  VT_SAsl,
  VT_FoRT,
  VT_Qpow,
  VT_Vmyo,
  VT_Vsr,
  VT_Vsl,
  VT_Vjunc,
  VT_Cli,
  VT_Clo,
  VT_Bmax_SLlowsl,
  VT_Bmax_SLlowj,
  VT_Bmax_SLhighsl,
  VT_Bmax_SLhighj,
  VT_Bmax_Csqn,
  VT_ecl,
  VT_V_init,
  VT_Na_j_init,
  VT_Na_sl_init,
  VT_K_i_init,
  VT_Ca_j_init,
  VT_Ca_sl_init,
  VT_m_init,
  VT_h_init,
  VT_j_init,
  VT_ml_init,
  VT_hl_init,
  VT_x_kr_init,
  VT_x_ks_init,
  VT_Na_i_init,
  VT_x_kur_init,
  VT_y_kur_init,
  VT_x_to_f_init,
  VT_y_to_f_init,
  VT_d_init,
  VT_f_init,
  VT_f_Ca_Bj_init,
  VT_f_Ca_Bsl_init,
  VT_Ry_Rr_init,
  VT_Ry_Ro_init,
  VT_Ry_Ri_init,
  VT_Ca_sr_init,
  VT_Ca_i_init,
  VT_Na_Bj_init,
  VT_Na_Bsl_init,
  VT_Tn_CL_init,
  VT_Tn_CHc_init,
  VT_Tn_CHm_init,
  VT_CaM_init,
  VT_Myo_c_init,
  VT_Myo_m_init,
  VT_SRB_init,
  VT_SLL_j_init,
  VT_SLL_sl_init,
  VT_SLH_j_init,
  VT_SLH_sl_init,
  VT_Csqn_b_init,
  VT_powncx,
  VT_powcal,
  VT_powsrcap,
  VT_powslcap,
  VT_powKmPCa,
  VT_Ikiconst,
  vtLast
};
}  // namespace NS_GrandiEtAlAtriumParameters

using namespace NS_GrandiEtAlAtriumParameters;

class GrandiEtAlAtriumParameters : public vbNewElphyParameters {
 public:
  GrandiEtAlAtriumParameters(const char *, ML_CalcType);
  ~GrandiEtAlAtriumParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);
  ML_CalcType K_i[RTDT];
  ML_CalcType mss[RTDT];
  ML_CalcType exptaum[RTDT];
  ML_CalcType mlss[RTDT];
  ML_CalcType exptauml[RTDT];
  ML_CalcType xrss[RTDT];
  ML_CalcType exptauxr[RTDT];
  ML_CalcType xsss[RTDT];
  ML_CalcType exptauxs[RTDT];
  ML_CalcType xkurss[RTDT];
  ML_CalcType ykurss[RTDT];
  ML_CalcType exptauxkur[RTDT];
  ML_CalcType exptauykur[RTDT];
  ML_CalcType xtoss[RTDT];
  ML_CalcType ytoss[RTDT];
  ML_CalcType exptauxtof[RTDT];
  ML_CalcType exptauytof[RTDT];
  ML_CalcType dss[RTDT];
  ML_CalcType exptaud[RTDT];
  ML_CalcType fss[RTDT];
  ML_CalcType exptauf[RTDT];
  ML_CalcType exptauh[RTDT];
  ML_CalcType exptauj[RTDT];
  ML_CalcType hss[RTDT];
  ML_CalcType hlss[RTDT];
  ML_CalcType jss[RTDT];
  ML_CalcType fnak[RTDT];
  ML_CalcType rkr[RTDT];
  ML_CalcType kp_kp[RTDT];
  ML_CalcType I_Clbk[RTDT];
};  // class GrandiEtAlAtriumParameters

#endif  // ifndef GRANDIETALATRIUMPARAMETERS_H
