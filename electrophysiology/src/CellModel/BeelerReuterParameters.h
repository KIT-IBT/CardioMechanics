/*
 * File: BeelerReuterParameters.h
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


#ifndef BEELER_REUTER_PARAMETERS
#define BEELER_REUTER_PARAMETERS

#include <ParameterLoader.h>

namespace NS_BeelerReuterParameters {
enum varType {
  VT_C_m = vtFirst,
  VT_g_Na,
  VT_g_NaC,
  VT_E_Na,
  VT_g_s,
  VT_Vol,
  VT_Amp,
  VT_dC_m,
  VT_Init_Ca_i,
  VT_Init_m,
  VT_Init_h,
  VT_Init_j,
  VT_Init_d,
  VT_Init_f,
  VT_Init_x1,
  VT_Init_Vm,
  VT_RC_ax1C1,
  VT_RC_ax1C2,
  VT_RC_ax1C3,
  VT_RC_ax1C4,
  VT_RC_ax1C5,
  VT_RC_ax1C6,
  VT_RC_ax1C7,
  VT_RC_bx1C1,
  VT_RC_bx1C2,
  VT_RC_bx1C3,
  VT_RC_bx1C4,
  VT_RC_bx1C5,
  VT_RC_bx1C6,
  VT_RC_bx1C7,
  VT_RC_amC1,
  VT_RC_amC2,
  VT_RC_amC3,
  VT_RC_amC4,
  VT_RC_amC5,
  VT_RC_amC6,
  VT_RC_amC7,
  VT_RC_bmC1,
  VT_RC_bmC2,
  VT_RC_bmC3,
  VT_RC_bmC4,
  VT_RC_bmC5,
  VT_RC_bmC6,
  VT_RC_bmC7,
  VT_RC_ahC1,
  VT_RC_ahC2,
  VT_RC_ahC3,
  VT_RC_ahC4,
  VT_RC_ahC5,
  VT_RC_ahC6,
  VT_RC_ahC7,
  VT_RC_bhC1,
  VT_RC_bhC2,
  VT_RC_bhC3,
  VT_RC_bhC4,
  VT_RC_bhC5,
  VT_RC_bhC6,
  VT_RC_bhC7,
  VT_RC_ajC1,
  VT_RC_ajC2,
  VT_RC_ajC3,
  VT_RC_ajC4,
  VT_RC_ajC5,
  VT_RC_ajC6,
  VT_RC_ajC7,
  VT_RC_bjC1,
  VT_RC_bjC2,
  VT_RC_bjC3,
  VT_RC_bjC4,
  VT_RC_bjC5,
  VT_RC_bjC6,
  VT_RC_bjC7,
  VT_RC_adC1,
  VT_RC_adC2,
  VT_RC_adC3,
  VT_RC_adC4,
  VT_RC_adC5,
  VT_RC_adC6,
  VT_RC_adC7,
  VT_RC_bdC1,
  VT_RC_bdC2,
  VT_RC_bdC3,
  VT_RC_bdC4,
  VT_RC_bdC5,
  VT_RC_bdC6,
  VT_RC_bdC7,
  VT_RC_afC1,
  VT_RC_afC2,
  VT_RC_afC3,
  VT_RC_afC4,
  VT_RC_afC5,
  VT_RC_afC6,
  VT_RC_afC7,
  VT_RC_bfC1,
  VT_RC_bfC2,
  VT_RC_bfC3,
  VT_RC_bfC4,
  VT_RC_bfC5,
  VT_RC_bfC6,
  VT_RC_bfC7,
  VT_RC_ix1C1,
  VT_RC_ix1C2,
  VT_RC_ix1C3,
  VT_RC_ix1C4,
  VT_RC_ix1C5,
  VT_RC_ix1C6,
  VT_RC_ix1C7,
  vtLast
};
}  // namespace NS_BeelerReuterParameters

using namespace NS_BeelerReuterParameters;

class BeelerReuterParameters : public vbNewElphyParameters {
 public:
  BeelerReuterParameters(const char *, ML_CalcType);
  ~BeelerReuterParameters() {}

  void PrintParameters();

  // virtual inline int GetSize(void) { return (&RC[6][12]-&C_m+1)*sizeof(ML_CalcType);};
  // virtual inline ML_CalcType* GetBase(void) { return v(VT_C_m);};
  // virtual int GetNumParameters() { return 105; };

  void Init(const char *, ML_CalcType);
  void Calculate();
  void InitTable(ML_CalcType);

  double x1_inf[RTDT];
  double exptau_x1[RTDT];
  double m_inf[RTDT];
  double exptau_m[RTDT];
  double h_inf[RTDT];
  double exptau_h[RTDT];
  double j_inf[RTDT];
  double exptau_j[RTDT];
  double d_inf[RTDT];
  double exptau_d[RTDT];
  double f_inf[RTDT];
  double exptau_f[RTDT];
  double i_K_1[RTDT];
  double i_x_1[RTDT];
};  // class BeelerReuterParameters

#endif  // ifndef BEELER_REUTER_PARAMETERS
