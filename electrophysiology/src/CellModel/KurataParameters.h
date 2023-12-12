/*
 * File: KurataParameters.h
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



#ifndef KURATAPARAMETERS_H
#define KURATAPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_KurataParameters {
enum varType {
  VT_F = vtFirst,
  VT_R,
  VT_Tx,
  VT_nvar,
  VT_Nao,
  VT_Ko,
  VT_Cao,
  VT_C_m,
  VT_Vcell,
  VT_Vi,
  VT_Vrel,
  VT_Vup,
  VT_gCaLmax,
  VT_KmfCa,
  VT_ECaL,
  VT_gKrmax,
  VT_gKsmax,
  VT_PNaKs,
  VT_gtomax,
  VT_PNato,
  VT_gNamax,
  VT_PKNa,
  VT_gK1,
  VT_gbNa,
  VT_gbCa,
  VT_gbK,
  VT_iNaKmax,
  VT_KmNap,
  VT_nNa,
  VT_KmKp,
  VT_kNaCa,
  VT_KmNaex,
  VT_KmCaex,
  VT_rNaCa,
  VT_ksat,
  VT_ipCamax,
  VT_KmCap,
  VT_Prel,
  VT_nrel,
  VT_Pup,
  VT_Kup,
  VT_Pleak,
  VT_Ttr,
  VT_ConcTC,
  VT_kfTC,
  VT_kbTC,
  VT_ConcCM,
  VT_KdCM,
  VT_ConcCQ,
  VT_KdCQ,
  VT_Init_Vm,
  VT_Init_gdL,
  VT_Init_gfL,
  VT_Init_gpa,
  VT_Init_gn,
  VT_Init_gq,
  VT_Init_gh,
  VT_Init_Cai,
  VT_Init_Carel,
  VT_Init_Caup,
  VT_Init_gdR,
  VT_Init_gfR,
  VT_Init_Rtc,
  VT_Init_Nai,
  VT_Init_Ki,
  VT_RTdF,
  VT_FdRT,
  VT_RhoNaK,
  VT_CmdF,
  VT_CmdTwoF,
  VT_Amp,
  vtLast
};
}  // namespace NS_KurataParameters

using namespace NS_KurataParameters;

class KurataParameters : public vbNewElphyParameters {
 public:
  KurataParameters(const char *);
  ~KurataParameters();
  void PrintParameters();
  void Calculate();
  void InitTable();
  void Init(const char *);
  ML_CalcType gdLs[RTDT];
  ML_CalcType gfLs[RTDT];
  ML_CalcType gdLt[RTDT];
  ML_CalcType gfLt[RTDT];
  ML_CalcType gpas[RTDT];
  ML_CalcType gpat[RTDT];
  ML_CalcType gpis[RTDT];
  ML_CalcType gns[RTDT];
  ML_CalcType gnt[RTDT];
  ML_CalcType grs[RTDT];
  ML_CalcType gqs[RTDT];
  ML_CalcType gqt[RTDT];
  ML_CalcType gms[RTDT];
  ML_CalcType ghs[RTDT];
  ML_CalcType ght[RTDT];
  ML_CalcType gdRs[RTDT];
  ML_CalcType gfRs[RTDT];
};  // class KurataParameters
#endif  // ifndef KURATAPARAMETERS_H
