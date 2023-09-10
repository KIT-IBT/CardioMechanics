/**@file HybridModelParameters.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef HYBRIDMODEL_PARAMETERS
#define HYBRIDMODEL_PARAMETERS

#include <ParameterLoader.h>

namespace NS_HybridModelParameters {
enum varType {
  VT_AMATPa,
  VT_MATPa,
  VT_MADPPa,
  VT_AwMADPPa,
  VT_AsMADPPa,
  VT_AsMADPa,
  VT_AMADPa,
  VT_TCaa,
  VT_TMona,
  VT_MADPa,
  VT_Ma,
  VT_k_1,
  VT_k_m1,
  VT_k_2,
  VT_k_3,
  VT_k_m3,
  VT_k_4,
  VT_k_m4,
  VT_k_5,
  VT_k_m5,
  VT_k_6,
  VT_k_m6,
  VT_k_7,
  VT_k_8,
  VT_k_m8,
  VT_k_on,
  VT_k_off,
  VT_k_9,
  VT_k_10,
  VT_k_11,
  VT_k_12,
  VT_k_13,
  VT_k_14,
  VT_tm_on,
  VT_tm_off,
  VT_TCaMax,
  VT_TCaMin,
  VT_Fmax,
  VT_dFmax,
  VT_N_v,
  VT_v50,
  VT_TCa_stretch,
  VT_TMon_coop,
  VT_TMon_pow,
  VT_ATP,
  VT_F_physiol,
  VT_k7_base,
  VT_k7_stretch,
  VT_k7_force,
  VT_k5_stretch,
  VT_k5_xb,
  VT_detach_vel,
  vtLast
};
}  // namespace NS_HybridModelParameters

using namespace NS_HybridModelParameters;

class HybridModelParameters : public vbNewForceParameters {
 public:
  HybridModelParameters(const char *);
  ~HybridModelParameters() {}

  // virtual inline int GetSize(void){return (&detach_vel-&AMATPa)*sizeof(T);};
  // virtual inline T* GetBase(void){return AMATPa;};
  // virtual int GetNumParameters() { return 51; };
  void Init(const char *);
  void Calculate();
  void PrintParameters();
};

#endif  // ifndef HYBRIDMODEL_PARAMETERS
