/*
 * File: MitchellSchaefferParameters.h
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


#ifndef MITCHELLSCHAEFFERPARAMETERS_H
#define MITCHELLSCHAEFFERPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_MitchellSchaefferParameters {
enum varType {
  VT_tau_in = vtFirst,
  VT_tau_open,
  VT_tau_close,
  VT_V_gate,
  VT_tau_out,
  VT_Vm_init,
  VT_u_init,
  VT_h_init,
  VT_Amp,
  VT_InitTableDone,
  vtLast
};
}

using namespace NS_MitchellSchaefferParameters;

class MitchellSchaefferParameters : public vbNewElphyParameters {
 public:
  MitchellSchaefferParameters(const char *);
  ~MitchellSchaefferParameters();
  void PrintParameters();
  void Calculate();
  void InitTable();
  void Init(const char *);
};

#endif // ifndef MITCHELLSCHAEFFERPARAMETERS_H
