/* File: MitchellSchaefferParameters.h
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

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
