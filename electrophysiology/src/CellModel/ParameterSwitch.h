/*
 * File: ParameterSwitch.h
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


#ifndef PARAMETERSWITCH
#define PARAMETERSWITCH

#include <ElphyModelBasis.h>

class ParameterSwitch {
 public:
  ParameterSwitch(vbNewElphyParameters *s, unsigned int);
  ML_CalcType getValue(int vt);
  bool addDynamicParameter(Parameter);

 private:
  vbNewElphyParameters *stat;
  ML_CalcType *dyn;
  unsigned int cnt;
  unsigned int vtLAST;
  bool useDynamicValues;
};


#endif  // ifndef PARAMETERSWITCH
