/*
 *  ParameterSwitch.h
 *  CellModel
 *
 *  Created by dw.local on 24.01.07.
 *  Copyright 2007 IBT Universität Karlsruhe. All rights reserved.
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
