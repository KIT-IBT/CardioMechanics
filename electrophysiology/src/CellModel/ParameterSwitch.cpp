/*
 *  ParameterSwitch.cpp
 *  CellModel
 *
 *  Created by dw.local on 24.01.07.
 *  Copyright 2007 IBT Universitt Karlsruhe. All rights reserved.
 *
 */

#include <ParameterSwitch.h>

ParameterSwitch::ParameterSwitch(vbNewElphyParameters *s, unsigned int vtLASTEntry) {
  // cerr<<"ParaSwitch() mit vbNewElphyParameters, last="<<vtLASTEntry<<"\n";
  stat             = s;
  dyn              = NULL;
  cnt              = 0;
  vtLAST           = vtLASTEntry;
  useDynamicValues = false;
}

ML_CalcType ParameterSwitch::getValue(int vt) {
  // ElphyParameter v=stat->P[vt];
  /*
      int dynVar=stat->P[vt].dynamicVar;
      ML_CalcType value=stat->P[vt].value;
      return useDynamicValues?(dynVar<0?value:dyn[dynVar]):value;
   */
  return useDynamicValues ? (stat->P[vt].dynamicVar <
                             0 ? stat->P[vt].value : dyn[stat->P[vt].dynamicVar]) : stat->P[vt].value;

  // cerr<<"getValue: stat->K_o"<<stat->K_o<<", *stat->K_o="<<*(stat->K_o)<<endl;
  // cerr<<"dynamicVar="<<stat->P[vt].dynamicVar<<endl;
  if (stat->P[vt].dynamicVar < 0) {
    // cerr<<"static ...\n";
    return stat->P[vt].value;
  } else {
    cerr<<"dynamic ...\t(staticValue="<<stat->P[vt].value<<")\n";
    if ((int)cnt < stat->P[vt].dynamicVar) {
      cerr<<"undefined in this parameterset - using static value ...\n";
      return stat->P[vt].value;
    } else {
      return dyn[stat->P[vt].dynamicVar];
    }
  }
}

bool ParameterSwitch::addDynamicParameter(Parameter pDynPara) {
  // cerr<<"bisher sind "<<cnt<<" dynamische Parameter angelegt ...\n";
  ML_CalcType *tmp;

  if (cnt > 0) {
    // copy entries from dyn to tmp
    tmp = new ML_CalcType[cnt];
    memcpy(tmp, dyn, sizeof(ML_CalcType)*cnt);
    delete[]dyn;
  }
  dyn = new ML_CalcType[cnt+1];
  if (cnt > 0) {
    memcpy(dyn, tmp, sizeof(ML_CalcType)*cnt);
    delete[]tmp;
  }
  dyn[cnt] = pDynPara.value;
  unsigned int index = vtFirst;
  for (unsigned int y = vtFirst; y < vtLAST; y++) {
    if (stat->P[y].name == pDynPara.name) {
      stat->P[y].dynamicVar = cnt;

      // cerr<<y<<" -> "<<cnt<<endl;
      index = y;
    }
  }
  cnt++;

  if (index == vtFirst) {
    throw kaBaseException("Parameter '%s' was not defined in the current implementation of the model!",
                          pDynPara.name.c_str());
  }

  // cerr<<"jetzt sind "<<cnt<<" dynamische Parameter angelegt ...\n";
  // for (int x=0;x<cnt+1;x++){
  // cerr<<"\t"<<x<<": "<<(float)dyn[x]<<endl;
  // }
  // if (!useDynamicValues)
  //    cerr<<"using dynamic values is set!\n";
  useDynamicValues = true;
  /*if (dyn[stat->P[index].dynamicVar]!=pDynPara.value)
      cerr<<dyn[stat->P[index].dynamicVar]<<"!="<<pDynPara.value<<endl;*/
  return dyn[stat->P[index].dynamicVar] == pDynPara.value;
}  // ParameterSwitch::addDynamicParameter
