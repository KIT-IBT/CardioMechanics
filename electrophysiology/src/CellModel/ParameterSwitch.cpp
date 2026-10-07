/*
 * File: ParameterSwitch.cpp
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


#include <ParameterSwitch.h>

#include <set>
#include <typeinfo>

namespace {
set<string> warnedMissingParameters;
}  // namespace

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
    string warnKey = string(typeid(*stat).name()) + ":" + pDynPara.name;
    if (warnedMissingParameters.insert(warnKey).second) {
      cerr << "Warning: Parameter '" << pDynPara.name
           << "' was not defined in the current implementation of the cell model '"
           << typeid(*stat).name() << "'!" << endl;
    }
    return false;
  }
  useDynamicValues = true;
  return dyn[stat->P[index].dynamicVar] == pDynPara.value;
}  // ParameterSwitch::addDynamicParameter
