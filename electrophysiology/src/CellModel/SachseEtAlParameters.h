/*
 * File: SachseEtAlParameters.h
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


#ifndef SachseEtAl_H
#define SachseEtAl_H

#include <ParameterLoader.h>

namespace NS_SachseEtAlParameters {
enum varType {
  VT_Tx = vtFirst,
  VT_Vm,
  VT_Cm,
  VT_Vfibro,
  VT_Ki,
  VT_Ko,
  VT_Cai,
  VT_Cao,
  VT_C0Shaker,
  VT_C1Shaker,
  VT_C2Shaker,
  VT_C3Shaker,
  VT_C4Shaker,
  VT_OShaker,
  VT_PShaker,
  VT_Shakerkv,
  VT_Shakerkvm,
  VT_Shakerko,
  VT_Shakerkom,
  VT_Shakerzv,
  VT_Shakerzvm,
  VT_Shakerzom,
  VT_GKir,
  VT_aKir,
  VT_bKir,
  VT_Gb,
  VT_Eb,
  VT_Gstretch,
  VT_Estretch,
  VT_GBK,
  VT_L0BK,
  VT_zLBK,
  VT_J0BK,
  VT_zJBK,
  VT_KDBK,
  VT_CBK,
  VT_DBK,
  VT_EBK,
  VT_Amp,
  vtLast
};
}  // namespace NS_SachseEtAlParameters

using namespace NS_SachseEtAlParameters;

class SachseEtAlParameters : public vbNewElphyParameters {
 public:
  SachseEtAlParameters(const char *initFile);
  ~SachseEtAlParameters() {}

  void PrintParameters();
  void Init(const char *initFile);

  void InitTable() {}

  void Calculate() {}
};
#endif  // ifndef SachseEtAl_H
