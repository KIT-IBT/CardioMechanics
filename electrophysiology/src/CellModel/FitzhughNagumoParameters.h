/*
 * File: FitzhughNagumoParameters.h
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


#ifndef FITZHUGHNAGUMOPARAMETERS_H
#define FITZHUGHNAGUMOPARAMETERS_H

#include <ParameterLoader.h>

namespace NS_FitzhughNagumoParameters {
enum varType {
  VT_alpha = vtFirst,
  VT_gamma,
  VT_epsilon,
  VT_v_init,
  VT_w_init,
  VT_Vm_init,
  VT_AMP,
  VT_InitTableDone,
  vtLast
};
}

using namespace NS_FitzhughNagumoParameters;

class FitzhughNagumoParameters : public vbNewElphyParameters {
 public:
  FitzhughNagumoParameters(const char *);
  ~FitzhughNagumoParameters();
  void PrintParameters();
  void Calculate();
  void InitTable();
  void Init(const char *);
};

#endif // ifndef FITZHUGHNAGUMOPARAMETERS_H
