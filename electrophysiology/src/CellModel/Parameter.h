/*
 * File: Parameter.h
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

#ifndef PARAMETER
#define PARAMETER

#include <kaMachineOS.h>

enum PrintParameterModus { PrintParameterModeOn = 1, PrintParameterModeOff = 0 };
extern PrintParameterModus PrintParameterMode;

typedef double MYTYPEDEF;
static const int vtFirst = 0;

class Parameter {
 public:
  Parameter();
  Parameter(string, MYTYPEDEF);
  string name;
  MYTYPEDEF value;
  bool readFromFile;
  int dynamicVar;
};
#endif  // ifndef PARAMETER
