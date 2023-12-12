/*
 * File: Parameter.cpp
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


#include <Parameter.h>

Parameter::Parameter() {
  // cerr<<"loading Parameter()\n";
  name         = "UNDEFINED";
  value        = 0;
  readFromFile = true;
  dynamicVar   = -1;
}

Parameter::Parameter(string n, MYTYPEDEF v) {
  // cerr<<"loading Parameter("<<n.c_str()<<","<<v<<")\n";
  name         = n;
  value        = v;
  dynamicVar   = -1;
  readFromFile = true;
}
