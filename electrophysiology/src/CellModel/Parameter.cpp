/*
 *  Parameter.cpp
 *  CellModel
 *
 *  Created by dw.local on 24.01.07.
 *  Copyright 2007 IBT Universitt Karlsruhe. All rights reserved.
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
