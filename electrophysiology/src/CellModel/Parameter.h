/*
 *  Parameter.h
 *  CellModel
 *
 *  Created by dw.local on 24.01.07.
 *  Copyright 2007 IBT Universität Karlsruhe. All rights reserved.
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
