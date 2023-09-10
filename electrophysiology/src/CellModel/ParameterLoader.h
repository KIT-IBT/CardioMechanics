/*
 *  ParameterLoader.h
 *  CellModel
 *
 *  Created by dw.local on 05.02.07.
 *  Copyright 2007 IBT Universität Karlsruhe. All rights reserved.
 *
 */


#ifndef PARAMETERLOADER
#define PARAMETERLOADER

#include <ElphyModelBasis.h>
#include <ForceModelBasis.h>
#include <ParameterSwitch.h>

class ParameterLoader {
 private:
  string currentConfigFile;
  string defaultConfigFile;
  CellModelValues cmvCurrentConfig;
  CellModelValues cmvFullDefault;

 public:
  ParameterLoader(const char *initFile, ElphyModelType emt);
  ParameterLoader(const char *initFile, ForceModelType fmt);
  ML_CalcType getParameterValue(const char *desc, bool readFromFile);

  inline char *getOverlapString(void) {return cmvCurrentConfig.getOLString();}  // TODO (ew095): kann man prüfen, ob

  // cmvCurrentConfig hier eine gültige
  // Zeichenfolge liefert und ggf. zu
  // cmvFullDefault wechseln? Muss das
  // geprüft werden?
};
#endif  // ifndef PARAMETERLOADER
