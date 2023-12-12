/*
 * File: ParameterLoader.h
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

  inline char *getOverlapString(void) {return cmvCurrentConfig.getOLString();}  // TODO (ew095): kann man prÃ¼fen, ob

  // cmvCurrentConfig hier eine gÃ¼ltige
  // Zeichenfolge liefert und ggf. zu
  // cmvFullDefault wechseln? Muss das
  // geprÃ¼ft werden?
};
#endif  // ifndef PARAMETERLOADER
