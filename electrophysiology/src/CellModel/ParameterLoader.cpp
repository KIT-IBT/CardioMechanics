/*
 * File: ParameterLoader.cpp
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


#include "ParameterLoader.h"

ParameterLoader::ParameterLoader(const char *initFile, ElphyModelType emt) {
  currentConfigFile = initFile;
  cmvCurrentConfig.HeaderCheck(initFile, ElphyModelDescriptor[emt]);
  cmvCurrentConfig.Load(initFile);

  cmvFullDefault.HeaderCheck(FindElphyModelFileName(emt).c_str(), ElphyModelDescriptor[emt]);
  cmvFullDefault.Load(FindElphyModelFileName(emt).c_str());

  //    cerr<<"using "<<currentConfigFile<<" as main, " << FindElphyModelFileName(emt) << " as default configuration
  // file\n";
  defaultConfigFile = FindElphyModelFileName(emt);
}

ParameterLoader::ParameterLoader(const char *initFile, ForceModelType fmt) {
  currentConfigFile = initFile;
  cmvCurrentConfig.HeaderCheck(initFile, ForceModelDescriptor[fmt]);
  cmvCurrentConfig.Load(initFile, MT_FORCE);

  cmvFullDefault.HeaderCheck(FindForceModelFileName(fmt).c_str(), ForceModelDescriptor[fmt]);
  cmvFullDefault.Load(FindForceModelFileName(fmt).c_str(), MT_FORCE);

  defaultConfigFile = FindForceModelFileName(fmt);
}

ML_CalcType ParameterLoader::getParameterValue(const char *desc, bool readFromFile) {
  // cerr<<"loading value '"<<desc<<"' from "<<currentConfigFile.c_str()<<" ...\n";
  for (int x = 0; x < cmvCurrentConfig.getnumElements(); x++) {
    if (!strcmp(desc, cmvCurrentConfig.getconstName(x)))
      return cmvCurrentConfig.getconstValue(x);
  }

  //    if (readFromFile)
  //        cerr<<"value '"<<desc<<"' not found in "<<currentConfigFile<<", searching in "<<defaultConfigFile.c_str()<<"
  // ...\n";
  //  for (int x = 0; x < cmvFullDefault.getnumElements(); x++) {
  //    if (!strcmp(desc, cmvFullDefault.getconstName(x)))
  //      return cmvFullDefault.getconstValue(x);
  //  }
  if (readFromFile) {
    // throw kaBaseException("Value for the necessary parameter '%s' not found in parameter file!",desc);
    cerr<<"Value for necessary parameter "<<desc<<" not found in parameter file!\n";
  }
  return -1;
}
