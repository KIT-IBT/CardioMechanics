/*
 * File: FitzhughNagumoParameters.cpp
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


#include <FitzhughNagumoParameters.h>

FitzhughNagumoParameters::FitzhughNagumoParameters(const char *initFile) {
  P = new Parameter[vtLast];
  Init(initFile);
}

FitzhughNagumoParameters::~FitzhughNagumoParameters() {}

void FitzhughNagumoParameters::PrintParameters() {
  cout << "FitzhughNagumoParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void FitzhughNagumoParameters::Init(const char *initFile) {
#if KADEBUG
  cerr << "Loading the FitzhughNagumo parameter from " << initFile << " ...\n";
#endif // if KADEBUG

  P[VT_alpha].name = "alpha";
  P[VT_gamma].name = "gamma";
  P[VT_epsilon].name = "epsilon";
  P[VT_v_init].name = "v_init";
  P[VT_w_init].name = "w_init";
  P[VT_Vm_init].name = "Vm_init";
  P[VT_AMP].name = "AMP";
  P[VT_InitTableDone].name = "InitTableDone";

  P[VT_InitTableDone].readFromFile = false;

  ParameterLoader EPL(initFile, EMT_FitzhughNagumo);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  Calculate();
  InitTable();

#if KADEBUG
  cerr << "#Init() done ... \n";
#endif // if KADEBUG
} // FitzhughNagumoParameters::Init

void FitzhughNagumoParameters::Calculate() {
#if KADEBUG
  cerr << "#FitzhughNagumoParameters - Calculate ..." << endl;

#endif // if KADEBUG
}

void FitzhughNagumoParameters::InitTable() {}
