/* File: MitchellSchaefferParameters.cpp
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#include <MitchellSchaefferParameters.h>
MitchellSchaefferParameters::MitchellSchaefferParameters(const char *initFile) {
  P = new Parameter[vtLast];
  Init(initFile);
}

MitchellSchaefferParameters::~MitchellSchaefferParameters() {}

void MitchellSchaefferParameters::PrintParameters() {
  cout << "MitchellSchaefferParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void MitchellSchaefferParameters::Init(const char *initFile) {
#if KADEBUG
  cerr << "Loading the MitchellSchaeffer parameter from " << initFile << " ...\n";
#endif // if KADEBUG

  P[VT_tau_in].name = "tau_in";
  P[VT_tau_open].name = "tau_open";
  P[VT_tau_close].name = "tau_close";
  P[VT_V_gate].name = "V_gate";
  P[VT_tau_out].name = "tau_out";
  P[VT_Vm_init].name = "Vm_init";
  P[VT_u_init].name = "u_init";
  P[VT_h_init].name = "h_init";
  P[VT_Amp].name = "Amp";
  P[VT_InitTableDone].name =       "InitTableDone";

  P[VT_InitTableDone].readFromFile = false;

  ParameterLoader EPL(initFile, EMT_MitchellSchaeffer);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  Calculate();
  InitTable();

#if KADEBUG
  cerr << "#Init() done ... \n";
#endif // if KADEBUG
} // MitchellSchaefferParameters::Init

void MitchellSchaefferParameters::Calculate() {
#if KADEBUG
  cerr << "#MitchellSchaefferParameters - Calculate ..." << endl;
#endif // if KADEBUG
}

void MitchellSchaefferParameters::InitTable() {}
