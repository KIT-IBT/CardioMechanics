/*! \file ForceSteadyState.cpp
   \brief Main functions for program ForceSteadyState

   \author fs, CVRTI - University of Utah
 */

#include <CellModelLayer.h>

PrintParameterModus PrintParameterMode = PrintParameterModeOff;

void printHelp(char *name) {
  cerr<<name<<" "<<endl
      <<"\t[-fvfile <file>]"<<endl;
  for (int fmt = FMT_Dummy+1; fmt != FMT_Last; fmt++)
    cerr << "\t[-" << FindForceModelDescriptor((ForceModelType)fmt) << "]"  << endl;

  exit(-1);
}

typedef double CalcType;

int main(int argc, char **argv) {
  ForceModelType fmt;
  string initFVFile = "";

  vbForceModel<CalcType> *pfm       = NULL;
  vbForceParameters<CalcType> *ppfm = NULL;

  if (argc < 2)
    printHelp(argv[0]);

  bool verbose = false;

  try {
    double stretchdefault = 1.;

    for (int ac = 1; ac < argc; ac++) {
      const char *parg = argv[ac];
      if (!strcasecmp("-h", parg) || !strcasecmp("-help", parg)) {
        printHelp(argv[0]);
      } else if (!strcasecmp("-verbose", parg) || !strcasecmp("-v", parg)) {
        verbose = true;
      } else if (!strcasecmp("-fvfile", parg) && (ac < argc) ) {
        initFVFile = argv[++ac];
      } else if ((*parg == '-') && (FindForceModelType(parg+1) != FMT_Dummy) ) {
        fmt = FindForceModelType(parg+1);
      } else if (!strcasecmp("-stretch", parg) && (ac+1 < argc) ) {
        stretchdefault = atof(argv[++ac]);
      } else if (!strcasecmp("-sl", parg) && (ac < argc) ) {
        double sl = atof(argv[++ac]);
        stretchdefault = sl/2.0;
      } else {
        throw kaBaseException("Unknown or incomplete option %s", parg);
      }
    }

    if (initFVFile.size())
      fmt = FindForceModelFileType(initFVFile.c_str());

    if (!initFVFile.size())
      initFVFile = FindForceModelFileName(fmt);
    if (verbose)
      cerr << argv[0] << ": ForceModelType " << fmt << " used in " << initFVFile << endl;

    initForceParameters(&ppfm, initFVFile.c_str(), fmt);
    initForceModel(&pfm, ppfm, fmt);

    double pCa = -1.5, pCaInc = 0.1, pCaEnd = 1.5+pCaInc/2.;
    for (; pCa <= pCaEnd; pCa += pCaInc) {
      double x = pow(10., pCa);
      pfm->steadyState(x, stretchdefault);
    }
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    delete pfm;
    delete ppfm;
    exit(-1);
  }

  delete pfm;
  delete ppfm;
  return 0;
}  // main
