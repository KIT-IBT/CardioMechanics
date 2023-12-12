/*
 * File: GetEVFile.cpp
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


#include <CellModel.h>

PrintParameterModus PrintParameterMode = PrintParameterModeOff;

void printHelp(char *name) {
  cerr<<"GetEVFile returns the *.ev file that corresponds to a certain cell model\n\n";
  cerr<<"\tusage: "<<name<<
    " <value>\t(value can be a model number (int),\n\t\t\t\t\t\ta model descriptor, or an *.ev file describing a model)"
      <<
    endl;
  exit(-1);
}

int main(int argc, char **argv) {
  if (argc != 2)
    printHelp(argv[0]);
  try {
    char *caption      = argv[1];
    ElphyModelType emt = EMT_Dummy;
    if (atoi(caption) > 0) {
      emt = (ElphyModelType)atoi(caption);
    } else {
      nskaGlobal::FileCheck FC;
      FC.Check(caption);
      if (FC.Access() && FC.Exist()) {
        emt = FindElphyModelFileType(caption);
        if (emt != EMT_Last) {
          // valid ev-file - return it
          cout << caption << endl;
          return 0;
        } else {
          emt = FindElphyModelType(caption);
        }
      } else {
        emt = FindElphyModelType(caption);
        if (emt == EMT_Last)
          throw kaBaseException("no model number or valid model name given and file %s does not exist!", caption);
      }
    }
    cout << FindElphyModelFileName(emt)<<endl;
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
