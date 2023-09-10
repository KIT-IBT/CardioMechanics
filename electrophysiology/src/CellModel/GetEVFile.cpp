/*
 *  GetEVFile.cpp
 *  CellModel
 *
 *  Created by dw on 20.07.07.
 *  Copyright 2007 IBT Uni Karlsruhe (TH) All rights reserved.
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
