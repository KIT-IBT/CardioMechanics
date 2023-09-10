/*! \file Drug2EV.cpp
   \brief Functions for program Drug2EV

   \version 1.0.1

   \date Created Christian Rombach (24.06.09) \n

   \author      Christian Rombach\n
   Institute of Biomedical Engineering\n
   Universitaet Karlsruhe (TH)\n
   http://www.ibt.uni-karlsruhe.de\n
 */

// doxygen manpage Drug2EV
/*! \page Drug2EV Drug2EV
   Modify the evfile of a cellmodel to simulate the impact of drugs on the conductance of the ion channels

   \section SYNOPSIS_Drug2EV SYNOPSIS
   Drug2EV \<evfile\> \<datafile\>

   \section OPTIONS_Drug2EV OPTIONS
   \param "<evfile>" Original evfile with the initvalues and standardparameters of any Cellmodel.
   \param "<datafile>" File with the specific data of the applied drug. Each line contains a conductance parameter, the
      corresponding IC50 value and the Hillcoefficient nH.
   \param "-concentration" Applied concentration of the drug.
   \param "-bound" Percentage of the drug that is bound to plasma proteins.
   \param "-outfile" Name of the modified evfile.
   \param "-verbose" Additional informations about the performed modifications.

   \section DESCRIPTION_Drug2EV DESCRIPTION
   Drug2EV manipulates the conductance values of any cellmodel. The conductances will be multiplied with a hillfunction
      that
   contains the parameters given in the datafile and the concentration of the specific drug.\n
   The parameters in the datafile should look like this: \<conductance parameter\> \<ic50\> \<nh\>.\n
   Lines that should be ignored can be commented with "//"\n.

   \section SOURCE_Drug2EV SOURCE
   Drug2EV.cpp

 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <kaParser.h>
#include <vector>
#include <kaVersion.h>
#include "Drug2EV.h"

using namespace std;

// define hill equation
double hill(double ic50, double nh, double concentration) {
  return 1/(1+pow((ic50/concentration), nh));
}

int main(int argc, char **argv) {
  //! read the original evfile and create a new one
  kaParser CL(argc, argv);
  double   conc_old, value, value_mod;
  double   bound         = 0;
  double   concentration = 1;
  int num_parameters     = 0;
  int is_value[500];
  char inputfile[64], outputfile[64], parameterfile[64], zeile[4096];
  string temp_string, var_name;
  char   f_open;
  int para_set;
  int line = 0, found = 0, j, k = 0;
  bool verbose;

  if (argc < 3) {
    cout << "Drug2EV <evfile> <datafile>" << endl;
    kaVersion v(1, 0, 1, 2009, 8, 7, argc, argv);
    v.addOption("", "evfile", OT_required, "",
                "\t\t\t\t\t original evfile with the initvalues and standardparameters of any Cellmodel");
    v.addOption("", "datafile", OT_required, "",
                "\t\t\t\t\t file with the specific data of the applied drug. Each line contains a conductance parameter, the corresponding\n\t\t\t\t\t\t\t IC50 value and the Hillcoefficient nH (<g_X> <IC50> <nH>)");
    v.addOption("concentration", "", OT_optional, "1", "\t\t\t applied concentration of the drug");
    v.addOption("bound", "", OT_optional, "0", "\t\t\t\t percentage of the drug that is bound to plasma proteins");
    v.addOption("outfile", "", OT_optional, "<inputfile>_mod.ev", "\t name of the modified evfile.");
    v.addOption("verbose", "", OT_optional, "", "\t\t\t\t\t additional informations about the performed modifications");
    v.printHelp();
    return 1;
  } else {
    strcpy(inputfile, argv[1]);
    strcpy(parameterfile, argv[2]);
  }


  // check inputparameters
  if (CL.IsOption("-verbose")) {
    verbose = true;
  }
  if (CL.IsOption("-outfile")) {
    sscanf(CL.GetOption("-outfile", 1), "%s", outputfile);
  } else {
    strcpy(outputfile, argv[1]);
    outputfile[strlen(inputfile)-3] = '_';
    outputfile[strlen(inputfile)-2] = 'm';
    outputfile[strlen(inputfile)-1] = 'o';
    outputfile[strlen(inputfile)+0] = 'd';
    outputfile[strlen(inputfile)+1] = '.';
    outputfile[strlen(inputfile)+2] = 'e';
    outputfile[strlen(inputfile)+3] = 'v';
    outputfile[strlen(inputfile)+4] = '\0';
  }
  if (CL.IsOption("-concentration")) {
    sscanf(CL.GetOption("-concentration", 1), "%lf", &concentration);
  }
  if (CL.IsOption("-bound")) {
    sscanf(CL.GetOption("-bound", 1), "%lf", &bound);
  }

  // check if inputfile exists
  ifstream eingabe_ev(inputfile, ios::in);
  if (!eingabe_ev.good()) {
    cout << "file " << inputfile << " not found" << endl;
    return 1;
  } else {}
  eingabe_ev.close();

  // check if datafile exists
  ifstream eingabe_para(parameterfile, ios::in);
  if (!eingabe_para.good()) {
    cout << "file " << parameterfile << " not found" << endl;
    return 1;
  } else {}
  eingabe_para.close();

  // check if outputfile allready exists
  ofstream ausgabe;
  ausgabe.open(outputfile, ios_base::in);
  ausgabe.close();
  if (ausgabe.good()) {
    cout << "overwrite " <<  outputfile << " (y/n)? : ";
    cin >> f_open;
    if (f_open != 'y')
      return 1;
  } else {}

  // check parameters in datafile
  eingabe_para.open(parameterfile, ios_base::in);
  eingabe_para.seekg(0L, ios::beg);

  while (!eingabe_para.eof()) {
    eingabe_para.getline(zeile, 4096);
    if ((zeile[0] == '/') && (zeile[1] == '/') ) {
      is_value[line] = 0;
    } else {
      num_parameters++;
      is_value[line] = 1;
    }
    line++;
  }
  eingabe_para.close();

  // initiate arrays with parameternames, ic50 and nh values
  vector<string> parameter_name(num_parameters, "");
  vector<double> ic50(num_parameters);
  vector<double> nh(num_parameters);

  eingabe_para.open(parameterfile, ios_base::in);
  eingabe_para.seekg(0L, ios::beg);

  line = 0;

  while (getline(eingabe_para, temp_string)) {
    istringstream ss(temp_string);
    if (is_value[line] == 1) {
      ss >> parameter_name[k];
      ss >> ic50[k];
      ss >> nh[k];
      k++;
    } else {}
    line++;
  }

  eingabe_para.close();

  // output
  if (verbose == true) {
    cout << "evfile: " << inputfile << endl;
    cout << "datafile: " << parameterfile << endl;
    cout << "\n";
    for (int i = 0; i < num_parameters; i++) {
      cout << parameter_name[i] << ": IC50 = " <<  ic50[i] << ", nH = " << nh[i] << endl;
    }
    cout << "\n";
    cout << "concentration: " << concentration << endl;
    cout << "bound to plasma proteins (%): " << bound << endl;
    cout << "effective concentration: " << concentration * ((100-bound)/100) << endl;
    cout << "\n";
  }

  // set effective concentration
  concentration = concentration * ((100-bound)/100);

  // create modified ev file
  eingabe_ev.open(inputfile, ios_base::in);
  ausgabe.open(outputfile, ios_base::out);
  eingabe_ev.seekg(0L, ios::beg);
  while (getline(eingabe_ev, temp_string)) {
    if (line > 2) {
      istringstream ss(temp_string);
      ss >> var_name;
      ss >> value;
      for (int i = 0; i < num_parameters; i++) {
        int k = var_name.compare(parameter_name[i]);
        if (k == 0) {
          found = 1;
          j     = i;
          break;
        } else {}
      }
      if (found == 1) {
        if (ic50[j] != 0) {
          value_mod = value*(1-hill(ic50[j], nh[j], concentration));
          if (verbose == true) {
            cout << "modify " << var_name << ": " << value << " -> " << value_mod << ", reduced to: " << value_mod/
              value*100 << "%" << endl;
          }
        } else {
          value_mod = value;
        }
        ausgabe << var_name << " " << value_mod << endl;
      } else {
        ausgabe << temp_string << endl;
      }
      found = 0;
    } else {
      ausgabe << temp_string << endl;
    }
    line++;
  }
  if (verbose == true) {
    cout << "\noutputfile: " << outputfile << endl;
  }
  eingabe_ev.close();
  ausgabe.close();

  return 0;
}  // main
