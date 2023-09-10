/*! \file evReplaceValue.cpp
   \brief Replace values in for given token at start of line

   \author fs, CVRTI - University of Utah, 25 Nov 05
 */

#include <kaExceptions.h>
#include <vector>

using std::vector;

class evoperation {
 public:
  evoperation() {}

  ~evoperation() {}

  void set(int modeArg, char *nameArg, char *valueArg) {
    mode  = modeArg;
    name  = nameArg;
    value = valueArg;
  }

  int mode;
  char *name;
  char *value;
};

int main(int argc, char *argv[]) {
  if (argc <= 3) {
    cerr << argv[0] << " <evFileName>" << endl;
    cerr << "\t[-replace <ParameterName> <Value>]" << endl;
    cerr << "\t[-add <ParameterName> <Value>]" << endl;
    cerr << "\t[-mult <ParameterName> <Value>]" << endl;
    cerr << "\t[-replaceNumber <ParameterNumber> <Value>]" << endl;
    cerr << "\t[-multNumber <ParameterNumber> <Value>]" << endl;
    cerr << "\t[-addNumber <ParameterNumber> <Value>]" << endl;
    cerr << "\t[-col 2|3]" << endl;
    exit(-1);
  }

  try {
    evoperation evo;
    vector<evoperation> vevo;

    char *evFileName = argv[1];
    int   col        = 2;

    for (int ac = 2; ac < argc; ac++) {
      const char *parg = argv[ac];
      if (!strcasecmp("-replace", parg) && (ac+2 < argc)) {
        evo.set(0, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-add", parg) && (ac+2 < argc)) {
        evo.set(1, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-mult", parg) && (ac+2 < argc)) {
        evo.set(2, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-replaceNumber", parg) && (ac+2 < argc)) {
        evo.set(3, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-addNumber", parg) && (ac+2 < argc)) {
        evo.set(4, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-multNumber", parg) && (ac+2 < argc)) {
        evo.set(5, argv[ac+1], argv[ac+2]);
        vevo.push_back(evo);
        ac += 2;
      } else if (!strcasecmp("-col", parg) && (ac+1 < argc)) {
        col = atof(argv[++ac]);
      } else {
        throw kaBaseException("Unknown or incomplete option %s", parg);
      }
    }

    ifstream in(evFileName);
    if (!in.is_open())
      throw kaBaseException("Can't open file %s", evFileName);
    for (int line = 1; !in.eof(); line++) {
      char buf[256], pn[256], val1[256], val2[256];
      in.getline(buf, sizeof(buf));
      int rc = sscanf(buf, "%s %s %s", pn, val1, val2);
      if (rc >= col) {
        for (int i = 0; i < vevo.size(); i++) {
          bool replace = false;
          switch (vevo[i].mode) {
            case 0:
            case 1:
            case 2:
              replace = (!strcmp(pn, vevo[i].name));
              break;
            case 3:
            case 4:
            case 5:
              replace = (line == atoi(vevo[i].name));
              break;
          }

          if (replace) {
            switch (vevo[i].mode) {
              case 0:
              case 3:
                if (col == 2)
                  sprintf(buf, "%s %s", pn, vevo[i].value);
                else if (col == 3)
                  sprintf(buf, "%s %s %s", pn, val1, vevo[i].value);
                break;
              case 1:
              case 4:
                if (col == 2)
                  sprintf(buf, "%s %lg", pn, atof(vevo[i].value)+atof(val1));
                else if (col == 3)
                  sprintf(buf, "%s %s %lg", pn, val1, atof(vevo[i].value)+atof(val2));
                break;
              case 2:
              case 5:
                if (col == 2)
                  sprintf(buf, "%s %lg", pn, atof(vevo[i].value)*atof(val1));
                else if (col == 3)
                  sprintf(buf, "%s %s %lg", pn, val1, atof(vevo[i].value)*atof(val2));
                break;
            }
          }
        }
      }
      printf("%s\n", buf);
    }
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
