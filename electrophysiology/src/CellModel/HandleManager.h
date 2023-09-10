/**@file HandleManager.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef HANDLEMANAGER
#define HANDLEMANAGER
#include <midLayer.h>
#include <ElphyModelBasis.h>

#define MAX_EVENTRIES 20

template<class X>
class HandleManager {
 public:
  HandleManager(const char *initFile);
  ~HandleManager();

  void **getHandle() {return handle;}

  vbElphyParameters<X> **getpEP() {return pEP;}

  int num;
  int numCurrents;
  int numConcentrations;

  int indexNaConcentration;
  int indexKConcentration;
  int indexCaConcentration;

 private:
  vbElphyParameters<X> **pEP;
  void **handle;
  string *FN;
  string *EV;
  string *sFN;
  string *sEV;
  bool loadInitFile(const char *initFile);
  bool loadHandles(void);
};  // class HandleManager

template<class X>
HandleManager<X>::HandleManager(const char *initFile) {
  numCurrents          = 0;
  numConcentrations    = 0;
  num                  = 0;
  indexNaConcentration = 0;
  indexKConcentration  = 0;
  indexCaConcentration = 0;
#if KADEBUG
  cerr<<"using "<<initFile<<endl;
#endif  // if KADEBUG
  if (loadInitFile(initFile))
    if (!(loadHandles()))
      throw kaBaseException("error during load of the handles!");
}

template<class X>
HandleManager<X>::~HandleManager() {
  if (FN) {
    FN = sFN;
    delete[] FN;
  }
  if (EV) {
    EV = sEV;
    delete[] EV;
  }
  if (pEP)
    delete[] pEP;
  if (handle) {
    //          for (int x = 0;x<numCurrents+numConcentrations;x++){
    //                  cerr<<"x="<<x<<": delete handle[x] ="<<handle[x]<<endl;
    //            delete handle[x];
    //          }
    delete[] handle;
  }
}

template<class X>
bool HandleManager<X>::loadInitFile(const char *initFile) {
  nskaGlobal::FileCheck FC;
  static kaRootDir localroot("bin/macosx/bundles/");
  CellModelValues  cmv;

  cmv.HeaderCheck(initFile);
  ifstream fp(initFile);
  if (fp) {
    char instr[640];
    EV  = new string[MAX_EVENTRIES];
    FN  = new string[MAX_EVENTRIES];
    sEV = EV;
    sFN = FN;
    for (int x = 0; x < cmv.getnumElements()+3; x++) {
      fp.getline(instr, sizeof(instr));
    }
    while (!fp.eof()) {
      if ((EV-sEV) < MAX_EVENTRIES) {
        fp.getline(instr, sizeof(instr));
        char evstr[200];
        char bdlstr[200];
        int  cnt = sscanf(instr, "%s %s", evstr, bdlstr);
        if (cnt > 0) {
          FC.Check(evstr);
          if (FC.Access() && FC.Exist()) {
            *EV = evstr;
            if (cnt == 2)
              *FN = bdlstr;
#if KADEBUG
            cerr<<"ev="<<EV<<", evstr="<<EV->c_str()<<", cnt="<<EV-sEV<<endl;
#endif  // if KADEBUG
            if (!(*FN->c_str())) {
              cmv.HeaderCheck(EV->c_str());
              string tmp = localroot.GetRoot();
              tmp += cmv.getModelName();
              tmp += ".bundle/Contents/MacOS/";
              tmp += cmv.getModelName();
              *FN  = tmp;
            } else {
              FC.Check(FN->c_str());
              if (FC.Access() && FC.Exist()) {
                cerr<<"using "<<FN->c_str()<<endl;
              } else {
                string tmp = localroot.GetRoot();
                tmp += FN->c_str();
                tmp += ".bundle/Contents/MacOS/";
                tmp += FN->c_str();
                *FN  = tmp;
              }
            }
            FC.Check(FN->c_str());
            if (!(FC.Access() && FC.Exist()))
              throw kaBaseException("can't access %s!", FN->c_str());
#if KADEBUG
            cerr<<"EV: "<<EV->c_str()<<"\tFN: "<<FN->c_str()<<endl;
#endif  // if KADEBUG
            EV++;
            FN++;
          } else {
            throw kaBaseException("can't access %s!", evstr);
          }
        }
      } else {
        fp.getline(instr, sizeof(instr));
        throw kaBaseException("can't use %s, MAX_EVENTRIES defined to %i in HandleManager.h!\n", instr, MAX_EVENTRIES);
      }
    }
    fp.close();
    return true;
  } else {
    throw kaBaseException("Error !fp\n");
  }
}  // >::loadInitFile

template<class X>
bool HandleManager<X>::loadHandles(void) {
  num    = EV-sEV;
  handle = new void *[num];
  pEP    = new vbElphyParameters<X> *[num];
  EV     = sEV;
  FN     = sFN;
  while (EV-sEV < num) {
    handle[EV-sEV] = LoadLibrary((char *)FN->c_str());
    Description_Type pDesc = (Description_Type)LoadSymbol(handle[EV-sEV], "Description");
    lParaType pPara        = (lParaType)LoadSymbol(handle[EV-sEV], "LoadParameter");
    so_Type   ot           = (so_Type)LoadSymbol(handle[EV-sEV], "SOType");
    if (ot() == ionicCurrent) {
      numCurrents++;
      /*                        c_Type ct=(c_Type)LoadSymbol(handle[EV-sEV],"CurrentType");
                              if (ct()==INaCa)
                                      cerr<<"INaCa in "<<EV-sEV<<endl;
       */
    } else {
      numConcentrations++;
      conc_Type cot = (conc_Type)LoadSymbol(handle[EV-sEV], "ConcentrationType");
      if (cot() == NaConc)
        indexNaConcentration = EV-sEV;
      else if (cot() == KConc)
        indexKConcentration = EV-sEV;
      else if (cot() == CaConc)
        indexCaConcentration = EV-sEV;
      else
        throw kaBaseException("unknown ConcentrationType '%i'!", cot());
    }
    pEP[EV-sEV] = pPara(*EV, 0);
#if KADEBUG
    cerr<<EV-sEV<<" - "<<pDesc()<<"\tParameters "<<EV->c_str()<<" in "<<FN->c_str()<<", adress="<<handle[EV-sEV]<<endl;
#endif  // if KADEBUG
    EV++; FN++;
  }
#if KADEBUG
  cerr<<"considering "<<numCurrents<<" ionic currents and "<<numConcentrations<<" concentrations ...\n";
#endif  // if KADEBUG
  return (numCurrents+numConcentrations) == num;

  // delete[] EV;
  // delete[] FN;
}  // >::loadHandles

#endif  // ifndef HANDLEMANAGER
