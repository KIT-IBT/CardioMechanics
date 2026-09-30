/*
 * File: ElphyModelBasis.h
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


#ifndef ELPHYMODELBASIS_H
#define ELPHYMODELBASIS_H

#include <Parameter.h>
typedef MYTYPEDEF ML_CalcType;  // don't change to float!

#ifdef osMac_Vec

// #define USE_DYNAMIC_MULTIPLIER

# ifdef CVODE
 #  include "SharedObjects/CVODE_opt.h"
# endif  // ifdef CVODE

# include <kaSharedObjectHandling.h>
#endif  // ifdef osMac_Vec

namespace ElphyModelConstants {
const double F  = 96.485;     // kC/mol
const double R  = 8.314;      // J/(K mol)
const double k  = 1.3807e-23; // J/K
const double h  = 6.6261e-31; // kJ s
const double qe = 1.6022e-19;  // C
}

#include <CellModelValues.h>

const double d1d6 = 1.0/6.0;
const double d1d3 = 1.0/3.0;
const double d1d8 = 1.0/8.0;

#ifdef STIMULATION
const int RangeTab    = 5000;
const int DivisionTab = 10;
#else  // ifdef STIMULATION
const int RangeTab    = 400;
const int DivisionTab = 10;
#endif  // ifdef STIMULATION

const int RTDT            = RangeTab*DivisionTab;
const int RangeTabhalf    = RangeTab/2;
const double dDivisionTab = 1.0/DivisionTab;

typedef unsigned short Gates;
const double gatemult = 1.525878906e-5;
const double gatediv  = 1.0/gatemult;

enum ElphyModelType {
  EMT_Dummy = 0,
  EMT_BeelerReuter,
  EMT_CourtemancheEtAl,
  EMT_TenTusscher,
  EMT_SachseEtAl,
  EMT_TenTusscher2,
  EMT_MaleckarEtAl,
  EMT_KoivumaekiEtAl,
  EMT_GrandiEtAlVentricle,
  EMT_GrandiEtAlAtrium,
  EMT_FabbriEtAl,
  EMT_HimenoEtAl,
  EMT_OHaraRudy,
  EMT_OHaraRudyIso,
  EMT_Tomek,
  EMT_TWorld,
  EMT_MitchellSchaeffer,
  EMT_FitzhughNagumo,
  EMT_Last
};

static const char *ElphyModelDescriptor[] = {
  "ElphyDummy",
  "BeelerReuter",
  "CourtemancheEtAl",
  "TenTusscherEtAl",
  "SachseEtAl",
  "TenTusscher2",
  "MaleckarEtAl",
  "KoivumaekiEtAl",
  "GrandiEtAlVentricle",
  "GrandiEtAlAtrium",
  "FabbriEtAl",
  "HimenoEtAl",
  "OHaraRudy",
  "OHaraRudyIso",
  "Tomek",
  "TWorld",
  "MitchellSchaeffer",
  "FitzhughNagumo",
  "Last"
};


template<class T> class vbElphySharedObject;
template<class T> class SO_ionicCurrent;

inline const char *FindElphyModelDescriptor(ElphyModelType emt) {
  return ElphyModelDescriptor[emt];
}

inline ElphyModelType FindElphyModelType(string cmnr_str) {
  for (int emt = (int)EMT_Dummy; emt < (int)EMT_Last; emt++)
    if (cmnr_str == ElphyModelDescriptor[emt])
      return (ElphyModelType)emt;

  return EMT_Last;  //!< return EMT_Last if string is not valid
}

inline string FindElphyModelFileName(ElphyModelType emt) {
  if (emt == EMT_Dummy)
    return "";

  if (emt == EMT_Last)
    throw kaBaseException("Standard Elphy value file for %s is missing", FindElphyModelDescriptor(emt));

  static kaRootDir ElphyFiles("");
  string filename = ElphyFiles.GetData();
  filename.append(FindElphyModelDescriptor(emt));
  filename.append(".ev");
  return filename;
}

inline string FindElphyPictureFileName(ElphyModelType emt) {
  if (emt == EMT_Last)
    throw kaBaseException("Standard Elphy picture file for %s is missing", FindElphyModelDescriptor(emt));

  static kaRootDir ElphyFiles("");
  string filename = ElphyFiles.GetData();
  filename.append(FindElphyModelDescriptor(emt));
  filename.append(".xpm");
  return filename.c_str();
}

inline ElphyModelType FindElphyModelFileType(string cmnr_str) {
  CellModelValues mv;

  mv.HeaderCheck(cmnr_str.c_str());
  string ModelName = mv.getModelName();
  return FindElphyModelType(ModelName);
}

const ElphyModelType emt_Default = EMT_Dummy;

template<class T>
class vbElphyParameters : public CellModelMembers<T>, public nskaIPC::IPCShm {
 public:
  virtual int  GetSize(void)      = 0;
  virtual T   *GetBase(void)      = 0;
  virtual int  GetNumParameters() = 0;
  virtual void PrintParameters()  = 0;

  //  virtual void Init(const char * initFile)=0;
  // virtual void InitTable(void)=0; // Auskommentiert wegen Parameter tinc in Courtemanche
  virtual void Calculate(void) = 0;

  virtual void InitTableWithArrayStart(T *arrayStart) {}

  bool ReInitRequired;
};

class vbNewElphyParameters : public vbElphyParameters<ML_CalcType> {
 public:
  Parameter *P;

  virtual int GetSize(void) {return 0;}

  virtual ML_CalcType *GetBase(void) {return NULL;}

 private:
  virtual int GetNumParameters() {return -1;}
};

template<class T>
class vbElphyModel {
 public:
  vbElphyModel(void) {}

  virtual ~vbElphyModel(void) {}

  virtual void Init(void)                                 = 0;
  virtual void Print(ostream &tempstr, double t, T V)     = 0;
  virtual void LongPrint(ostream &tempstr, double t, T V) = 0;
  virtual void GetParameterNames(vector<string> &)        = 0;
  virtual void GetLongParameterNames(vector<string> &)    = 0;
  virtual T    GetVm(void)                                = 0;
  virtual T    GetCai(void)                               = 0;
  virtual T    GetCao(void)                               = 0;
  virtual T    GetNai(void)                               = 0;
  virtual T    GetNao(void)                               = 0;
  virtual T    GetKi(void)                                = 0;
  virtual T    GetKo(void)                                = 0;
  virtual int  GetSize(void)                              = 0;
  virtual T   *GetBase(void)                              = 0;
  virtual T    Volume(void)                               = 0; // in m^3
  virtual T    GetAmplitude(void)                         = 0; // in nA
  virtual T    GetStimTime(void)                          = 0; // in s

  virtual T GetCm(void) const  // Membrane capacitance per unit area in F/m^2
  {throw kaBaseException("GetCm not implemented in this CellModel.");}

  virtual T GetBeta(void) const  // Surface-to-volume ratio in 1/m
  {throw kaBaseException("GetBeta not implemented in this CellModel.");}

  virtual void SetCai(T) {}

  virtual void SetCao(T) {}

  virtual void SetNai(T) {}

  virtual void SetNao(T) {}

  virtual void SetKi(T) {}

  virtual void SetKo(T) {}

  virtual void SetTinc(T) {}

  virtual inline bool OutIsTCa() {return false;}

  virtual inline T GetTCa() {return 0;}

  virtual inline unsigned char getSpeed(T) {return 0;}

  virtual inline unsigned char setSpeed(T dVm, unsigned char speed) {
    const T adVm        = fabs(dVm);
    unsigned char stemp = getSpeed(adVm);

    if (stemp > speed)
      speed++;
    else
      speed = stemp;
    return (unsigned char)(speed == 5 ? 255 : (speed == 4 ? 108 : (speed == 3 ? 27 : (speed == 2 ? 6 : 1))));
  }

  inline void WriteStatus(FILE *writeto) {
    kaWrite(GetBase(), (int)sizeof(ML_CalcType), (size_t)GetSize()/sizeof(ML_CalcType), writeto);
  }

  inline void ReadStatus(FILE *readfrom) {
    kaRead(GetBase(), (int)sizeof(ML_CalcType), (size_t)GetSize()/sizeof(ML_CalcType), readfrom);
  }

  virtual inline void Set(vbElphyModel &from) {
    memcpy(GetBase(), from.GetBase(), from.GetSize());
  }

  virtual int GetNumStatus() {throw kaBaseException("GetNumStatus() not implemented.");}

  virtual void GetStatus(double *p) const {throw kaBaseException("GetStatus() not implemented.");}

  virtual void SetStatus(const double *p) {throw kaBaseException("SetStatus() not implemented.");}

  virtual T Calc(double tinc, T V, T i_external = 0., T stretch = 1., int euler = 1) = 0;

  virtual T CalcRungeKutta2(double tinc, T V, T i_external = 0., T stretch = 1., int euler = 1) {
    double status[255], status2[255], status3[255];

    GetStatus(status);
    double dV = Calc(tinc*.5, V, i_external, stretch, euler);
    GetStatus(status2);
    dV = Calc(tinc, V+dV, i_external, stretch, euler);
    GetStatus(status3);
    for (int i = 0; i < GetNumStatus(); i++)
      status[i] += (status3[i]-status2[i]);
    SetStatus(status);

    return dV;
  }

  virtual T CalcRungeKutta4(double tinc, T V, T i_external = 0., T stretch = 1., int euler = 1) {
    double status1[255], status2[255], status4[255], status6[255], status7[255], status8[255];

    GetStatus(status1);
    double dV1 = Calc(tinc*.5, V, i_external, stretch, euler);
    GetStatus(status2);
    double dV2 = Calc(tinc*.5, V+dV1, i_external, stretch, euler);
    GetStatus(status4);
    double dV3 = Calc(tinc, V+dV2, i_external, stretch, euler);
    GetStatus(status6);
    for (int i = 0; i < GetNumStatus(); i++)
      status7[i] = status1[i]+(status6[i]-status4[i]);
    SetStatus(status7);
    double dV4 = Calc(tinc, V+dV3, i_external, stretch, euler);
    GetStatus(status8);

    for (int i = 0; i < GetNumStatus(); i++)
      status1[i] += (status2[i]-status1[i])/3.
        + (status4[i]-status2[i])/1.5
        + (status6[i]-status4[i])/3.
        + (status8[i]-status7[i])/6.;
    SetStatus(status1);

    return dV1/3.+dV2/1.5+dV3/3.+dV4/6.;
  }

  inline void checkGatingVariable(T &gatingvariable) {
    if (gatingvariable < 0.0)
      gatingvariable = 0.0;
    else if (gatingvariable > 1.0)
      gatingvariable = 1.0;
  }

  virtual inline SO_ionicCurrent<T> *GetSO_ionicCurrent(int index) {
    throw kaBaseException("virtual inline SO_ionicCurrent<T> *GetSO_ionicCurrent HAS TO BE OVERWRITTEN!\n");
    return NULL;
  }

  virtual inline bool AddHeteroValue(string desc, double val) {
    cerr<<"AddHeteroValue is not defined in the corresponding cell model!\n";
    return false;
  }
};  // class vbElphyModel

#ifdef osMac_Vec

/*!<  Size of the array used to store global values such as Nai, Ki, RTdF, etc.
          Do not forget to increment the size if you add values below or everything will fail!!! */
const int ArraySize = 15;

template<class T>
class vbElphySharedObject {
 public:
  T *Vm;
  T *tinc_1000;
  T *Nai;
  T *Nao;
  T *Ki;
  T *Ko;
  T *Cai;
  T *Cao;
  T *RTdF;
  T *IKv14Na;
  T *CaSS;
  T *ICaK;
  int *Vi;
  T *C_m;
  T *CmdF;
  T *F;

  T tinc_Mul;
  int Multiplier;
  T tinc;

# ifdef USE_DYNAMIC_MULTIPLIER
  int oldStatusSince;
  bool recalcRequired;
  T tinc_old;
  T tinc_max;
  T tinc_min;
  int cnt;

  void checkMultiplier(T newValue, T bckValue, bool showText = false) {
    if (!recalcRequired) {
      if (abs((newValue-bckValue)/bckValue) > 10) {
        oldStatusSince = 0;
        if (showText)
          cerr<<"cnt="<<cnt<<": "<<newValue<<"\t"<<bckValue<<"\t"<<abs((newValue-bckValue)/bckValue)<<endl;
        recalcRequired = true;
      }
    }
    cnt++;
  }

  bool useCurrentResult(int minOldStatusCnt = 1000, T tincdiv = 10, T tincmul = 1.1) {
    if (!recalcRequired)
      oldStatusSince++;
    if (recalcRequired || ((tinc < tinc_max) && (oldStatusSince > minOldStatusCnt))) {
      bool useThisResult = false;
      if (recalcRequired) {
        useThisResult = (tinc == tinc_min);
        if (!useThisResult) {
          // cerr<<"runter:\ttinc from "<<tinc<<" to "<<tinc/tincdiv<<endl;
          tinc /= tincdiv;

          // tinc=tinc_min;
          // tinc/=2 -> 4 sec;
        }
      } else {
        // cerr<<"hoch:\ttinc from "<<tinc<<" to "<<tinc*tincmul<<endl;
        useThisResult = true;
        tinc         *= tincmul;

        // tinc*=1.1 -> 4 sec;
      }
      return useThisResult;
    }
    return true;
  }

# endif  // ifdef USE_DYNAMIC_MULTIPLIER

  vbElphySharedObject(void) {
# ifdef USE_DYNAMIC_MULTIPLIER
    MultiplierSet  = false;
    tinc_min       = 1.23E-10;
    tinc_old       = -1;
    tinc           = 0;
    oldStatusSince = 0;
# endif  // ifdef USE_DYNAMIC_MULTIPLIER
# if KADEBUG
    cerr<<"vbElphySharedObject()\n";
# endif  // if KADEBUG
  }

  virtual ~vbElphySharedObject() {
# if KADEBUG
    cerr<<"~vbElphySharedObject()\n";
# endif  // if KADEBUG
  }

  virtual void Init() = 0;

  virtual T *GetParentBase(void) {
    // cerr<<"vbElphySharedObject.GetBase\n";
    return Vm;
  }

  virtual int GetParentSize(void) {
    // cerr<<"vbElphySharedObject.GetSize\n";
    return (ICaK-Vm+1)*sizeof(T)+sizeof(int);
  }

  inline void SetParent(vbElphySharedObject *from) {
    // cerr<<"vbElphySharedObject.Set: copying " <<from->GetParentSize()<<" bytes from "<<from->GetParentBase()<<endl;
    memcpy(GetParentBase(), from->GetParentBase(), from->GetParentSize());
  }

  inline void SetMultiplier(T min_tinc, bool showText = TRUE, string src = "unknown") {
# ifdef USE_DYNAMIC_MULTIPLIER
    if (tinc == 0)
      tinc = min_tinc;
    if (tinc != tinc_old) {
      if (tinc < tinc_min)
        tinc = tinc_min;
      tinc_old = tinc;
# else  // ifdef USE_DYNAMIC_MULTIPLIER
    if (!MultiplierSet) {
# endif  // ifdef USE_DYNAMIC_MULTIPLIER
      Multiplier = (int)((*tinc_1000)/(min_tinc*1000));
      Multiplier = (Multiplier < 1 ? 1 : Multiplier);
      tinc_Mul   = *tinc_1000/Multiplier;
      if (showText)
        cerr<<src.c_str()<<":\tsetting dynamic multiplier to "<<Multiplier<<", tinc="<<tinc_Mul/1000<<" ...\n";
# if KADEBUG

      // cerr<<min_tinc<<"\t"<<*tinc_1000<<"\tMultiplier set to " <<Multiplier<<endl;
# endif  // if KADEBUG
      MultiplierSet = true;
    }
  }

  virtual void GetParameterNames(vector<string> &getpara) {
    getpara.push_back("GPN.ToDo4ThisObject");
  }

  virtual void GetLongParameterNames(vector<string> &getpara) {
    getpara.push_back("GLPN.ToDo4ThisObject");
  }

  virtual void Print(ostream &tempstr) {
    tempstr<<"Print:-69"<<' ';
  }

  virtual void LongPrint(ostream &tempstr) {
    tempstr<<"LongPrint:-69"<<' ';
  }

  virtual T CalcNext() = 0;

  /*! \fn void InitWithNewTinc( X newTinc )
   *  \brief Funiction to inform the shared ionic current about the tinc value set via the "-tinc" Parameter at command
   * line
   *  \param newTinc the new tinc value
   */
  virtual void InitWithNewTinc(T newTinc) {
    // cout << "vbElphySharedObject::InitWithNewTinc( " << newTinc << " ) called" << endl;
  }

  void SetPublicValues(T *ArrayStart, int *pVi) {
# if KADEBUG
    cerr<<"setting public values\n";
# endif  // if KADEBUG
    Vi        = pVi;
    Vm        = ArrayStart++; // 0
    tinc_1000 = ArrayStart++;  // 1
    Nai       = ArrayStart++; // 2
    Nao       = ArrayStart++; // 3
    Ki        = ArrayStart++; // 4
    Ko        = ArrayStart++; // 5
    Cai       = ArrayStart++; // 6
    Cao       = ArrayStart++; // 7
    RTdF      = ArrayStart++; // 8
    IKv14Na   = ArrayStart++; // 9
    CaSS      = ArrayStart++; // 10
    ICaK      = ArrayStart++; // 11
    C_m       = ArrayStart++; // 12
    CmdF      = ArrayStart++; // 13
    F         = ArrayStart;   // 14
  }

 private:
  bool MultiplierSet;
};  // class vbElphySharedObject


template<class T>
class vbElphySharedIonicCurrent : public vbElphySharedObject<T> {
 public:
  vbElphySharedIonicCurrent(void) {
# if KADEBUG
    cerr<<"vbElphySharedIonicCurrent()\n";
# endif  // if KADEBUG
  }

  ~vbElphySharedIonicCurrent() {
# if KADEBUG
    cerr<<"~vbElphySharedIonicCurrent()\n";
# endif  // if KADEBUG
  }

  void SetEquilibrium(T *value) {
    eq = value;
# if KADEBUG
    cerr<<"setting eq to " << *eq << endl;
# endif  // if KADEBUG
  }

  virtual T *GetBase(void) {
    // cerr<<"vbElphySharedIonicCurrent.GetBase\n";
    return 0;
  }

  virtual int GetSize(void) {
    // cerr<<"vbElphySharedIonicCurrent.GetSize\n";
    return 0;
  }

  void Set(vbElphySharedIonicCurrent<T> *from) {
    vbElphySharedObject<T>::SetParent(from);
    if (GetSize() == from->GetSize())
      memcpy(GetBase(), from->GetBase(), from->GetSize());
    else
      throw kaBaseException("ElphyModelBasis.h: GetSize!=from->GetSize\n");
  }

  T *eq;
};  // class vbElphySharedIonicCurrent

class vbElphySharedIonicCurrentNew : public vbElphySharedIonicCurrent<ML_CalcType> {
  ML_CalcType *dwTest;
};

template<class T>
class vbElphyConcentrationHandling : public vbElphySharedObject<T> {
 public:
  T *dVm;

  vbElphyConcentrationHandling(void) {
# if KADEBUG
    cerr<<"vbElphyConcentrationHandling()\n";
# endif  // if KADEBUG
  }

  ~vbElphyConcentrationHandling() {
# if KADEBUG
    cerr<<"~vbElphyConcentrationHandling()\n";
# endif  // if KADEBUG
  }
};

template<class T>
class vbElphyCaHandling : public vbElphyConcentrationHandling<T> {
 public:
  T *I_Ca;
  T *I_Cab;
  T *I_NaCa;
  T *I_pCa;

  vbElphyCaHandling(void) {
# if KADEBUG
    cerr<<"vbElphyCaHandling()\n";
# endif  // if KADEBUG
  }

  ~vbElphyCaHandling() {
# if KADEBUG
    cerr<<"~vbElphyCaHandling()\n";
# endif  // if KADEBUG
  }
};

template<class T>
class vbElphyNaHandling : public vbElphyConcentrationHandling<T> {
 public:
  T *I_Na;
  T *I_Nab;
  T *I_NaCa;
  T *I_NaK;

  vbElphyNaHandling(void) {
# if KADEBUG
    cerr<<"vbElphyNaHandling()\n";
# endif  // if KADEBUG
  }

  ~vbElphyNaHandling() {
# if KADEBUG
    cerr<<"~vbElphyNaHandling()\n";
# endif  // if KADEBUG
  }
};

template<class T>
class vbElphyKHandling : public vbElphyConcentrationHandling<T> {
 public:
  T *I_Kr;
  T *I_Ks;
  T *I_to;
  T *I_K1;
  T *I_NaK;
  T *I_CaK;
  T *I_Kb;
  T *I_Katp;
  T *i_Stim;
//i_SAC missing here? lh326

  vbElphyKHandling(void) {
# if KADEBUG
    cerr<<"vbElphyKHandling()\n";
# endif  // if KADEBUG
  }

  ~vbElphyKHandling() {
# if KADEBUG
    cerr<<"~vbElphyKHandling()\n";
# endif  // if KADEBUG
  }
};

enum SO_Type {
  ionicCurrent = 1,
  concHandling = 2
};


enum equilibrium_Type {
  eqE_Na = 1,
  eqE_to = 2,
  eqE_K  = 3,
  eqE_Ks = 4,
  eqE_Ca = 5
};
enum current_Type {
  Ito    = 1,
  INa    = 2,
  IKr    = 3,
  IKs    = 4,
  Iconst = 5,
  IK1    = 6,
  ICa    = 7,
  ICab   = 8,
  INaCa  = 9,
  INaK   = 10,
  INab   = 11,
  IKb    = 12,
  IpCa   = 13,
  iStim  = 14
};

enum concentration_Type {
  NaConc = 1,
  KConc  = 2,
  CaConc = 3
};

typedef current_Type (*c_Type)();
typedef concentration_Type (*conc_Type)();
typedef const char * (*Description_Type)();
typedef equilibrium_Type (*eq_Type)();
typedef SO_Type (*so_Type)();
typedef vbElphyParameters<ML_CalcType> *(*lParaType)(string evfile, ML_CalcType *ArrayBegin);
typedef vbElphySharedIonicCurrent<ML_CalcType> *(*lNew_Type)(string evfile, ML_CalcType *ArrayBegin, int *pVi);
typedef vbElphySharedIonicCurrent<ML_CalcType> *(*lNew_TypePara)(vbElphyParameters<ML_CalcType> *p,
                                                                 ML_CalcType *ArrayBegin, int *pVi);

typedef vbElphyCaHandling<ML_CalcType> *(*lCaHandling)(string evfile, ML_CalcType *ArrayBegin, int *pVi,
                                                       ML_CalcType *dVm, ML_CalcType *ICa, ML_CalcType *ICab,
                                                       ML_CalcType *INaCa, ML_CalcType *IpCa);
typedef vbElphyNaHandling<ML_CalcType> *(*lNaHandling)(string evfile, ML_CalcType *ArrayBegin, int *pVi,
                                                       ML_CalcType *dVm, ML_CalcType *INa, ML_CalcType *INab,
                                                       ML_CalcType *INaCa, ML_CalcType *INaK);
typedef vbElphyKHandling<ML_CalcType> *(*lKHandling)(string evfile, ML_CalcType *ArrayBegin, int *pVi, ML_CalcType *dVm,
                                                     ML_CalcType *IKr, ML_CalcType *IKs, ML_CalcType *Ito,
                                                     ML_CalcType *IK1, ML_CalcType *INaK, ML_CalcType *ICaK,
                                                     ML_CalcType *IKb, ML_CalcType *iStim);

typedef vbElphyNaHandling<ML_CalcType> *(*lNaHandlingPara)(vbElphyParameters<ML_CalcType> *p, ML_CalcType *ArrayBegin,
                                                           int *pVi, ML_CalcType *dVm, ML_CalcType *INa,
                                                           ML_CalcType *INab, ML_CalcType *INaCa, ML_CalcType *INaK);
typedef vbElphyKHandling<ML_CalcType> *(*lKHandlingPara)(vbElphyParameters<ML_CalcType> *p, ML_CalcType *ArrayBegin,
                                                         int *pVi, ML_CalcType *dVm, ML_CalcType *IKr, ML_CalcType *IKs,
                                                         ML_CalcType *Ito, ML_CalcType *IK1, ML_CalcType *INaK,
                                                         ML_CalcType *ICaK, ML_CalcType *IKb, ML_CalcType *IKatp,
                                                         ML_CalcType *iStim);
typedef vbElphyCaHandling<ML_CalcType> *(*lCaHandlingPara)(vbElphyParameters<ML_CalcType> *p, ML_CalcType *ArrayBegin,
                                                           int *pVi, ML_CalcType *dVm, ML_CalcType *ICa,
                                                           ML_CalcType *ICab, ML_CalcType *INaCa, ML_CalcType *IpCa);

class SharedObjectType {
 public:
  SharedObjectType(string fileName);
  ~SharedObjectType() { /*cerr<<"~SharedObjectType()\n";*/}

  so_Type GetSharedObjectType;

 private:
  void *handle;
};

// Definition SO-Object
template<class X>
class SO_Object {
 public:
  SO_Object(string bundleFileName, X *ArrayStart, int *pVi);
  SO_Object(void *handle);
  ~SO_Object(void);
  virtual X    CalcNext()                                     = 0;
  virtual void Print(ostream &tempstr)                        = 0;
  virtual void LongPrint(ostream &tempstr)                    = 0;
  virtual void GetParameterNames(vector<string> &getpara)     = 0;
  virtual void GetLongParameterNames(vector<string> &getpara) = 0;
  void *handle;
  Description_Type pDesc;
};


// Ionenstrom

template<class X>
class SO_ionicCurrent : public SO_Object<X> {
 public:
  SO_ionicCurrent(string bundleFileName, string evfile, X *ArrayStart, int *pVi);
  SO_ionicCurrent(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi);
  virtual ~SO_ionicCurrent();
  eq_Type GetEquilibriumType;
  void SetEquilibrium(X *value);
  void InitWithNewTinc(X newTinc);
  X    CalcNext();
  void Print(ostream &tempstr);
  void LongPrint(ostream &tempstr);
  void GetParameterNames(vector<string> &getpara);
  void GetLongParameterNames(vector<string> &getpara);

  virtual void Set(SO_ionicCurrent<X> *from) {
    pParaIC->Set(from->getParaIC());
  }

  vbElphySharedIonicCurrent<X> *getParaIC(void) {
    return pParaIC;
  }

 private:
  vbElphySharedIonicCurrent<X> *pParaIC;
};

// Ca-Handling
template<class X>
class SO_CaConc_Handling : public SO_Object<X> {
 public:
  SO_CaConc_Handling(string bundleFileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *ICa, X *ICab, X *INaCa,
                     X *IpCa);
  SO_CaConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm, X *ICa, X *ICab, X *INaCa,
                     X *IpCa);
  virtual ~SO_CaConc_Handling();
  X    CalcNext();
  void Print(ostream &tempstr);
  void LongPrint(ostream &tempstr);
  void GetParameterNames(vector<string> &getpara);
  void GetLongParameterNames(vector<string> &getpara);

 private:
  vbElphyCaHandling<X> *pParaConc;
};

// Na-Handling
template<class X>
class SO_NaConc_Handling : public SO_Object<X> {
 public:
  SO_NaConc_Handling(string bundleFileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *INa, X *INab, X *INaCa,
                     X *INaK);
  SO_NaConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm, X *INa, X *INab, X *INaCa,
                     X *INaK);
  virtual ~SO_NaConc_Handling();
  X    CalcNext();
  void Print(ostream &tempstr);
  void LongPrint(ostream &tempstr);
  void GetParameterNames(vector<string> &getpara);
  void GetLongParameterNames(vector<string> &getpara);

 private:
  vbElphyNaHandling<X> *pParaConc;
};

// K-Handling
template<class X>
class SO_KConc_Handling : public SO_Object<X> {
 public:
  SO_KConc_Handling(string bundleFileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *IKr, X *IKs, X *Ito,
                    X *IK1, X *INaK, X *ICaK, X *IKb, X *iStim);
  SO_KConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm, X *IKr, X *IKs, X *Ito,
                    X *IK1, X *INaK, X *ICaK, X *IKb, X *IKatp, X *iStim);
  virtual ~SO_KConc_Handling();
  X    CalcNext();
  void Print(ostream &tempstr);
  void LongPrint(ostream &tempstr);
  void GetParameterNames(vector<string> &getpara);
  void GetLongParameterNames(vector<string> &getpara);

 private:
  vbElphyKHandling<X> *pParaConc;
};

template<class X>
SO_Object<X>::~SO_Object(void) {
# if KADEBUG

  // printf("~SO_Object()\n");
# endif  // if KADEBUG
}

template<class X>
SO_Object<X>::SO_Object(string bundleFileName, X *ArrayStart, int *pVi) {
# if KADEBUG
  printf("SO_Object(): loading shLib from %s\n", bundleFileName.c_str());
# endif  // if KADEBUG
  handle = LoadLibrary((char *)bundleFileName.c_str());
  pDesc  = (Description_Type)LoadSymbol(this->handle, "Description");
}

template<class X>
SO_Object<X>::SO_Object(void *extHandle) {
  handle = extHandle;
  pDesc  = (Description_Type)LoadSymbol(this->handle, "Description");
# if KADEBUG
  cerr<<"SO_Object(): using shLib "<<extHandle<<" for " << pDesc()<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_ionicCurrent<X>::SO_ionicCurrent(string bundleFileName, string evfile, X *ArrayStart, int *pVi) : SO_Object<X>(
    bundleFileName, ArrayStart, pVi) {
# if KADEBUG
  printf("SO_ionicCurrent(): initializing %s\n", bundleFileName.c_str());
# endif  // if KADEBUG
  c_Type pcType = (c_Type)LoadSymbol(this->handle, "CurrentType");
  GetEquilibriumType = (eq_Type)LoadSymbol(this->handle, "EquilibriumType");
  lNew_Type plNewType = (lNew_Type)LoadSymbol(this->handle, "LoadNew");
# if KADEBUG
  printf("symbols loaded, evfile =%s ...\n", evfile.c_str());
# endif  // if KADEBUG
  pParaIC = plNewType(evfile, ArrayStart, pVi);
# if KADEBUG
  printf("ionicCurrent (%s): current loaded and initialized...\n", bundleFileName.c_str());

  // cerr<<"cV="<<*pParaIC->Vm<<"\t\tcVi="<<*pParaIC->Vi<<"\t\ttinc1000="<<*pParaIC->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_ionicCurrent<X>::SO_ionicCurrent(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi) : SO_Object<X>(
    handle) {
# if KADEBUG

  //    printf("SO_ionicCurrent(): initializing from handle %d and parameters\n",handle);
# endif  // if KADEBUG
  c_Type pcType = (c_Type)LoadSymbol(this->handle, "CurrentType");
  GetEquilibriumType = (eq_Type)LoadSymbol(this->handle, "EquilibriumType");
  lNew_TypePara plNewType = (lNew_TypePara)LoadSymbol(this->handle, "LoadNewPara");
# if KADEBUG
  cerr<<"symbols loaded, handle = " << handle<<endl;
# endif  // if KADEBUG
  pParaIC = plNewType(pP, ArrayStart, pVi);
# if KADEBUG
  printf("ionicCurrent (): current loaded and initialized...\n");

  // cerr<<"cV="<<*pParaIC->Vm<<"\t\tcVi="<<*pParaIC->Vi<<"\t\ttinc1000="<<*pParaIC->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_ionicCurrent<X>::~SO_ionicCurrent() {
# if KADEBUG

  // printf("~SO_ionicCurrent()\n");
# endif  // if KADEBUG
  delete pParaIC;
}

template<class X>
X SO_ionicCurrent<X>::CalcNext() {
  // cerr<<"Ausgabe SO_ionicCurrent \n";
  return pParaIC->CalcNext();
}

template<class X>
void SO_ionicCurrent<X>::Print(ostream &tempstr) {
  pParaIC->Print(tempstr);
}

template<class X>
void SO_ionicCurrent<X>::LongPrint(ostream &tempstr) {
  pParaIC->LongPrint(tempstr);
}

template<class X>
void SO_ionicCurrent<X>::GetParameterNames(vector<string> &getpara) {
  pParaIC->GetParameterNames(getpara);
}

template<class X>
void SO_ionicCurrent<X>::GetLongParameterNames(vector<string> &getpara) {
  pParaIC->GetLongParameterNames(getpara);
}

template<class X>
void SO_ionicCurrent<X>::SetEquilibrium(X *value) {
  pParaIC->SetEquilibrium(value);
}

/*! \fn void SO_ionicCurrent::InitWithNewTinc( X newTinc )
 *  \brief Funiction to inform the shared ionic current about the tinc value set via the "-tinc" Parameter at command
 * line
 *  \param newTinc the new tinc value
 */

template<class X>
void SO_ionicCurrent<X>::InitWithNewTinc(X newTinc) {
  pParaIC->InitWithNewTinc(newTinc);
}

/*  Ca - Handling  */

template<class X>
SO_CaConc_Handling<X>::SO_CaConc_Handling(string bundleFileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *ICa,
                                          X *ICab, X *INaCa, X *IpCa) : SO_Object<X>(bundleFileName, ArrayStart, pVi) {
# if KADEBUG

  // printf("SO_concentration(): initializing %s\n",fileName.c_str());
# endif  // if KADEBUG
  lCaHandling plNewType = (lCaHandling)LoadSymbol(this->handle, "Load");
  pParaConc = plNewType(evfile, ArrayStart, pVi, dVm, ICa, ICab, INaCa, IpCa);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<" fr "<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_CaConc_Handling<X>::SO_CaConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm,
                                          X *ICa, X *ICab, X *INaCa, X *IpCa) : SO_Object<X>(handle) {
# if KADEBUG

  // printf("SO_concentration(): initializing %s\n",fileName.c_str());
# endif  // if KADEBUG
  lCaHandlingPara plNewType = (lCaHandlingPara)LoadSymbol(this->handle, "LoadNewPara");
  pParaConc = plNewType(pP, ArrayStart, pVi, dVm, ICa, ICab, INaCa, IpCa);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<" fr "<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_CaConc_Handling<X>::~SO_CaConc_Handling() {
# if KADEBUG

  // printf("~SO_CaConc_Handling()\n");
# endif  // if KADEBUG
  delete pParaConc;
}

template<class X>
X SO_CaConc_Handling<X>::CalcNext() {
  return pParaConc->CalcNext();
}

template<class X>
void SO_CaConc_Handling<X>::Print(ostream &tempstr) {
  pParaConc->Print(tempstr);
}

template<class X>
void SO_CaConc_Handling<X>::LongPrint(ostream &tempstr) {
  pParaConc->LongPrint(tempstr);
}

template<class X>
void SO_CaConc_Handling<X>::GetParameterNames(vector<string> &getpara) {
  pParaConc->GetParameterNames(getpara);
}

template<class X>
void SO_CaConc_Handling<X>::GetLongParameterNames(vector<string> &getpara) {
  pParaConc->GetLongParameterNames(getpara);
}

/*
   virtual void Print(ostream &tempstr, double tArg, T V) {
    tempstr<<tArg<<' '<<Vm<<' '
           ; // -18
   };
 */

/*  Na - Handling  */
template<class X>
SO_NaConc_Handling<X>::SO_NaConc_Handling(string fileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *INa,
                                          X *INab, X *INaCa, X *INaK) : SO_Object<X>(fileName, ArrayStart, pVi) {
# if KADEBUG

  // printf("SO_concentration(): initializing %s\n",fileName.c_str());
# endif  // if KADEBUG
  lNaHandling plNewType = (lNaHandling)LoadSymbol(this->handle, "Load");
  pParaConc = plNewType(evfile, ArrayStart, pVi, dVm, INa, INab, INaCa, INaK);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_NaConc_Handling<X>::SO_NaConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm,
                                          X *INa, X *INab, X *INaCa, X *INaK) : SO_Object<X>(handle) {
# if KADEBUG
  cerr<<"SO_NaConc_Handling: handle="<<handle<<endl;
# endif  // if KADEBUG
  lNaHandlingPara plNewType = (lNaHandlingPara)LoadSymbol(this->handle, "LoadNewPara");
  pParaConc = plNewType(pP, ArrayStart, pVi, dVm, INa, INab, INaCa, INaK);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_NaConc_Handling<X>::~SO_NaConc_Handling() {
# if KADEBUG

  // printf("~SO_NaConc_Handling()\n");
# endif  // if KADEBUG
  delete pParaConc;
}

template<class X>
X SO_NaConc_Handling<X>::CalcNext() {
  return pParaConc->CalcNext();
}

template<class X>
void SO_NaConc_Handling<X>::Print(ostream &tempstr) {
  pParaConc->Print(tempstr);
}

template<class X>
void SO_NaConc_Handling<X>::LongPrint(ostream &tempstr) {
  pParaConc->LongPrint(tempstr);
}

template<class X>
void SO_NaConc_Handling<X>::GetParameterNames(vector<string> &getpara) {
  pParaConc->GetParameterNames(getpara);
}

template<class X>
void SO_NaConc_Handling<X>::GetLongParameterNames(vector<string> &getpara) {
  pParaConc->GetLongParameterNames(getpara);
}

/*  K - Handling  */
template<class X>
SO_KConc_Handling<X>::SO_KConc_Handling(string fileName, string evfile, X *ArrayStart, int *pVi, X *dVm, X *IKr, X *IKs,
                                        X *Ito, X *IK1, X *INaK, X *ICaK, X *IKb, X *iStim) : SO_Object<X>(fileName,
                                                                                                           ArrayStart,
                                                                                                           pVi) {
# if KADEBUG

  // printf("SO_concentration(): initializing %s\n",fileName.c_str());
# endif  // if KADEBUG
  lKHandling plNewType = (lKHandling)LoadSymbol(this->handle, "Load");
  pParaConc = plNewType(evfile, ArrayStart, pVi, dVm, IKr, IKs, Ito, IK1, INaK, ICaK, IKb, iStim);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_KConc_Handling<X>::SO_KConc_Handling(void *handle, vbElphyParameters<X> *pP, X *ArrayStart, int *pVi, X *dVm, X *IKr,
                                        X *IKs, X *Ito, X *IK1, X *INaK, X *ICaK, X *IKb, X *IKatp,
                                        X *iStim) : SO_Object<X>(handle) {
# if KADEBUG

  // printf("SO_concentration(): initializing %s\n",fileName.c_str());
# endif  // if KADEBUG
  lKHandlingPara plNewType = (lKHandlingPara)LoadSymbol(this->handle, "LoadNewPara");
  pParaConc = plNewType(pP, ArrayStart, pVi, dVm, IKr, IKs, Ito, IK1, INaK, ICaK, IKb, IKatp, iStim);
# if KADEBUG

  // printf("concentration handling (%s): loaded and initialized...\n",pDesc().c_str());
  // cerr<<"cV="<<*pParaConc->Vm<<"\t\tcVi="<<*pParaConc->Vi<<"\t\ttinc1000="<<*pParaConc->tinc_1000<<endl;
# endif  // if KADEBUG
}

template<class X>
SO_KConc_Handling<X>::~SO_KConc_Handling() {
# if KADEBUG

  // printf("~SO_KConc_Handling()\n");
# endif  // if KADEBUG
  delete pParaConc;
}

template<class X>
X SO_KConc_Handling<X>::CalcNext() {
  return pParaConc->CalcNext();
}

template<class X>
void SO_KConc_Handling<X>::Print(ostream &tempstr) {
  pParaConc->Print(tempstr);
}

template<class X>
void SO_KConc_Handling<X>::LongPrint(ostream &tempstr) {
  pParaConc->LongPrint(tempstr);
}

template<class X>
void SO_KConc_Handling<X>::GetParameterNames(vector<string> &getpara) {
  pParaConc->GetParameterNames(getpara);
}

template<class X>
void SO_KConc_Handling<X>::GetLongParameterNames(vector<string> &getpara) {
  pParaConc->GetLongParameterNames(getpara);
}

#endif  // ifdef osMac_Vec

#endif  // ifndef ELPHYMODELBASIS_H
