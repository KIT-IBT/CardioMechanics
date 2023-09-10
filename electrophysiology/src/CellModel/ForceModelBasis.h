/**@file ForceModelBasis.h
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

#ifndef FORCE_MODEL_BASIS_H
#define FORCE_MODEL_BASIS_H

#include <Parameter.h>
#include <CellModelValues.h>

typedef MYTYPEDEF ML_CalcType;  // don't change to float!

enum ForceModelType {
  FMT_Dummy            = 0,
  FMT_Hybrid           = 1,
  FMT_Land17           = 2,
  FMT_Last             = 3
};


static const char *ForceModelDescriptor[] = {
  "ForceDummy",
  "Hybrid",
  "Land17",
  "Last"
};

inline const char *FindForceModelDescriptor(ForceModelType fmt) {
  return ForceModelDescriptor[fmt];
}

inline ForceModelType FindForceModelType(string cmnr_str) {
  for (int fmt = (int)FMT_Dummy; fmt < (int)FMT_Last; fmt++)
    if (cmnr_str == ForceModelDescriptor[fmt])
      return (ForceModelType)fmt;

  return FMT_Last;  //!< return FMT_Last if string is not valid
}

inline string FindForceModelFileName(ForceModelType fmt) {
  if (fmt == FMT_Dummy)
    return "";

  if (fmt == FMT_Last)
    throw kaBaseException("Standard Elphy value file for %s is missing", FindForceModelDescriptor(fmt));

  static kaRootDir ForceFiles("");
  string filename = ForceFiles.GetData();
  filename.append(FindForceModelDescriptor(fmt));
  filename.append(".fv");
  return filename;
}

inline string FindForcePictureFileName(ForceModelType fmt) {
  switch (fmt) {
    case FMT_Last:
      throw kaBaseException("Standard force picture file for %s is missing", FindForceModelDescriptor(fmt));
    default:
      break;
  }

  static kaRootDir ForceFiles("");
  string filename = ForceFiles.GetData();
  filename.append(FindForceModelDescriptor(fmt));
  filename.append(".xpm");
  return filename;
}

inline ForceModelType FindForceModelFileType(string fmnr_str) {
  CellModelValues mv;

  mv.HeaderCheck(fmnr_str.c_str());
  string ModelName = mv.getModelName();
  return FindForceModelType(ModelName);
}

const ForceModelType fmt_Default = FMT_Dummy;


template<class T>
class overlapparameters {
 public:
  inline int getOverlapID() {return (int)*overlapID;}

  inline T *getOverlapParameters() {return *olparams;}

  inline void SetOLShmaddr(T *ShmAddr) {
    overlapID = ShmAddr;
    for (int i = 1; i < OVERLAP_LEN; i++) {
      olparams[i-1] = overlapID+i;
    }
  }

  inline void setOverlapParameters(const char *olpar) {
    float olID, olp[OVERLAP_LEN-1] = {.0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0};
    int   rc = sscanf(olpar, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
                      , &olID, &olp[0], &olp[1], &olp[2], &olp[3], &olp[4], &olp[5], &olp[6],
                      &olp[7], &olp[8], &olp[9], &olp[10], &olp[11], &olp[12], &olp[13], &olp[14]);

    // cerr<<olpar<<" "<<rc<<endl;

    if (rc == 0)
      *overlapID = 0;
    else
      *overlapID = olID;

    int i = 0;
    for (i = 0; i < OVERLAP_LEN-1; i++)
      *olparams[i] = .0;

    if ((*overlapID > 3) || (*overlapID == 0)) {
      for (i = 0; i < OVERLAP_LEN-1; i++)
        *olparams[i] = olpar[i];
    } else {
      *olparams[0] = (olp[1]-olp[3])/(olp[0]-olp[2]);
      *olparams[1] = olp[1]-*olparams[0]*olp[0];
      *olparams[2] = olp[2];
      *olparams[3] = (olp[3]-olp[5])/(olp[2]-olp[4]);
      *olparams[4] = olp[3]-*olparams[3]*olp[2];
      for (i = 3; i < rc-1; i++) {
        *olparams[3*i-4] = olp[(i-1)*2];
        *olparams[3*i-3] = (olp[(i-1)*2+1]-olp[i*2+1])/(olp[(i-1)*2]-olp[i*2]);
        *olparams[3*i-2] = olp[(i-1)*2+1] - *olparams[3*i-3]*olp[(i-1)*2];
      }
    }
  }  // setOverlapParameters

 private:
  T *overlapID;
  T *olparams[OVERLAP_LEN-1];
};  // class overlapparameters

template<class T>
class vbForceParameters : public overlapparameters<T>, public CellModelMembers<T>, public nskaIPC::IPCShm {
 public:
  virtual int  GetSize(void)          = 0;
  virtual T   *GetBase(void)          = 0;
  virtual int  GetNumParameters(void) = 0;
  virtual void PrintParameters()      = 0;
};

class vbNewForceParameters : public vbForceParameters<ML_CalcType> {
 public:
  vbNewForceParameters() {SetOLShmaddr(overlap);}  // SharedModel Reste

  Parameter *P;

  virtual int GetSize(void) {return 0;}

  virtual ML_CalcType *GetBase(void) {return NULL;}

 private:
  virtual int GetNumParameters() {return -1;}

  ML_CalcType overlap[OVERLAP_LEN];  // Feld 0 ist overlapID, die folgenden sind olparams
};


template<class T>
class vbForceModel : public CellModelMembers<T> {
 public:
  vbForceModel(void) {}

  virtual ~vbForceModel(void) {}

  virtual T Calc(double tinc, T stretch, T velocity, T &Ca, int euler = 1) = 0;

  virtual T CalcTrop(double tinc, T stretch, T velocity, T TCa, int euler = 1) {return 0;}

  virtual void steadyState(T Ca, T stretch) {
    T time, tinc = 0.0001;
    double force;

    for (time = 0.0; time < 1.5; time += tinc) {
      T x = Ca;
      force = Calc(tinc, stretch, 0, x);
      if (std::isnan(force))
        throw kaBaseException("Force is NaN");
    }

    cout<<Ca<<' '<<force<<endl;
  }

  virtual void GetParameterNames(vector<string> &) = 0;
  virtual void Print(ostream &tempstr)             = 0;

  inline void Print(ostream &tempstr, double t, T Ca, T Force) {
    tempstr<<t<<' '<<Ca<<' '<<Force<<' ';
    Print(tempstr);
  }

  inline void Print(ostream &tempstr, double t, T Vm, T Ca, T Force) {
    tempstr<<t<<' '<<Vm<<' '<<Ca<<' '<<Force<<' ';
    Print(tempstr);
  }

  virtual void Init() = 0;

  virtual inline bool InIsTCa(void) {return false;}

  virtual inline T getCa50(void) {return 0;}

  virtual inline T Overlap(T stretch, int overlapID, T olparams[OVERLAP_LEN-1]) {
    T sp;

    switch (overlapID) {
      case 0:

        // Geradenglweichung y = ax + b.
        return stretch < 1.0 ? (1.0-1.4666667*(1.0-stretch)) : (stretch <= 1.1 ? 1.0 : (1.0-1.4666667*(stretch-1.1)));

      case 1:

        // Annaehrung durch 2 Geradenstuecke.
        return stretch <= olparams[2] ? (olparams[0]*stretch+olparams[1]) : (olparams[3]*stretch+olparams[4]);

      case 2:

        // Annaehrung durch 3 Geradenstuecke.
        return stretch <=
               olparams[2] ? (olparams[0]*stretch+olparams[1]) : (stretch <=
                                                                  olparams[5] ? (olparams[3]*stretch+
                                                                                 olparams[4]) : (olparams[6]*stretch+
                                                                                                 olparams[7]));

      case 3:

        // Annaehrung durch 4 Geradenstuecke.
        return stretch <=
               olparams[2] ? (olparams[0]*stretch+olparams[1]) : (stretch <=
                                                                  olparams[5] ? (olparams[3]*stretch+
                                                                                 olparams[4]) : (stretch <=
                                                                                                 olparams
                                                                                                 [8] ? (olparams[6]*
                                                                                                        stretch+
                                                                                                        olparams[7]) : (
                                                                                                   olparams[9]*stretch+
                                                                                                   olparams[10])));

      case 4:

        // Interpolation 1.Ordnung. y = ax + b
        return olparams[0]*stretch+olparams[1];

      case 5:

        // Interpolation 2.Ordnung. y = ax^2 + bx + c
        return olparams[0]*stretch*stretch+olparams[1]*stretch+olparams[2];

      case 6:

        // Interpolation 3.Ordnung. y = ax^3 + bx^2 + cx + d
        sp = stretch*stretch;
        return olparams[0]*sp*stretch+olparams[1]*sp+olparams[2]*stretch+olparams[3];

      case 7:

        // Interpolation 4.Ordnung. y = ax^4 + bx^3 + cx^2 + dx + e
        sp = stretch*stretch;
        return olparams[0]*sp*sp+olparams[1]*sp*stretch+olparams[2]*sp+olparams[3]*stretch+olparams[4];

      case 8:

        // original overlap function Winslow2
        return stretch <= 1.15 ? 1.0-1.4*(1.15-stretch) : 1.0;

      default:
        return 1.;
    }  // switch
  }  // Overlap

  inline void WriteStatus(FILE *writeto) {
    kaWrite(this->GetBase(), sizeof(T), this->GetSize()/sizeof(T), writeto);
  }

  inline void ReadStatus(FILE *readfrom) {
    kaRead(this->GetBase(), sizeof(T), this->GetSize()/sizeof(T), readfrom);
  }
};  // class vbForceModel
#endif  // ifndef FORCE_MODEL_BASIS_H
