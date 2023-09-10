/* File: MaleckarEtAl.h
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#ifndef MALECKARETAL
#define MALECKARETAL

#include <MaleckarEtAlParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_MaleckarEtAlParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_MaleckarEtAlParameters::a].value
#endif  // ifdef HETERO

class MaleckarEtAl : public vbElphyModel<ML_CalcType> {
 public:
  MaleckarEtAlParameters *ptTeaP;
  ML_CalcType Ca_i;
  ML_CalcType Na_c;
  ML_CalcType Na_i;
  ML_CalcType m;
  ML_CalcType h1;
  ML_CalcType h2;
  ML_CalcType Ca_d;
  ML_CalcType d_L;
  ML_CalcType f_L1;
  ML_CalcType f_L2;
  ML_CalcType K_c;
  ML_CalcType K_i;
  ML_CalcType r;
  ML_CalcType s;
  ML_CalcType a_ur;
  ML_CalcType i_ur;
  ML_CalcType n;
  ML_CalcType pa;
  ML_CalcType Ca_c;
  ML_CalcType O_C;
  ML_CalcType O_TC;
  ML_CalcType O_TMgC;
  ML_CalcType O_TMgMg;
  ML_CalcType O;
  ML_CalcType Ca_rel;
  ML_CalcType Ca_up;
  ML_CalcType O_Calse;
  ML_CalcType F1;
  ML_CalcType F2, yKur;

  MaleckarEtAl(MaleckarEtAlParameters *pp);
  ~MaleckarEtAl();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return (v(VT_Vol_i))*1e-12*(v(VT_V_corrcell));}

  virtual inline  ML_CalcType GetVm() {return v(VT_V_init);}

  virtual inline  ML_CalcType GetCai() {return Ca_i;}

  virtual inline  ML_CalcType GetCao() {return Ca_c;}

  virtual inline  ML_CalcType GetNai() {return Na_i;}

  virtual inline  ML_CalcType GetNao() {return Na_c;}

  virtual inline  ML_CalcType GetKi() {return K_i;}

  virtual inline  ML_CalcType GetKo() {return K_c;}

  virtual inline ML_CalcType *GetBase(void) {return &Ca_i;}

  virtual inline  ML_CalcType GetIto() {return 0.0;}

  virtual inline  ML_CalcType GetIKr() {return 0.0;}

  virtual inline  ML_CalcType GetIKs() {return 0.0;}

  virtual inline int GetSize(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp);}

  virtual inline  ML_CalcType GetStimTime() {return v(VT_stim_duration);}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch                                  = 1., int euler = 2);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class MaleckarEtAl
#endif  // ifndef MALECKARETAL
