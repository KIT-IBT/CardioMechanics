/* File: MitchellSchaeffer.h
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#ifndef MITCHELLSCHAEFFER
#define MITCHELLSCHAEFFER

#include <MitchellSchaefferParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_MitchellSchaefferParameters::a)
#else // ifdef HETERO
# define v(a) ptTeaP->P[NS_MitchellSchaefferParameters::a].value
#endif // ifdef HETERO

class MitchellSchaeffer : public vbElphyModel<ML_CalcType>{
 public:
  MitchellSchaefferParameters *ptTeaP;
  ML_CalcType u, h;

  MitchellSchaeffer(MitchellSchaefferParameters *pp);
  ~MitchellSchaeffer();

#ifdef HETERO
  ParameterSwitch *PS;
#endif // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return 1.15606568e-12;}

  virtual inline  ML_CalcType GetVm() {return v(VT_Vm_init); }

  virtual inline  ML_CalcType GetCai() {return 0.0;}

  virtual inline  ML_CalcType GetCao() {return 0.0;}

  virtual inline  ML_CalcType GetNai() {return 0.0;}

  virtual inline  ML_CalcType GetNao() {return 0.0;}

  virtual inline  ML_CalcType GetKi() {return 0.0;}

  virtual inline  ML_CalcType GetKo() {return 0.0;}

  virtual inline ML_CalcType *GetBase(void) {return &u;}

  virtual inline  ML_CalcType GetIto() {return 0.0;}

  virtual inline  ML_CalcType GetIKr() {return 0.0;}

  virtual inline  ML_CalcType GetIKs() {return 0.0;}

  virtual inline int GetSize(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return 0.0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp);}

  virtual inline  ML_CalcType GetStimTime() {return 0.002;}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void Init();
  virtual  ML_CalcType Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
}; // class MitchellSchaeffer
#endif // ifndef MITCHELLSCHAEFFER
