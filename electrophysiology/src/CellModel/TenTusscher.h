/*
 * File: TenTusscher.h
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


#ifndef TENTUSSCHERETAL
#define TENTUSSCHERETAL

#include <TenTusscherEtAlParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_TenTusscherEtAlParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_TenTusscherEtAlParameters::a].value
#endif  // ifdef HETERO

class TenTusscherEtAl : public vbElphyModel<ML_CalcType> {
 public:
  TenTusscherEtAlParameters *ptTeaP;

  // state variables (voltage & time)
  ML_CalcType Ca_i, CaSR, Na_i, K_i;
  ML_CalcType m, h, j;      // I_Na
  ML_CalcType xr1, xr2, xs;  // IKr & IKs
  ML_CalcType r, s;         // Ito1
  ML_CalcType d, f, fCa;    // ICa
  ML_CalcType g;            // IRel
  TenTusscherEtAl(TenTusscherEtAlParameters *pp);
  ~TenTusscherEtAl();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return 2.064403e-13*5.6;}

  virtual inline  ML_CalcType GetVm() {return v(VT_V_init);}

  virtual inline  ML_CalcType GetCai() {return Ca_i;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o);}

  virtual inline  ML_CalcType GetNai() {return Na_i;}

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o);}

  virtual inline  ML_CalcType GetKi() {return K_i;}

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o);}

  virtual inline  ML_CalcType GetIto() {return 0.0;}

  virtual inline  ML_CalcType GetIKr() {return 0.0;}

  virtual inline  ML_CalcType GetIKs() {return 0.0;}

  virtual inline int GetSize(void);

  virtual inline ML_CalcType *GetBase(void) {return &Ca_i;}

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp); /*52.*/}

  virtual inline  ML_CalcType GetStimTime() {return 0.001;}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch                                  = 1., int euler = 2);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class TenTusscherEtAl
#endif  // ifndef TENTUSSCHERETAL
