/**@file FabbriEtAl.h
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

#ifndef Fabbri_H
#define Fabbri_H

#include <ParameterLoader.h>
#include <FabbriParameters.h>

#undef v

#define v(a) pCmP->P[NS_FabbriParameters::a].value


class Fabbri : public vbElphyModel<ML_CalcType>{
 public:
  FabbriParameters  *pCmP;
  ML_CalcType Na_i;
  ML_CalcType K_i, Ca_i, Mg_i;
  ML_CalcType Ca_jsr, Ca_nsr, Ca_sub;
  ML_CalcType y, m, m_mut, h, h_mut, dL, fL;
  ML_CalcType x;
  ML_CalcType fCa, dT, fT, r_Kur, s;
  ML_CalcType r_to, q, paS, paF, piy, n, a;
  ML_CalcType R_Ca_rel, O_Ca_rel, I_Ca_rel, RI_Ca_rel;
  ML_CalcType fTMM, fCMi, fCMs;
  ML_CalcType fTC, fTMC, fCQ, fBAPTA, fBAPTA_sub;

  Fabbri(FabbriParameters *pp);
  ~Fabbri();

  virtual inline  ML_CalcType Volume() {return 3.14159265359*1e-6*v(VT_R_cell)*1e-6*v(VT_R_cell)*1e-6*v(VT_L_cell);}

  virtual inline  ML_CalcType GetAmplitude() {return 0.0;}

  virtual inline  ML_CalcType GetStimTime() {return 0.003;}

  virtual inline  ML_CalcType GetVm() {return v(VT_Init_Vm); }

  virtual inline  ML_CalcType GetCai() {return Ca_i;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Cao); }

  virtual inline  ML_CalcType GetNao() {return v(VT_Nao); }

  virtual inline  ML_CalcType GetKo() {return v(VT_Ko); }

  virtual inline  ML_CalcType GetNai() {return Na_i;}

  virtual inline  ML_CalcType GetKi() {return K_i;}

  virtual inline int GetSize(void);

  virtual inline ML_CalcType *GetBase(void) {return &Na_i;}

  virtual inline void SetCai(ML_CalcType val) {Ca_i = val;}

  virtual inline void SetNai(ML_CalcType val) {Na_i = val;}

  virtual inline void SetKi(ML_CalcType val) {K_i = val;}

  virtual void Init();
  virtual  ML_CalcType Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler);
  virtual void Print(ostream &tempstr, double t,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double t,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
}; // class Fabbri
#endif // ifndef Fabbri_H
