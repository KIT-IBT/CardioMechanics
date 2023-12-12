/*
 * File: Himeno.h
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


#ifndef HIMENO
#define HIMENO

#include <HimenoParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_HimenoParameters::a)
#else  // ifdef HETERO
# define v(a) pHimP->P[NS_HimenoParameters::a].value
#endif  // ifdef HETERO


class Himeno : public vbElphyModel<ML_CalcType> {
 public:
  HimenoParameters *pHimP;
  ML_CalcType CaMCa;
  ML_CalcType TnChCa, SRCa;
  ML_CalcType p_O_NaT, p_I_2_NaT, p_I_s_NaT, p_O_NaL, p_I_1_NaL, p_I_2_NaL, p_I_s_NaL;
  ML_CalcType chi_r_fast, chi_r_slow;
  ML_CalcType para_Xs1, para_Xs2, i_fast, i_slow;
  ML_CalcType P_7, P_8_13, P_1_6;
  ML_CalcType p_E_1_NCX_blk, p_I_1_NCX_blk, p_I_2_NCX_blk, p_E_1_NCX_iz, p_I_1_NCX_iz, p_I_2_NCX_iz;
  ML_CalcType Y_ooo, Y_ooc, Y_coo, Y_coc, Y_cco, Y_oco, Y_occ;
  ML_CalcType Y_co_iz, Y_oo_iz, Y_oc_iz, Y_co_blk, Y_oo_blk, Y_oc_blk;
  ML_CalcType Ca_2_tot_jnc, Ca_2_tot_iz, Ca_2_tot_blk, Ca_2_tot_SRrl;
  ML_CalcType Ca_2_jnc, Ca_2_iz, Ca_2_blk, Ca_2_SRup, Ca_2_SRrl;
  ML_CalcType Ca_2_nd_00, Ca_2_nd_L0, Ca_2_nd_0R, Ca_2_nd_LR;
  ML_CalcType Nai, Ki;
  ML_CalcType TSCa_3, TSCa_3W, TSCa_3S, TS_S, TS_W;
  ML_CalcType hw, hp, Pb_spm, a;
  ML_CalcType L_bound_iz, H_bound_iz, L_free_jnc;
  ML_CalcType H_free_jnc;


  Himeno(HimenoParameters *pp);
  ~Himeno();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return v(VT_vol)*1e-18;}

  virtual inline  ML_CalcType GetVm() {return v(VT_V_init);}

  virtual inline  ML_CalcType GetCai() {return Ca_2_tot_blk;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Cao);}

  virtual inline  ML_CalcType GetNai() {return Nai;}

  virtual inline  ML_CalcType GetNao() {return v(VT_Nao);}

  virtual inline  ML_CalcType GetKi() {return Ki;}

  virtual inline  ML_CalcType GetKo() {return v(VT_Ko);}

  virtual inline  ML_CalcType GetIto() {return 0.0;}

  virtual inline  ML_CalcType GetIKr() {return 0.0;}

  virtual inline  ML_CalcType GetIKs() {return 0.0;}

  virtual inline int          GetSize(void);
  virtual inline ML_CalcType *GetBase(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp);}

  virtual ML_CalcType GetStimTime() {return v(VT_stim_duration);}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                    int euler);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class Himeno
#endif  // ifndef HIMENO
