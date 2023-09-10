/**@file TenTusscher2.h
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

#ifndef TENTUSSCHER2
#define TENTUSSCHER2

#include <TenTusscher2Parameters.h>
//#include <HandleManager.h> idea //lh326



#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_TenTusscher2Parameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_TenTusscher2Parameters::a].value
#endif  // ifdef HETERO

#ifdef KO
        # define vK_o K_o
#else  // ifdef KO
        # define vK_o v(VT_K_o)
#endif  // ifdef KO
#ifdef GCAL
        # define vg_CaL g_CaL
#else  // ifdef GCAL
        # define vg_CaL v(VT_g_CaL)
#endif  // ifdef GCAL
#ifdef GNA
        # define vg_Na g_Na
        # define vkNaCa kNaCa // phase 1b
#else  // ifdef GNA
        # define vg_Na v(VT_g_Na)
        # define vkNaCa v(VT_kNaCa)
#endif  // ifdef GNA
#ifdef ATP
        # define vatpi atpi
        # define vadpi adpi
        # define vknak knak // phase 1b
        # define vVmaxup Vmaxup // phase 1b
        # define vVrel Vrel // phase 1b
#else  // ifdef ATP
        # define vatpi v(VT_atpi)
        # define vadpi v(VT_adpi)
        # define vknak v(VT_knak)
        # define vVmaxup v(VT_Vmaxup)
        # define vVrel v(VT_Vrel)
#endif  // ifdef ATP
#ifdef MGI
        # define vMgi Mgi
#else  // ifdef MGI
        # define vMgi v(VT_Mgi)
#endif  // ifdef MGI

#define pa(value) ((ptTeaP->value[Vi]));

// **********************************************************************************
// *                                                            Ischemia
//                                                                                *
// **********************************************************************************
// the rest moved to TenTusscher2Parameters.h for performance and memory reasons

class TenTusscherEtAl2 : public vbElphyModel<ML_CalcType> {
 public:
  TenTusscher2Parameters *ptTeaP;
  ML_CalcType Ca_i, CaSR, CaSS, Na_i, K_i;
ML_CalcType Array[18]; //lh326 for longPrint arra

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA

  // things that should be included into the backup for ischemia
  ML_CalcType time_bp;
 y
    
  // The following variables are normally of type bool however, the size of bool differs from 32 to 64 bit.
  // Therefore, ML_CalcType has been chosen which is equal in 32 and 64 bit
  // for details see: http://developer.apple.com/macosx/64bit.html
  // ML_CalcType printed0_bp, printed1_bp, printed2_bp;

#  ifdef KO
  ML_CalcType K_o;
#  endif  // ifdef KO
#  ifdef dVmNa
  ML_CalcType dVm_Na;
#  endif  // ifdef dVmNa
#  ifdef GCAL
  ML_CalcType g_CaL;
#  endif  // ifdef GCAL
#  ifdef GNA
  ML_CalcType g_Na, kNaCa;
#  endif  // ifdef GNA
#  ifdef MGI
  ML_CalcType Mgi;
#  endif  // ifdef MGI
#  ifdef ATP
  ML_CalcType atpi, adpi, knak, Vmaxup, Vrel;
#  endif  // ifdef ATP

# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

#ifdef MARKOV_I_NA
  ML_CalcType MINALC3, MINALC2, MINALC1, MINALO, MINAUC3, MINAUC2, MINAUC1, MINAUO, MINAUIC3, MINAUIC2, MINAUIF,
              MINAUIM1, MINAUIM2;
#endif  // ifdef MARKOV_I_NA

  ML_CalcType m, h, j;       // I_Na
  ML_CalcType xr1, xr2, xs;  // IKr & IKs
  ML_CalcType r, s;          // Ito1
  ML_CalcType d, f, f2, fCa;  // ICa
  ML_CalcType Rq;            // new

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA
  ML_CalcType time, stage0, stage1, stage2;
  bool zonefactor_calculated, restore;

#  ifdef ISCHEMIA_VERBOSE_OUPUT
  bool printed0, printed1, printed2;
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT

  // **************************************************************************************************************
  // *** definition of the individual zoneFactors for the ischemia effects. These MUST be stored here in the cell model
  // *** instead of the parameter set because the VT_ZoneFactor will be overwritten for the individual cells by a
  // *** heterogeneous lattice (AddHeteroValue call) and thus, the value of zonefactor_X changes from cell 2 cell

#  ifdef ACIDOSIS
  ML_CalcType zonefactor_fpH;
#  endif  // ifdef ACIDOSIS
#  ifdef HYPERKALEMIA
  ML_CalcType zonefactor_Ko;
#  endif  // ifdef HYPERKALEMIA
#  ifdef HYPOXIA
  ML_CalcType zonefactor_pO;
#  endif  // ifdef HYPOXIA

  // *** end of: definition of the individual zoneFactors for the ischemia effects.
  // **************************************************************************************************************


  // **************************************************************************************************************
  // *** definition of the change of the ischemia effects over time from stage0 via stage1 to stage2
  // *** this values also depend on the individual zonefactor_X values and must therefore also be
  // *** stored here in the cell model

#  ifdef KO
  ML_CalcType KoAdder1, KoAdder2;
#  endif  // ifdef KO

#  ifdef GCAL
  ML_CalcType CaLAdder1, CaLAdder2;
#  endif  // ifdef GCAL

#  ifdef GNA
  ML_CalcType NaAdder1, NaAdder2, kNaCaAdder;
#  endif  // ifdef GNA

#  ifdef dVmNa
  ML_CalcType dVmNaAdder1, dVmNaAdder2;
#  endif  // ifdef dVmNa

#  ifdef ATP
  ML_CalcType ATPAdder1, ATPAdder2, ADPAdder1, ADPAdder2, knakAdder, VmaxupAdder, VrelAdder;
#  endif  // ifdef ATP

#  ifdef MGI
  ML_CalcType MgiAdder1, MgiAdder2;
#  endif  // ifdef MGI

  // *** end of: definition of the change of the ischemia effects over time from stage0 via stage1 to stage2
  // **************************************************************************************************************
# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO

  TenTusscherEtAl2(TenTusscher2Parameters *pp);
  ~TenTusscherEtAl2();
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return v(VT_Vcell);}

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

  virtual inline int          GetSize(void);
  virtual inline ML_CalcType *GetBase(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp);}

  virtual ML_CalcType GetStimTime() {return v(VT_stim_duration);}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA
  virtual void InitIschemiaTimeCourse();
  virtual void CalcIschemiaTimeCourse(double mytinc);
# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
#ifdef MARKOV_I_NA
  virtual void CalcMarkovINa(int, ML_CalcType);
#endif  // ifdef MARKOV_I_NA
  virtual  ML_CalcType Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler);
  virtual void         Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void         LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void         GetParameterNames(vector<string> &getpara);
  virtual void         GetLongParameterNames(vector<string> &getpara);
};  // class TenTusscherEtAl2
#endif  // ifndef TENTUSSCHER2
