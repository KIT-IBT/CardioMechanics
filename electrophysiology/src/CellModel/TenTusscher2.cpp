/*      File: TenTusscher2.cpp
    automatically created by ExtractParameterClass.pl - done by dw (19.07.2007)
    Institute of Biomedical Engineering, Universität Karlsruhe (TH)
    send comments to dw@ibt.uka.de      */

#include <TenTusscher2.h>

#define useCaTab

// #define ACTIVATE_ISAC_CHANNEL
// #define SAC_KUIJPERS
// #define SAC_SACHS

TenTusscherEtAl2::TenTusscherEtAl2(TenTusscher2Parameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_TenTusscher2Parameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

TenTusscherEtAl2::~TenTusscherEtAl2() {}

#ifdef HETERO

inline bool TenTusscherEtAl2::AddHeteroValue(string desc, double val) {
  // cerr << "desc = " << desc.c_str() << ", " << val << endl;
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool TenTusscherEtAl2::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

/*! \fn int TenTusscherEtAl2::GetSize(void)
 *  \brief Return the size of the memory that has to be included into the backup.
 *  \return Size of memory to be used for backup.
 *  This function returns the size of memory to use for the backup snapshot. Be sure to exclude all variables that are
 * not to be backuped.
 *  These are especially the variables used for modeling ischemia and the InitTableDone variable. If this was not done a
 * cell model inited
 *  with pcafile option and a memory snapshot would not be able to use values from an ev-file or through the init
 * function.
 */
inline int TenTusscherEtAl2::GetSize(void) {
  return (&Rq - &Ca_i + 1) * sizeof(ML_CalcType);
}

inline ML_CalcType *TenTusscherEtAl2::GetBase(void) {
  return &Ca_i;
}

inline unsigned char TenTusscherEtAl2::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

#ifdef MARKOV_I_NA

inline void TenTusscherEtAl2::CalcMarkovINa(int Vi, ML_CalcType tinc_int) {
  const ML_CalcType a11_i        = ptTeaP->a11[Vi];
  const ML_CalcType a12_i        = ptTeaP->a12[Vi];
  const ML_CalcType a13_i        = ptTeaP->a13[Vi];
  const ML_CalcType b11_i        = ptTeaP->b11[Vi];
  const ML_CalcType b12_i        = ptTeaP->b12[Vi];
  const ML_CalcType b13_i        = ptTeaP->b13[Vi];
  const ML_CalcType a3_i         = ptTeaP->a3[Vi];
  const ML_CalcType b3_i         = ptTeaP->b3[Vi];
  const ML_CalcType alpha6       = v(VT_alpha6);
  const ML_CalcType beta6        = v(VT_beta6);
  const ML_CalcType transition1  = MINAUIC3*a11_i-MINAUIC2*b11_i;
  const ML_CalcType transition2  = MINAUIC2*a12_i-MINAUIF*b12_i;
  const ML_CalcType transition3  = MINAUIF*ptTeaP->a4[Vi]-MINAUIM1*ptTeaP->b4[Vi];
  const ML_CalcType transition4  = MINAUIM1*ptTeaP->a5[Vi]-MINAUIM2*ptTeaP->b5[Vi];
  const ML_CalcType transition5  = MINAUIC3*a3_i-MINAUC3*b3_i;
  const ML_CalcType transition6  = MINAUIC2*a3_i-MINAUC2*b3_i;
  const ML_CalcType transition7  = MINAUIF*a3_i-MINAUC1*b3_i;
  const ML_CalcType transition8  = MINAUO*ptTeaP->a2[Vi]-MINAUIF*ptTeaP->b2[Vi];
  const ML_CalcType transition9  = MINAUC3*a11_i-MINAUC2*b11_i;
  const ML_CalcType transition10 = MINAUC2*a12_i-MINAUC1*b12_i;
  const ML_CalcType transition11 = MINAUC1*a13_i-MINAUO*b13_i;
  const ML_CalcType transition12 = MINAUC3*alpha6-MINALC3*beta6;
  const ML_CalcType transition13 = MINAUC2*alpha6-MINALC2*beta6;
  const ML_CalcType transition14 = MINAUC1*alpha6-MINALC1*beta6;
  const ML_CalcType transition15 = MINAUO*alpha6-MINALO*beta6;
  const ML_CalcType transition16 = MINALC3*a11_i-MINALC2*b11_i;
  const ML_CalcType transition17 = MINALC2*a12_i-MINALC1*b12_i;
  const ML_CalcType transition18 = MINALC1*a13_i-MINALO*b13_i;

  MINAUIC3 += tinc_int*(-transition1-transition5); checkGatingVariable(MINAUIC3);
  MINAUIC2 += tinc_int*(transition1-transition2-transition6); checkGatingVariable(MINAUIC2);
  MINAUIF  += tinc_int*(transition2-transition3-transition7+transition8); checkGatingVariable(MINAUIF);
  MINAUIM1 += tinc_int*(transition3-transition4); checkGatingVariable(MINAUIM1);
  MINAUIM2 += tinc_int*(transition4); checkGatingVariable(MINAUIM2);
  MINAUC3  += tinc_int*(transition5-transition9-transition12); checkGatingVariable(MINAUC3);
  MINAUC2  += tinc_int*(transition9+transition6-transition10-transition13); checkGatingVariable(MINAUC2);
  MINAUC1  += tinc_int*(transition10+transition7-transition11-transition14); checkGatingVariable(MINAUC1);
  MINAUO   += tinc_int*(-transition8+transition11-transition15); checkGatingVariable(MINAUO);
  MINALC3  += tinc_int*(transition12-transition16); checkGatingVariable(MINALC3);
  MINALC2  += tinc_int*(transition16+transition13-transition17); checkGatingVariable(MINALC2);
  MINALC1  += tinc_int*(transition17+transition14-transition18);  checkGatingVariable(MINALC1);
  MINALO    = 1.0-MINAUIC3-MINAUIC2-MINAUIF-MINAUIM1-MINAUIM2-MINAUC3-MINAUC2-MINAUC1-MINAUO-MINALC3-MINALC2-MINALC1;
  checkGatingVariable(MINALO);
}  // TenTusscherEtAl2::CalcMarkovINa

#endif  // ifdef MARKOV_I_NA

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA

void TenTusscherEtAl2::InitIschemiaTimeCourse() {
  // ------ variables ------
  time = 0;
#  ifdef ISCHEMIA_VERBOSE_OUPUT
  printed0 = 0;
  printed1 = 0;
  printed2 = 0;
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
  stage0 = v(VT_IschemiaStart);
  stage1 = v(VT_IschemiaStage1);
  stage2 = v(VT_IschemiaStage2);

  // Recalculate ischemia stage of this cell according to their position in the tissue
  ML_CalcType diff_a = stage1 - stage0;
  ML_CalcType diff_b = stage2 - stage1;

  // change the 2 into 5 if diffusion takes too long
  // implementation done by Manuel Ifland and documented in his diploma thesis

  ////stage0 = stage0 + stage0 * ( ( 1 - v( VT_DiffusionFactor ) ) );
  ////stage1 = stage1 + diff_a * ( ( 1 - v( VT_DiffusionFactor ) ) );
  ////stage2 = stage2 + diff_b * ( ( 1 - v( VT_DiffusionFactor ) ) );

  // change of the DiffusionFactor implementation by dw
  // the "new" DiffusionFactor is really a factor that is used to compute the "local" stages from the global ones
  // so the value is easier to understand:
  // DiffusionFactor = 1: stage_local = stage_global
  // DiffusionFactor = 2: stage_local = 2 x stage_global
  // -> the larger the DiffusionFactor the larger take the ischemia effects to develop
  // SubEndo: DiffusionFactor = 1
  // SubEpi (in proximity to blood supply): DiffusionFactor = 2
  // however, for the IEEE publication the diffusionFactor was always 1. The sense and meaning of it needs to be
  // discussed. But it is already implemented ...

  stage0 = stage0 * v(VT_DiffusionFactor);
  stage1 = stage1 * v(VT_DiffusionFactor);
  stage2 = stage2 * v(VT_DiffusionFactor);

  // Ko_stage0 = v(VT_K_o);
  // Ko_stage1 = 8.7;
  // Ko_stage2 = 12.5;

#  ifdef dVmNa

  // dVmNa_stage0 = 0;
  // dVmNa_stage1 = 1.7;
  // dVmNa_stage2= 3.4;
  // see shaw97a p 270 - added by dw
#  endif  // ifdef dVmNa

  // gCaL_stage0 = v(VT_g_CaL);
  // gCaL_stage1 = v(VT_g_CaL) * 0.875;
  // gCaL_stage2 = v(VT_g_CaL) * 0.75;

  // gNa_stage0 = v(VT_g_Na);
  // gNa_stage1 = v(VT_g_Na) * 0.875;
  // gNa_stage2 = v(VT_g_Na) * 0.75;

  // Mgi_stage0 = v(VT_Mgi);
  // Mgi_stage1 = 3;
  // Mgi_stage2 = 6;

  // ATP_stage0 = v(VT_atpi);
  // ATP_stage1 = 5.7;
  // ATP_stage2 = 4.6;

  // ADP_stage0 = v(VT_adpi);
  // ADP_stage1 = 57;
  // ADP_stage2 = 99;
  // ------ variables ------

  // KoAdder1 = 0;
  // KoAdder2 = 0;
#  ifdef KO

  // consideration of Ko
  // see p. 1 in Shaw and Rudy: Electrophysiologic Effects of Acute Myocardial Ischemia
  KoAdder1 = (v(VT_K_o_stage1) - v(VT_K_o) ) / (stage1 - stage0);
  KoAdder2 = (v(VT_K_o_stage2) - v(VT_K_o_stage1) ) / (stage2 - stage1);
#  endif  // ifdef KO

#  ifdef dVmNa
  dVmNaAdder1 = (v(VT_dVmNa_stage1) - v(VT_dVmNa_stage0)) / (stage1 - stage0);
  dVmNaAdder2 = (v(VT_dVmNa_stage2) - v(VT_dVmNa_stage1)) / (stage2 - stage1);
#  endif  // ifdef dVmNa

#  ifdef GCAL

  // consideration of Ca and Na
  CaLAdder1 = (v(VT_gCaL_stage1) - v(VT_g_CaL) ) / (stage1 - stage0);
  CaLAdder2 = (v(VT_gCaL_stage2) - v(VT_gCaL_stage1) ) / (stage2 - stage1);
#  endif  // ifdef GCAL

#  ifdef GNA
  NaAdder1 = (v(VT_gNa_stage1) - v(VT_g_Na) ) / (stage1 - stage0);
  NaAdder2 = (v(VT_gNa_stage2) - v(VT_gNa_stage1) ) / (stage2 - stage1);

  kNaCaAdder = (v(VT_kNaCa_1b) - v(VT_kNaCa) ) / (stage2 - stage0);  // phase 1b
#  endif  // ifdef GNA

#  ifdef MGI

  // consideration of Mg2+ (see Carmeliet, p. 963)
  MgiAdder1 = (v(VT_Mgi_stage1) - v(VT_Mgi) ) / (stage1 - stage0);
  MgiAdder2 = (v(VT_Mgi_stage2) - v(VT_Mgi_stage1) ) / (stage2 - stage1);
#  endif  // ifdef MGI

#  ifdef ATP

  // consideration of ATP and ADP
  ATPAdder1 = (v(VT_ATP_stage1) - v(VT_atpi) ) / (stage1 - stage0);
  ATPAdder2 = (v(VT_ATP_stage2) - v(VT_ATP_stage1) ) / (stage2 - stage1);

  ADPAdder1 = (v(VT_ADP_stage1) - v(VT_adpi) ) / (stage1 - stage0);
  ADPAdder2 = (v(VT_ADP_stage2) - v(VT_ADP_stage1) ) / (stage2 - stage1);

  knakAdder   = (v(VT_knak_1b) - v(VT_knak) ) / (stage2 - stage0); // phase 1b
  VmaxupAdder = (v(VT_Vmaxup_1b) - v(VT_Vmaxup) ) / (stage2 - stage0);  // phase 1b
  VrelAdder   = (v(VT_Vrel_1b) - v(VT_Vrel) ) / (stage2 - stage0); // phase 1b
#  endif  // ifdef ATP
}  // TenTusscherEtAl2::InitIschemiaTimeCourse

# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

void TenTusscherEtAl2::Init() {
#if KADEBUG
  cerr << "#initializing Class: TenTusscherEtAl2 ... " << endl;
#endif  // if KADEBUG
  m    = v(VT_m_init);
  h    = v(VT_h_init);
  j    = v(VT_j_init);
  xr1  = v(VT_xr1_init);
  xr2  = v(VT_xr2_init);
  xs   = v(VT_xs_init);
  r    = v(VT_r_init);
  s    = v(VT_s_init);
  d    = v(VT_d_init);
  f    = v(VT_f_init);
  f2   = v(VT_f2_init); // new
  fCa  = v(VT_fCa_init);
  Rq   = v(VT_Rq_init);   // new
  Ca_i = v(VT_Cai_init);
  CaSR = v(VT_CaSR_init);
  CaSS = v(VT_CaSS_init);  // new
  Na_i = v(VT_Nai_init);
  K_i  = v(VT_Ki_init);

#ifdef MARKOV_I_NA
  MINALC3  = v(VT_initMINALC3);
  MINALC2  = v(VT_initMINALC2);
  MINALC1  = v(VT_initMINALC1);
  MINALO   = v(VT_initMINALO);
  MINAUC3  = v(VT_initMINAUC3);
  MINAUC2  = v(VT_initMINAUC2);
  MINAUC1  = v(VT_initMINAUC1);
  MINAUO   = v(VT_initMINAUO);
  MINAUIC3 = v(VT_initMINAUIC3);
  MINAUIC2 = v(VT_initMINAUIC2);
  MINAUIF  = v(VT_initMINAUIF);
  MINAUIM1 = v(VT_initMINAUIM1);
  MINAUIM2 = v(VT_initMINAUIM2);
#endif  // ifdef MARKOV_I_NA

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef MGI
  Mgi = v(VT_Mgi);
# endif  // ifdef MGI
# ifdef ATP
  atpi   = v(VT_atpi);
  adpi   = v(VT_adpi);
  knak   = v(VT_knak); // phase 1b
  Vmaxup = v(VT_Vmaxup);  // phase 1b
  Vrel   = v(VT_Vrel); // phase 1b
# endif  // ifdef ATP
# ifdef ISCHEMIA

  // if restore was wrongly set to one and no backup exists, we prevent errors by intializing the variables
#  ifdef KO
  K_o = v(VT_K_o);
#  endif  // ifdef KO
#  ifdef dVmNa
  dVm_Na = 0;
#  endif  // ifdef dVmNa
#  ifdef GCAL
  g_CaL = v(VT_g_CaL);
#  endif  // ifdef GCAL
#  ifdef GNA
  g_Na  = v(VT_g_Na);
  kNaCa = v(VT_kNaCa);  // phase 1b
#  endif  // ifdef GNA

#  ifdef ISCHEMIA_VERBOSE_OUPUT
#   ifdef HYPERKALEMIA
  cout << "HYPERKALEMIA is enabled ...\n";
#   endif  // ifdef HYPERKALEMIA
#   ifdef ACIDOSIS
  cout << "ACIDOSIS is enabled ...\n";
#   endif  // ifdef ACIDOSIS
#   ifdef HYPOXIA
  cout << "HYPOXIA is enabled ...\n";
#   endif  // ifdef HYPOXIA
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
  time_bp = 0;

  //    printed0_bp = 0;
  //    printed1_bp = 0;
  //    printed2_bp = 0;

  // do we have do restore from backup?
  if (1 == v(VT_RestoreIschemia) ) {
    restore = true;
  } else {
    restore = false;
  }

  // initialize Ischemia
  InitIschemiaTimeCourse();

  // the zonefactor has to be calculated once again in Calc() because it is not (in case of HETERO with hlf file) yet
  // available in Init()
  zonefactor_calculated = false;
# else  // ifdef ISCHEMIA
#  ifdef ISCHEMIA_VERBOSE_OUPUT
  cout << "I_Katp channel ON / Ischemia OFF" << endl;
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
# endif  // ifdef ISCHEMIA

#else  // ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA_VERBOSE_OUPUT
  cout << "I_Katp channel OFF / Ischemia OFF" << endl;
# endif  // ifdef ISCHEMIA_VERBOSE_OUPUT

  // initialize I_Katp to zero because it is deactivated (see ACTIVATE_IKATP_CHANNEL)
  // I_Katp = 0;
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

#ifdef ACTIVATE_ISAC_CHANNEL
  cout << "I_SAC channel ON" << endl;
# ifdef SAC_SACHS
  cout << "I_SAC SACHS" << endl;
# endif // SAC_SACHS
# ifdef SAC_KUIJPERS
  cout << "I_SAC KUIJPERS" << endl;
# endif // SAC_KUIJPERS
#endif //idef ACTIVATE_ISAC_CHANNEL
}  // TenTusscherEtAl2::Init

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA

void TenTusscherEtAl2::CalcIschemiaTimeCourse(double mytinc) {
  if (!zonefactor_calculated) {
    // we have to run init of the ischemia time course again once in calc because now the DiffusionFactor is available
    InitIschemiaTimeCourse();

    /* zones are defined as follows:
       NZ: 0
       BZ: 0 - 1
       CZ: 1
     */

    // ZoneFactor of Ko
    // Ko_begin = 0;
    // Ko_end = 1;
#  ifdef HYPERKALEMIA
    zonefactor_Ko = v(VT_ZoneFactor) - v(VT_Ko_ZoneFactor_Begin);

    if (v(VT_ZoneFactor) >= v(VT_Ko_ZoneFactor_End) ) {
      zonefactor_Ko = v(VT_Ko_ZoneFactor_End) - v(VT_Ko_ZoneFactor_Begin);
    }
    if (v(VT_ZoneFactor) <= v(VT_Ko_ZoneFactor_Begin) ) {
      zonefactor_Ko = 0;
    }
#  endif  // ifdef HYPERKALEMIA

    // cout << "Zonefactor Ko = " << zonefactor_Ko << endl;

#  ifdef dVmNa

    // Zonefactor for dVmNa
    // dVmNa_begin = 0.5;
    // dVmNa_end = 1;

    // zonefactor_dVmNa = v(VT_ZoneFactor) - v(VT_dVmNa_ZoneFactor_Begin);
    // if (v(VT_ZoneFactor) >= v(VT_dVmNa_ZoneFactor_End)) {
    //  zonefactor_dVmNa = v(VT_dVmNa_ZoneFactor_End) - v(VT_dVmNa_ZoneFactor_Begin);
    // }
    // if (v(VT_ZoneFactor) <= v(VT_dVmNa_ZoneFactor_Begin)){
    //  zonefactor_dVmNa = 0;
    // }
    // da dVmNa auch bei Acidosis: zonefactor_dVmNa = zonefactor_fpH
#  endif  // ifdef dVmNa

    // ZoneFactor representing the pH value (channel conductivity of I_Na and I_CaL)
    // fpH_begin = 0.5;
    // fpH_end = 1;

#  ifdef ACIDOSIS
    zonefactor_fpH = v(VT_ZoneFactor) - v(VT_fpH_ZoneFactor_Begin);

    if (v(VT_ZoneFactor) >= v(VT_fpH_ZoneFactor_End) ) {
      zonefactor_fpH = v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin);
    }
    if (v(VT_ZoneFactor) <= v(VT_fpH_ZoneFactor_Begin) ) {
      zonefactor_fpH = 0;
    }
#  endif  // ifdef ACIDOSIS

    // cout << "Zonefactor fpH = " << zonefactor_fpH << endl;

    // ZoneFactor representing the ATP and ADP concentration
    // pO_begin = 0;
    // pO_end = 0.1;
#  ifdef HYPOXIA
    zonefactor_pO = v(VT_ZoneFactor) - v(VT_pO_ZoneFactor_Begin);

    if (v(VT_ZoneFactor) >= v(VT_pO_ZoneFactor_End) ) {
      zonefactor_pO = v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin);
    }
    if (v(VT_ZoneFactor) <= v(VT_pO_ZoneFactor_Begin) ) {
      zonefactor_pO = 0;
    }
#  endif  // ifdef HYPOXIA

    // cout << "Zonefactor pO = " << zonefactor_pO << endl;

#  ifdef ISCHEMIA_VERBOSE_OUPUT
    if (!restore) {  // don't print this message if we are restoring from an ischemia backup
      cout << "I_Katp channel ON / Ischemia ON (S0: " << stage0 << ", S1: " << stage1 << ", S2: " << stage2 << ")" <<
        endl;
    }
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT

    // the ZoneFactor has been calculated, this only has to be done once a time
    zonefactor_calculated = true;
  }

  // load the backup if wanted
  if (restore) {
    time = time_bp;

#  ifdef ISCHEMIA_VERBOSE_OUPUT
    /*
                    if ( 0 == printed0_bp )
                    {
                            printed0 = false;
                    }
                    else
                    {
                            printed0 = true;
                    }

                    if ( 0 == printed1_bp )
                    {
                            printed1 = false;
                    }
                    else
                    {
                            printed1 = true;
                    }

                    if ( 0 == printed2_bp )
                    {
                            printed2 = false;
                    }
                    else
                    {
                            printed2 = true;
                    }
     */
    cout << "Restoring Ischemia (@: " << time << ", S0: " << stage0 << ", S1: " << stage1 << ", S2: " << stage2 <<
      ")" << endl;
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT

    // the backup has been read once
    restore = false;
  }

  // calculate time in milliseconds
  time   += (mytinc * 1000);
  time_bp = time;  // backup

  if (time >= stage2) {
    // nothing happens any more
#  ifdef ISCHEMIA_VERBOSE_OUPUT
    if (!printed2) {
      cout << "Stage2 reached, nothing happens any more, time = " << time << endl;
      printed2 = true;

      // printed2_bp = 1;       // backup
    }
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
  } else if (time >= stage1) {
#  ifdef ISCHEMIA_VERBOSE_OUPUT
    if (!printed1) {
      cout << "Stage1 reached, time = " << time << endl;

      // cout << "DOING NOTHING FROM NOW ON!" << endl;
      printed1 = true;

      // printed1_bp = 1;       // backup
    }
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
#  ifdef KO
    K_o += (KoAdder2 * mytinc * 1000 * zonefactor_Ko / (v(VT_Ko_ZoneFactor_End) - v(VT_Ko_ZoneFactor_Begin) ) );
#  endif  // ifdef KO
#  ifdef dVmNa
    dVm_Na +=
      (dVmNaAdder2 * mytinc * 1000 * zonefactor_fpH / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
#  endif  // ifdef dVmNa
#  ifdef GCAL
    g_CaL += (CaLAdder2 * mytinc * 1000 * zonefactor_fpH  / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
#  endif  // ifdef GCAL
#  ifdef GNA
    g_Na  += (NaAdder2 * mytinc * 1000 * zonefactor_fpH  / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
    kNaCa += (kNaCaAdder * mytinc * 1000 * zonefactor_fpH  / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );  //
                                                                                                                          //
                                                                                                                          // phase
                                                                                                                          //
                                                                                                                          // 1b
#  endif  // ifdef GNA
#  ifdef MGI
    Mgi += (MgiAdder2 * mytinc * 1000);
#  endif  // ifdef MGI
#  ifdef ATP
    atpi += (ATPAdder2 * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );
    adpi += (ADPAdder2 * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );
    knak += (knakAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) ); // phase
                                                                                                                   // 1b
    Vmaxup += (VmaxupAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );  //
                                                                                                                        //
                                                                                                                        // phase
                                                                                                                        //
                                                                                                                        // 1b
    Vrel += (VrelAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) ); // phase
                                                                                                                   // 1b
#  endif  // ifdef ATP

    // atpi += ( ATPAdder2 * mytinc * 1000 * zonefactor_pO );
    // adpi += ( ADPAdder2 * mytinc * 1000 * zonefactor_pO );
  } else if (time >= stage0) {
#  ifdef ISCHEMIA_VERBOSE_OUPUT
    if (!printed0) {
      cout << "Beginning with Ischemia, time = " << time << endl;
      printed0 = true;

      // printed0_bp = 1;       // backup
    }
#  endif  // ifdef ISCHEMIA_VERBOSE_OUPUT
#  ifdef KO
    K_o += (KoAdder1 * mytinc * 1000 * zonefactor_Ko / (v(VT_Ko_ZoneFactor_End) - v(VT_Ko_ZoneFactor_Begin) ) );
#  endif  // ifdef KO
#  ifdef dVmNa
    dVm_Na +=
      (dVmNaAdder1 * mytinc * 1000 * zonefactor_fpH / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
#  endif  // ifdef dVmNa
#  ifdef GCAL
    g_CaL += (CaLAdder1 * mytinc * 1000 * zonefactor_fpH / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
#  endif  // ifdef GCAL
#  ifdef GNA
    g_Na  += (NaAdder1 * mytinc * 1000 * zonefactor_fpH / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );
    kNaCa += (kNaCaAdder * mytinc * 1000 * zonefactor_fpH  / (v(VT_fpH_ZoneFactor_End) - v(VT_fpH_ZoneFactor_Begin) ) );  //
                                                                                                                          //
                                                                                                                          // phase
                                                                                                                          //
                                                                                                                          // 1b
#  endif  // ifdef GNA
#  ifdef MGI
    Mgi += (MgiAdder1 * mytinc * 1000);
#  endif  // ifdef MGI
#  ifdef ATP
    atpi += (ATPAdder1 * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );
    adpi += (ADPAdder1 * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );
    knak += (knakAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) ); // phase
                                                                                                                   // 1b
    Vmaxup += (VmaxupAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) );  //
                                                                                                                        //
                                                                                                                        // phase
                                                                                                                        //
                                                                                                                        // 1b
    Vrel += (VrelAdder * mytinc * 1000 * zonefactor_pO / (v(VT_pO_ZoneFactor_End) - v(VT_pO_ZoneFactor_Begin) ) ); // phase
                                                                                                                   // 1b
#  endif  // ifdef ATP
  }
}  // TenTusscherEtAl2::CalcIschemiaTimeCourse

# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

ML_CalcType TenTusscherEtAl2::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                                   int euler                                            = 2) {
  ML_CalcType svolt = V*1000;  // membrane voltage in mV
  const int   Vi    = (int)(DivisionTab*(RangeTabhalf+svolt)+.5); // array position

  // ECA Tabellisierung

  // convert from nA/cell to pA/pF with Cm: 2 mueF/cm^2=2e-2 F/m^2, S=0.2/mue m=0.2e6/m -fs
  //  i_external=i_external/(2e-2*Volume()*0.2e6*1e9);
  // hard-coded, because S is currently questionable in TenTusscher2.h
  // i_external = i_external/4.59312e-2;  // convert from nA to pA/pF

  // benchmark test: S=0.14/mue = 0.14e6/m , Cm: 1 mueF/cm^2=1e-2 F/m^2, Volume()=1.6404e-14 m^3
  // i_external=i_external/(1e-2*Volume()*0.14e6*1e9);
  // hard-coded for Benchmarkt test 2015
  i_external = i_external/2.29656e-2; // convert from nA to pA/pF

#ifdef ACTIVATE_IKATP_CHANNEL
# ifdef ISCHEMIA
  CalcIschemiaTimeCourse(tinc);
# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

#ifdef useCaTab
  static double CaiTabConst1 = -v(VT_CaiMin)/v(VT_StepCai)+0.5;
  static double CaiTabConst2 = 1./v(VT_StepCai);
  int Caitab                 = (int)(CaiTabConst2*Ca_i+CaiTabConst1);
#endif  // ifdef useCaTab
  const  ML_CalcType K_o_int = vK_o;
  const  ML_CalcType RTONF   = v(VT_RTONF);
  const  ML_CalcType EK      = RTONF*(log((K_o_int/K_i)));    // 25
  const  ML_CalcType VminEK  = svolt-EK;
  const  ML_CalcType ENa     = RTONF*(log((v(VT_Na_o)/Na_i))) // 25
#ifdef dVmNa
    +dVm_Na

    // VminENa = Vm-(ENa+dVmNa) = (Vm - dVmNa) - ENa = VmNa - ENa with VmNa = Vm - dVmNa - see p270 in shaw97a
#endif  // ifdef dVmNa
  ;
  const  ML_CalcType VminENa = svolt-ENa;
  const  ML_CalcType EKs     = RTONF*(log((K_o_int + v(VT_pKNa)*v(VT_Na_o))/(K_i+v(VT_pKNa)*Na_i)));       // 26
  const  ML_CalcType AK1     = 0.1/(1.+exp(0.06*(VminEK-200)));                                            // 77
  const  ML_CalcType BK1     = (3.*exp(0.0002*(VminEK+100))+exp(0.1*(VminEK-10)))/(1.+exp(-0.5*(VminEK))); // 78
  const  ML_CalcType rec_iK1 = AK1/(AK1+BK1);
  const  ML_CalcType I_pCa   = v(VT_g_pCa)*Ca_i/(v(VT_KpCa)+Ca_i);                                         // 82
  const  ML_CalcType I_pK    = ptTeaP->rec_ipK[Vi]*VminEK;                                                 // 83
  const  ML_CalcType I_to    = v(VT_g_to)*r*s*VminEK;                                                      // 55
                                                                                                           // Transient
                                                                                                           // Outward
                                                                                                           // Current
  const  ML_CalcType I_Kr = v(VT_g_Kr)*sqrt(K_o_int/5.4)*xr1*xr2*VminEK;                                // 67 Rapid
                                                                                                        // Delayed
                                                                                                        // Rectifier
                                                                                                        // Current
  const  ML_CalcType I_Ks = v(VT_g_Ks)*xs*xs*(svolt-EKs);                                               // 62 Slow
                                                                                                        // Delayed
                                                                                                        // Rectifier
                                                                                                        // Current
  const  ML_CalcType I_K  = I_Kr+I_Ks;
  const  ML_CalcType I_K1 = v(VT_g_K1)*sqrt(K_o_int/5.4)*rec_iK1*VminEK;                                // 76 Inward
                                                                                                        // Rectifier
                                                                                                        // K+ Current
  const  ML_CalcType I_bNa = v(VT_g_bNa)*VminENa;                                                        // 84
                                                                                                         // Background
                                                                                                         // Currents
#ifdef useCaTab
  const  ML_CalcType I_bCa = v(VT_g_bCa)*(svolt-(ptTeaP->ECA[Caitab]));
#else  // ifdef useCaTab
  const  ML_CalcType ECa   = 0.5*RTONF*(log((v(VT_Ca_o)/Ca_i))); // 25
  const  ML_CalcType I_bCa = v(VT_g_bCa)*(svolt-ECa);            // 85 Background Currents
#endif  // ifdef useCaTab

#ifdef MARKOV_I_NA
  CalcMarkovINa(Vi, tinc*1000);
  const  ML_CalcType I_Na = 5.0*vg_Na*(MINAUO+MINALO)*VminENa;
#else  // ifdef MARKOV_I_NA
  const  ML_CalcType I_Na = vg_Na*m*m*m*h*j*VminENa;  // 27 Fast Na+ Current
#endif  // ifdef MARKOV_I_NA

  const  ML_CalcType I_CaL  = vg_CaL * d*f*f2*fCa * (((ptTeaP->CaL_P1[Vi])*CaSS)-(ptTeaP->CaL_P2[Vi]));
  const  ML_CalcType I_NaCa = vkNaCa*(ptTeaP->NaCa_P1[Vi])*Na_i*Na_i*Na_i-vkNaCa*(ptTeaP->NaCa_P2[Vi])*Ca_i;
  const  ML_CalcType I_NaK  = vknak*(K_o_int/(K_o_int+v(VT_KmK)))*(Na_i/(Na_i+v(VT_KmNa)))*(ptTeaP->NaK_P1[Vi]);

  /*!< ------------------------------------------- ATP channel ---------------------------------------------------- */
  /*!< original rudy channel with modifications by Ferrero */

#ifdef ACTIVATE_IKATP_CHANNEL

  /// Maximum conductance of the ATP-sensitive K channel (mS/uF)
  /// - Nichols et al. determined a maximum conductance of 195 nS/cell in whole guinea-pig recordings at Ko = 4.0 mM
  ///   and physiologic concentrations of ADP, GDP and free Mg2+ (see p. 282 in Nichols et al)
  /// - assuming that the guinea pig ventricular cell is approximated by cylinder length 100 µm and circumference of 50
  /// µm, then the
  ///   surface area would be 5000 µm^2 (= 5 * 10^-5 cm^2), suggesting that there are approx. 25000 channels per cell
  /// (Nichols paper p. 286, footnote)
  // has been integrated into gamma: const double gkatp = P[VT_gkatp].value / P[VT_nicholsarea].value;

  /// this defines the relation between the single channel conduction and the K concentration in mmol/L
  /// the maximum conductance (gkatp) used here is different than those used by Ferrero et al. and measured by M. Kakei
  /// et al.
  /// Kakei defines: gamma = 23.6 * Ko^0.24 = 35.375 * (Ko / 5.4)^0.24 (see pages 448 in Kakei and 18 in Ferrero)
  const ML_CalcType gamma = v(VT_gammaconst) * pow(K_o_int * 0.25, 0.24);

  // intracellular Mg2+
  // KhMg has been integrated into f_M
  // KhMg[Vi] = ( ( 0.65 / sqrt( P[VT_Ko].value + 5 ) ) * exp( - 2 * 0.32 * P[VT_FdRT].value * V ) );
  /// (E11)
  const ML_CalcType f_M = 1 / (1 + vMgi / ( (0.65 / sqrt(K_o_int + 5) ) * exp(-2 * 0.32 * v(VT_inverseRTONF) * V) ) );

  // intracellular Na+ ions
  /// (E14)
  const ML_CalcType f_N = 1 / (1 + ( (Na_i / ptTeaP->KhNa[Vi]) * (Na_i / ptTeaP->KhNa[Vi]) ) );

  // Fraction of activated K_ATP channels
  /// maximum-inhibition constant (E18)
  // mw104: changed constants in calculation of K_m (original: K_m = ( 35.8 + 17.9 * pow( vadpi, 0.256 ) ) *
  // v(VT_Km_factor);)
  // this change is necessary to obtain the same APs in healthy cells as in the original ten Tusscher 2006 model
  // for this purpose, the factor K_m was adjusted to almost zero (0.001) in healthy cells (ADPi = 15 µmol), so that
  // i_Katp almost vanishes
  // without this correction, the APs are much shorter than in the original model and the plateau phase of epicardial
  // cells has a much lower amplitude than that of endo and mid, which leads to ST segment shifts in the ECG!!!
  // Corresponding values:      Stage 1 (orig. 57 µmol) --> 87.4593
  //                                                                                            Phase 1b (orig. 110
  // µmol) --> 101.5289
  const ML_CalcType K_m = (-151.0919 + 75.5379 * pow(vadpi, 0.256) ) * v(VT_Km_factor);

  /// Hill coefficient (E19)
  const ML_CalcType H = 1.3 + 0.74 * exp(-0.09 * vadpi);

  /// (E17) (see Kakei, p. 454, 456)
  /// convert K_m from µmol to mmol (by dividing by 1000, or multiplying by 0.001) because K_m calculation uses ADPi
  /// which is specified in µmol
  const ML_CalcType f_ATP = 1 / (1 + pow(vatpi / (K_m * 0.001), H) );

  /// anoxia formulation for I_CaL
  // I_CaL *= 1 / ( 1 + pow( ( 1.4 / atpi ), 2.6 ) );

  /// p_0 = 0.91                        open state probability (Kakei, p. 456)
  /// current
  const ML_CalcType I_Katp = gamma * 0.91 * f_ATP * f_M * f_N * v(VT_f_T) * VminEK;
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

  /*!< ------------------------------------------- ATP channel ---------------------------------------------------- */


  /*!< ------------------------------------------- SAC channel ---------------------------------------------------- */
#ifdef ACTIVATE_ISAC_CHANNEL
# ifdef SAC_SACHS
  const ML_CalcType I_SAC = v(VT_g_SAC) * (svolt - v(VT_ESAC) ) /
    (1.0 + v(VT_KSAC) * exp(-v(VT_alphaSAC) * (stretch-1.0) ) );
# endif //ifdef SAC_SACHS
# ifdef SAC_KUIJPERS
  const double gsac    = v(VT_Gsac) / (1.0 + v(VT_Ksac) * exp(-v(VT_alphasac) * (stretch - 1.0) ) );
  double VinverseRTONF = V*v(VT_inverseRTONF);
  const double Csac1   = gsac * VinverseRTONF * v(VT_F)/(1.0 - exp(-VinverseRTONF) );
  const double Csac2   = gsac *4.0*VinverseRTONF*v(VT_F)/(1.0 - exp(-2.0*VinverseRTONF) );
  const double I_sacNa = v(VT_PsacNa) * Csac1 * (Na_i-v(VT_Na_o)*exp(-v(VT_inverseRTONF)*svolt));
  const double I_sacK  = v(VT_PsacK) * Csac1 * (K_i-v(VT_K_o)*exp(-v(VT_inverseRTONF)*svolt));
  const double I_sacCa = v(VT_PsacCa) * Csac2 * (Ca_i-v(VT_Ca_o)*exp(-2.0*v(VT_inverseRTONF)*svolt));
  const double I_SAC   = I_sacNa + I_sacK + I_sacCa;

# endif  // ifdef SAC_KUIJPERS
#endif  // ifdef ACTIVATE_ISAC_CHANNEL

  /*!< ------------------------------------------- SAC channel ---------------------------------------------------- */

  const ML_CalcType I_tot = I_K+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaK+I_NaCa+I_pCa+I_pK
#ifdef ACTIVATE_IKATP_CHANNEL
    +I_Katp
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
#ifdef ACTIVATE_ISAC_CHANNEL
    +I_SAC
#endif  // ifdef ACTIVATE_ISAC_CHANNEL
    -i_external;

  // update concentrations
  const  ML_CalcType Vss        = v(VT_Vss);
  const  ML_CalcType Cm         = v(VT_C);
  const  ML_CalcType max_sr     = v(VT_max_sr);
  const  ML_CalcType Bufsr      = v(VT_Bufsr);
  const  ML_CalcType EC         = v(VT_EC);
  const  ML_CalcType Kbufsr     = v(VT_Kbufsr);
  const  ML_CalcType Bufss      = v(VT_Bufss);
  const  ML_CalcType Kbufss     = v(VT_Kbufss);
  const  ML_CalcType Kbufc      = v(VT_Kbufc);
  const  ML_CalcType inverseviF = v(VT_inverseviF);
  const  ML_CalcType Caisquare  = Ca_i*Ca_i;
  const  ML_CalcType CaSSsquare = CaSS*CaSS;                                  // new
  const  ML_CalcType kCaSR      = max_sr-((max_sr-v(VT_min_sr))/(1+(EC/CaSR)*(EC/CaSR))); // new n37
  const  ML_CalcType k2         = v(VT_ks2)*kCaSR;                            // new n36
  const  ML_CalcType k1         = v(VT_ks1)/kCaSR;                            // new n35
  ML_CalcType dRq               = 1000.0*((-k2*CaSS*Rq)+(v(VT_k4)*(1.0-Rq))); // new n34
  Rq += tinc*dRq;                                                             // new
  const  ML_CalcType O      = k1*CaSSsquare*Rq/(v(VT_k3)+k1*CaSSsquare);      // new n33
  const  ML_CalcType I_rel  = vVrel*O*(CaSR-CaSS);                            // new n31
  const  ML_CalcType I_leak = 0.00036*(CaSR-Ca_i);                            // new n29
  const  ML_CalcType I_xfer = v(VT_Vxfer)*(CaSS-Ca_i);                        // new n32
  const  ML_CalcType SERCA  = vVmaxup/(1.+(v(VT_Kupsquare)/Caisquare));       // 87 n30 I_up
  const  ML_CalcType dCaSR  = tinc*1000.0*(SERCA-I_rel-I_leak);               // new n41
  const  ML_CalcType CaCSQN = Bufsr*CaSR/(CaSR+Kbufsr);                       // 95
  const  ML_CalcType bjsr   = Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr; // CaSum_1;
  const  ML_CalcType cjsr   = Kbufsr*(CaCSQN+dCaSR+CaSR); // CaSum_1;
  CaSR = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2.;                                    // Lösung der quadr. Gleichung:
                                                                              // CaSR^2+bjsr*CaSR-cjsr=0 //vgl. Zeng
                                                                              // 1995
  const  ML_CalcType CaSSBuf = Bufss*CaSS/(CaSS+Kbufss);                      // new n42
  const  ML_CalcType dCaSS   = tinc*1000.0*
    (-I_xfer*(v(VT_Vc)/Vss)+I_rel*(v(VT_Vsr)/Vss)+(-I_CaL*v(VT_inversevssF2)*Cm));
  const  ML_CalcType bjss = Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;  // CaSum_3;                                //new
  const  ML_CalcType cjss = Kbufss*(CaSSBuf+dCaSS+CaSS);  // CaSum_3;                                            //new
  CaSS = (sqrt(bjss*bjss+4.*cjss)-bjss)/2.;                // new Lösung der quadr. Gleichung: CaSS^2+bjss*CaSS-cjss=0
  const  ML_CalcType CaBuf = v(VT_Bufc)*Ca_i/(Ca_i+Kbufc);  // 93
  const  ML_CalcType dCai  = tinc*1000.0*
    ((-(I_bCa+I_pCa-2*I_NaCa    // +I_CaL missing?? lh326
#ifdef SAC_KUIJPERS
        + I_sacCa
#endif  // ifdef SAC_KUIJPERS
        )*v(VT_inverseviF2)*Cm)-(SERCA-I_leak)*(1.0/v(VT_VcdVsr))+I_xfer);
  const  ML_CalcType bc = v(VT_BufcPKbufc)-CaBuf-dCai-Ca_i;
  const  ML_CalcType cc = Kbufc*(CaBuf+dCai+Ca_i);  // CaSum_2;
  Ca_i = (sqrt(bc*bc+4*cc)-bc)/2.;                                        // Lösung der quadr. Gleichung:
                                                                          // Ca_i^2+bc*Ca_i-cc=0 //94 //vgl. Zeng 1995
  const  ML_CalcType dNai = -(I_Na+I_bNa+3*I_NaK+3*I_NaCa
#ifdef SAC_KUIJPERS
                              + I_sacNa
#endif  // ifdef SAC_KUIJPERS
                              )*inverseviF*Cm;  // 97
  Na_i += tinc*dNai*1000;
  const  ML_CalcType dKi = -(-i_external+I_K1+I_to+I_K-2*I_NaK+I_pK
#ifdef ACTIVATE_IKATP_CHANNEL
                             +I_Katp
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
#ifdef SAC_KUIJPERS
                             + I_sacK
#endif  // ifdef SAC_KUIJPERS
                             )*inverseviF*Cm;  // new obwohl nicht in Paper, aber in TT2 Source Code
  K_i += tinc*dKi*1000.0;
  const  ML_CalcType CaSS005 = (CaSS/0.05);
  const  ML_CalcType FCa_INF = 0.6/(1.0+CaSS005*CaSS005)+0.4;
  const  ML_CalcType taufca  = 0.080/(1.0+CaSS005*CaSS005)+0.002;

  const  ML_CalcType M_INF   = pa(m_inf);
  const  ML_CalcType H_INF   = pa(h_inf);
  const  ML_CalcType J_INF   = pa(j_inf);
  const  ML_CalcType Xr1_INF = pa(Xr1_inf);
  const  ML_CalcType Xr2_INF = pa(Xr2_inf);
  const  ML_CalcType Xs_INF  = pa(Xs_inf);
  const  ML_CalcType S_INF   = pa(s_inf);
  const  ML_CalcType R_INF   = pa(r_inf);
  const  ML_CalcType D_INF   = pa(d_inf);
  const  ML_CalcType F_INF   = pa(f_inf);
  const  ML_CalcType F2_INF  = pa(f2_inf); // new

  m   = M_INF-(M_INF-m)*pa(exptau_m);
  h   = H_INF-(H_INF-h)*pa(exptau_h);
  j   = J_INF-(J_INF-j)*pa(exptau_j);
  xr1 = Xr1_INF-(Xr1_INF-xr1)*pa(exptau_Xr1);
  xr2 = Xr2_INF-(Xr2_INF-xr2)*pa(exptau_Xr2);
  xs  = Xs_INF-(Xs_INF-xs)*pa(exptau_Xs);
  s   = S_INF-(S_INF-s)*pa(exptau_s);
  r   = R_INF-(R_INF-r)*pa(exptau_r);
  d   = D_INF-(D_INF-d)*pa(exptau_d);
  f   = F_INF-(F_INF-f)*pa(exptau_f);
  f2  = F2_INF-(F2_INF-f2)*pa(exptau_f2);
  fCa = FCa_INF-(FCa_INF-fCa)*exp(-tinc/taufca);

  // lh326 Array for longPrint
  Array[0]  = I_Na;
  Array[1]  = I_CaL;
  Array[2]  = I_bCa;
  Array[3]  = I_pCa;
  Array[4]  = I_to;
  Array[5]  = I_Ks;
  Array[6]  = I_Kr;
  Array[7]  = I_K1;
  Array[8]  = I_pK;
  Array[9]  = I_bNa;
  Array[10] = I_NaK;
  Array[11] = I_NaCa;
  Array[12] = I_rel;
  Array[13] = I_leak;
  Array[14] = I_xfer;
  Array[15] = SERCA;
  Array[16] = I_tot;
#ifdef ACTIVATE_ISAC_CHANNEL
  Array[17] = I_SAC;
#endif //end Activate_ISAC


  return tinc*(-I_tot);
}  // TenTusscherEtAl2::Calc

void TenTusscherEtAl2::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  // Don't forget the blank (' ') at the end!!
  const  ML_CalcType max_sr     = v(VT_max_sr);
  const  ML_CalcType EC         = v(VT_EC);
  const  ML_CalcType kCaSR      = max_sr-((max_sr-v(VT_min_sr))/(1+(EC/CaSR)*(EC/CaSR))); // new n37
  const  ML_CalcType k1         = v(VT_ks1)/kCaSR; // new n35
  const  ML_CalcType CaSSsquare = CaSS*CaSS;       // new
  const  ML_CalcType O          = k1*CaSSsquare*Rq/(v(VT_k3)+k1*CaSSsquare);

  tempstr<<tArg<<' '<<V<<' '
         <<m<<' '<<h<<' '<<j<<' '<<d<<' '
         <<f<<' '<<f2<<' '<<fCa<<' '<<Rq<<' '<<O<<' '<<xr1<<' '<<xr2<<' '
         <<xs<<' '<<r<<' '<<s<<' '<<Ca_i<<' '<<CaSR<<' '<<CaSS<<' '<<Na_i<<' '<<K_i << ' ';

#ifdef MARKOV_I_NA
  tempstr<<MINALC3<<' '<<MINALC2<<' '
         <<MINALC1<<' '<<MINALO<<' '
         <<MINAUC3<<' '<<MINAUC2<<' '
         <<MINAUC1<<' '<<MINAUO<<' '
         <<MINAUIC3<<' '<<MINAUIC2<<' '
         <<MINAUIF<<' '<<MINAUIM1<<' '<<MINAUIM2<<' ';
#endif  // ifdef MARKOV_I_NA
}

void TenTusscherEtAl2::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);

  tempstr
    <<Array[0]<<' '
    <<Array[1]<<' '
    <<Array[2]<<' '
    <<Array[3]<<' '
    <<Array[4]<<' '
    <<Array[5]<<' '
    <<Array[6]<<' '
    <<Array[7]<<' '
    <<Array[8]<<' '
    <<Array[9]<<' '
    <<Array[10]<<' '
    <<Array[11]<<' '
    <<Array[12]<<' '
    <<Array[13]<<' '
    <<Array[14]<<' '
    <<Array[15]<<' '
    <<Array[16]<<' '
#ifdef ACTIVATE_ISAC_CHANNEL
    <<Array[17]<<' '
#endif //end Activate_ISAC
  ;
}  // TenTusscherEtAl2::LongPrint

void TenTusscherEtAl2::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {"m",       "h",                "j",                "d",                 "f",                  "f2",
   "fCa",
   "Rq",
   "O",
   "Xr1",     "Xr2",
   "Xs",      "r",                "s",                "Cai",               "CaSR",               "CaSS",
   "Nai",
   "Ki"
#ifdef MARKOV_I_NA
   ,          "MINALC3",          "MINALC2",          "MINALC1",           "MINALO",             "MINAUC3",
   "MINAUC2",
   "MINAUC1",
   "MINAUO",  "MINAUIC3",         "MINAUIC2",         "MINAUIF",           "MINAUIM1",           "MINAUIM2"
#endif  // ifdef MARKOV_I_NA
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}

void TenTusscherEtAl2::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {
    "I_Na",   "I_CaL",   "I_bCa",   "I_pCa",   "I_to",   "I_Ks", "I_Kr",  "I_K1", "I_pK", "I_bNa",
    "I_NaK",
    "I_NaCa", "I_rel",
    "I_leak", "I_xfer",  "SERCA"
#ifdef ACTIVATE_IKATP_CHANNEL
    ,         "I_Katp"
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
    ,         "I_tot"

#ifdef ACTIVATE_ISAC_CHANNEL
    , "I_SAC"
#endif  // ifdef ACTIVATE_ISAC_CHANNEL
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
