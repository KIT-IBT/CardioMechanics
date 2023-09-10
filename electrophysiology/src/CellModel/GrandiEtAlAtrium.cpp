/* File: GrandiEtAlAtrium.cpp
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#include <GrandiEtAlAtrium.h>

GrandiEtAlAtrium::GrandiEtAlAtrium(GrandiEtAlAtriumParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_GrandiEtAlAtriumParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

GrandiEtAlAtrium::~GrandiEtAlAtrium() {}

#ifdef HETERO

inline bool GrandiEtAlAtrium::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool GrandiEtAlAtrium::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int GrandiEtAlAtrium::GetSize(void) {
  return (&Csqn_b - &Na_j + 1) * sizeof(ML_CalcType);
}

inline unsigned char GrandiEtAlAtrium::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void GrandiEtAlAtrium::Init() {
#if KADEBUG
  cerr << "#initializing Class: GrandiEtAlAtrium ... " << endl;
        #endif  // if KADEBUG
  Na_j     = (v(VT_Na_j_init));
  Na_sl    = (v(VT_Na_sl_init));
  K_i      = (v(VT_K_i_init));
  Ca_j     = (v(VT_Ca_j_init));
  Ca_sl    = (v(VT_Ca_sl_init));
  m        = (v(VT_m_init));
  h        = (v(VT_h_init));
  j        = (v(VT_j_init));
  ml       = (v(VT_ml_init));
  hl       = (v(VT_hl_init));
  x_kr     = (v(VT_x_kr_init));
  x_ks     = (v(VT_x_ks_init));
  Na_i     = (v(VT_Na_i_init));
  x_kur    = (v(VT_x_kur_init));
  y_kur    = (v(VT_y_kur_init));
  x_to_f   = (v(VT_x_to_f_init));
  y_to_f   = (v(VT_y_to_f_init));
  d        = (v(VT_d_init));
  f        = (v(VT_f_init));
  f_Ca_Bj  = (v(VT_f_Ca_Bj_init));
  f_Ca_Bsl = (v(VT_f_Ca_Bsl_init));
  Ry_Rr    = (v(VT_Ry_Rr_init));
  Ry_Ro    = (v(VT_Ry_Ro_init));
  Ry_Ri    = (v(VT_Ry_Ri_init));
  Ca_sr    = (v(VT_Ca_sr_init));
  Ca_i     = (v(VT_Ca_i_init));
  Na_Bj    = (v(VT_Na_Bj_init));
  Na_Bsl   = (v(VT_Na_Bsl_init));
  Tn_CL    = (v(VT_Tn_CL_init));
  Tn_CHc   = (v(VT_Tn_CHc_init));
  Tn_CHm   = (v(VT_Tn_CHm_init));
  CaM      = (v(VT_CaM_init));
  Myo_c    = (v(VT_Myo_c_init));
  Myo_m    = (v(VT_Myo_m_init));
  SRB      = (v(VT_SRB_init));
  SLL_j    = (v(VT_SLL_j_init));
  SLL_sl   = (v(VT_SLL_sl_init));
  SLH_j    = (v(VT_SLH_j_init));
  SLH_sl   = (v(VT_SLH_sl_init));
  Csqn_b   = (v(VT_Csqn_b_init));
}  // GrandiEtAlAtrium::Init

ML_CalcType GrandiEtAlAtrium::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                   int euler) {
  ML_CalcType svolt = V*1000;
  ML_CalcType HT    = tinc*1000;
  const int   Vi    = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);

  const ML_CalcType Cmem = (v(VT_Cmem));

  i_external = i_external/(1e9*Cmem);  // convert from nA to pA/pF

  // calculating algebraic part

  const ML_CalcType FoRT         = (v(VT_FoRT));
  const ML_CalcType Fjunc        = (v(VT_Fjunc));
  const ML_CalcType K_o          = (v(VT_K_o));
  const ML_CalcType Na_o         = (v(VT_Na_o));
  const ML_CalcType Ca_o         = (v(VT_Ca_o));
  const ML_CalcType Frdy         = (v(VT_Frdy));
  const ML_CalcType Fsl          = (v(VT_Fsl));
  const ML_CalcType koff_na      = (v(VT_koff_na));
  const ML_CalcType kon_na       = (v(VT_kon_na));
  const ML_CalcType kon_sll      = (v(VT_kon_sll));
  const ML_CalcType koff_sll     = (v(VT_koff_sll));
  const ML_CalcType kon_slh      = (v(VT_kon_slh));
  const ML_CalcType koff_slh     = (v(VT_koff_slh));
  const ML_CalcType Bmax_myosin  = (v(VT_Bmax_myosin));
  const ML_CalcType Bmax_TnChigh = (v(VT_Bmax_TnChigh));
  const ML_CalcType Mgi          = (v(VT_Mgi));
  const ML_CalcType Vmyo         = (v(VT_Vmyo));
  const ML_CalcType Vsr          = (v(VT_Vsr));
  const ML_CalcType Vjunc        = (v(VT_Vjunc));
  const ML_CalcType Vsl          = (v(VT_Vsl));
  const ML_CalcType kom          = (v(VT_kom));
  const ML_CalcType kim          = (v(VT_kim));
  const ML_CalcType expVFRT      = exp(svolt*FoRT);
  const ML_CalcType exp2VFRT     = expVFRT*expVFRT;
  const ML_CalcType expnuVFRT    = exp(( (v(VT_nu))*svolt*FoRT));
  const ML_CalcType expnum1VFRT  = expnuVFRT/expVFRT;

  const ML_CalcType kCaSR = (v(VT_MaxSR)) - ((v(VT_MaxSR)) - (v(VT_MinSR)))/
    (1.00000+(pow(((v(VT_ec50SR))/Ca_sr), 2.50000)));
  const ML_CalcType koSRCa    = (v(VT_koCa))/kCaSR;
  const ML_CalcType kiSRCa    =  (v(VT_kiCa))*kCaSR;
  const ML_CalcType RI        = ((1.00000 - Ry_Rr) - Ry_Ro) - Ry_Ri;
  const ML_CalcType J_SRCarel =  (( (v(VT_ks))*Ry_Ro)/1.00000)*(Ca_sr - Ca_j);
  const ML_CalcType powCai    = pow((Ca_i/(v(VT_Kmf))), (v(VT_hillSRCaP)));
  const ML_CalcType powCasr   = pow((Ca_sr/(v(VT_Kmr))), (v(VT_hillSRCaP)));
  const ML_CalcType J_serca   = ( ((v(VT_powsrcap)))*(v(VT_Vmax_SRCaP))*((powCai) - (powCasr)))/
    (1.00000+(powCai)+(powCasr));
  const ML_CalcType J_SRleak =  (1.0+0.25*(v(VT_AF)))*5.34800e-06*(Ca_sr - Ca_j);

  const ML_CalcType dTn_CL  = (v(VT_kon_tncl))*Ca_i*((v(VT_Bmax_TnClow)) - Tn_CL) -  (v(VT_koff_tncl))*Tn_CL;
  const ML_CalcType dTn_CHc = (v(VT_kon_tnchca))*Ca_i*((Bmax_TnChigh - Tn_CHc) - Tn_CHm) -  (v(VT_koff_tnchca))*
    Tn_CHc;
  const ML_CalcType dTn_CHm = (v(VT_kon_tnchmg))*Mgi*((Bmax_TnChigh - Tn_CHc) - Tn_CHm) -  (v(VT_koff_tnchmg))*
    Tn_CHm;
  const ML_CalcType dCaM          = (v(VT_kon_cam))*Ca_i*((v(VT_Bmax_CaM)) - CaM) -  (v(VT_koff_cam))*CaM;
  const ML_CalcType dMyo_c        = (v(VT_kon_myoca))*Ca_i*((Bmax_myosin - Myo_c) - Myo_m) -  (v(VT_koff_myoca))*Myo_c;
  const ML_CalcType dMyo_m        = (v(VT_kon_myomg))*Mgi*((Bmax_myosin - Myo_c) - Myo_m) -  (v(VT_koff_myomg))*Myo_m;
  const ML_CalcType dSRB          = (v(VT_kon_sr))*Ca_i*((v(VT_Bmax_SR)) - SRB) -  (v(VT_koff_sr))*SRB;
  const ML_CalcType J_CaB_cytosol = dTn_CL+dTn_CHc+dTn_CHm+dCaM+dMyo_c+dMyo_m+dSRB;
  const ML_CalcType dNa_Bj_dt     =  kon_na*Na_j*((v(VT_Bmax_Naj)) - Na_Bj) -  koff_na*Na_Bj;
  const ML_CalcType dNa_Bsl_dt    =  kon_na*Na_sl*((v(VT_Bmax_Nasl)) - Na_Bsl) -  koff_na*Na_Bsl;

  const ML_CalcType ena_junc    =  (1.00000/FoRT)*(log((Na_o/Na_j)));
  const ML_CalcType I_naconst   = (v(VT_GNa))*m*m*m*h*j;
  const ML_CalcType I_Na_junc   =  Fjunc*I_naconst*(svolt - ena_junc);
  const ML_CalcType I_naLconst  = (v(VT_GNaL))*ml*ml*ml*hl;
  const ML_CalcType I_NaL_junc  = Fjunc*I_naLconst*(svolt - ena_junc);
  const ML_CalcType I_nabk_junc =  Fjunc*(v(VT_GNaB))*(svolt - ena_junc);
  ML_CalcType pownaj            = (v(VT_KmNaip))/Na_j;
  pownaj *= pownaj;
  pownaj *= pownaj;
  const ML_CalcType I_nakconst  = (v(VT_IbarNaK))*ptTeaP->fnak[Vi]*K_o/(K_o+(v(VT_KmKo)));
  const ML_CalcType I_nak_junc  = ((Fjunc*I_nakconst)/(1.00000+(pownaj)));
  const ML_CalcType dfconst     = svolt*Frdy*FoRT*d*f*((v(VT_powcal)))*0.450000;
  const ML_CalcType ibarnaconst = (v(VT_pNa))*dfconst*0.75/((expVFRT) - 1.00000);
  const ML_CalcType ibarna_j    = ibarnaconst*(Na_j*(expVFRT) - Na_o);
  const ML_CalcType I_CaNa_junc =  (v(VT_Fjunc_CaL))*ibarna_j*((1.00000 - f_Ca_Bj)+(v(VT_fcaCaj)));
  const ML_CalcType Ka_junc     = 1.00000/(1.00000+(pow(((v(VT_Kdact))/Ca_j), 2.00000)));
  const ML_CalcType powNaj      = Na_j*Na_j*Na_j;
  const ML_CalcType powNao      = Na_o*Na_o*Na_o;
  const ML_CalcType s1_junc     =  (expnuVFRT)*(powNaj)*Ca_o;
  const ML_CalcType s2_junc     =  (expnum1VFRT)*(powNao)*Ca_j;
  const ML_CalcType s3_junc     =  (v(VT_KmCai))*(powNao)*
    (1.00000+
     (pow((Na_j/(v(VT_KmNai))),
          3.00000)))+
    (pow((v(VT_KmNao)),
         3.00000))*Ca_j*(1.00000+Ca_j/(v(VT_KmCai)))+ (v(VT_KmCao))*(powNaj)+ (powNaj)*Ca_o+ (powNao)*Ca_j;
  const ML_CalcType incxconst     = (v(VT_IbarNCX))*(v(VT_powncx))/(1.00000+ (v(VT_ksat))*(expnum1VFRT));
  const ML_CalcType I_ncx_junc    = ((Fjunc*incxconst*Ka_junc*(s1_junc - s2_junc))/s3_junc);
  const ML_CalcType I_Na_tot_junc = I_Na_junc+I_nabk_junc+ 3.00000*I_ncx_junc+ 3.00000*I_nak_junc+I_CaNa_junc+
    I_NaL_junc;

  const ML_CalcType ena_sl    =  (1.00000/FoRT)*(log((Na_o/Na_sl)));
  const ML_CalcType I_Na_sl   =  Fsl*I_naconst*(svolt - ena_sl);
  const ML_CalcType I_NaL_sl  = Fsl*I_naLconst*(svolt - ena_sl);
  const ML_CalcType I_nabk_sl =  Fsl*(v(VT_GNaB))*(svolt - ena_sl);
  const ML_CalcType I_nak_sl  = ((Fsl*I_nakconst)/(1.00000+(pow(((v(VT_KmNaip))/Na_sl), 4.00000))));
  const ML_CalcType ibarna_sl = ibarnaconst*(Na_sl*(expVFRT) - Na_o);
  const ML_CalcType I_CaNa_sl =  (v(VT_Fsl_CaL))*ibarna_sl*((1.00000 - f_Ca_Bsl)+(v(VT_fcaCaMSL)));
  const ML_CalcType Ka_sl     = 1.00000/(1.00000+(pow(((v(VT_Kdact))/Ca_sl), 2.00000)));
  const ML_CalcType powNasl   = Na_sl*Na_sl*Na_sl;
  const ML_CalcType s1_sl     =  (expnuVFRT)*(powNasl)*Ca_o;
  const ML_CalcType s2_sl     =  (expnum1VFRT)*(powNao)*Ca_sl;
  const ML_CalcType s3_sl     =  (v(VT_KmCai))*(powNao)*
    (1.00000+
     (pow((Na_sl/(v(VT_KmNai))),
          3.00000)))+
    (pow((v(VT_KmNao)),
         3.00000))*Ca_sl*(1.00000+Ca_sl/(v(VT_KmCai)))+ (v(VT_KmCao))*(powNasl)+ (powNasl)*Ca_o+ (powNao)*Ca_sl;
  const ML_CalcType I_ncx_sl    = ((Fsl*incxconst*Ka_sl*(s1_sl - s2_sl))/s3_sl);
  const ML_CalcType I_Na_tot_sl = I_Na_sl+I_nabk_sl+ 3.00000*I_ncx_sl+ 3.00000*I_nak_sl+I_CaNa_sl+I_NaL_sl;

  const ML_CalcType dSLL_j  = kon_sll*Ca_j*((v(VT_Bmax_SLlowj)) - SLL_j) -  koff_sll*SLL_j;
  const ML_CalcType dSLL_sl = kon_sll*Ca_sl*((v(VT_Bmax_SLlowsl)) - SLL_sl) -  koff_sll*SLL_sl;
  const ML_CalcType dSLH_j  = kon_slh*Ca_j*((v(VT_Bmax_SLhighj)) - SLH_j) -  koff_slh*SLH_j;
  const ML_CalcType dSLH_sl = kon_slh*Ca_sl*((v(VT_Bmax_SLhighsl)) - SLH_sl) -  koff_slh*SLH_sl;

  const ML_CalcType J_CaB_junction = dSLL_j+dSLH_j;
  const ML_CalcType ibarcaconst    = (v(VT_pCa))*0.341000*4.00000/((exp2VFRT) - 1.00000)*dfconst;
  const ML_CalcType ibarca_j       = ibarcaconst*(Ca_j*(exp2VFRT) - Ca_o);
  const ML_CalcType I_Ca_junc      =  (v(VT_Fjunc_CaL))*ibarca_j*((1.00000 - f_Ca_Bj)+(v(VT_fcaCaj)));
  const ML_CalcType powCaj         = (pow(Ca_j, 1.60000));
  const ML_CalcType I_pca_junc     = (Fjunc*((v(VT_powslcap)))*(v(VT_IbarSLCaP))*powCaj)/((v(VT_powKmPCa))+powCaj);
  const ML_CalcType eca_junc       =  0.5/FoRT*log(Ca_o/Ca_j);
  const ML_CalcType I_cabk_junc    =  Fjunc*(v(VT_GCaB))*(svolt - eca_junc);
  const ML_CalcType I_Ca_tot_junc  = (I_Ca_junc+I_cabk_junc+I_pca_junc) -  2.00000*I_ncx_junc;

  const ML_CalcType J_CaB_sl    = dSLL_sl+dSLH_sl;
  const ML_CalcType ibarca_sl   = ibarcaconst*(Ca_sl*(exp2VFRT) - Ca_o);
  const ML_CalcType I_Ca_sl     =  (v(VT_Fsl_CaL))*ibarca_sl*((1.00000 - f_Ca_Bsl)+(v(VT_fcaCaMSL)));
  const ML_CalcType powCasl     = (pow(Ca_sl, 1.60000));
  const ML_CalcType I_pca_sl    = (Fsl*((v(VT_powslcap)))*(v(VT_IbarSLCaP))*powCasl)/((v(VT_powKmPCa))+powCasl);
  const ML_CalcType eca_sl      =  0.5/FoRT*log(Ca_o/Ca_sl);
  const ML_CalcType I_cabk_sl   =  Fsl*(v(VT_GCaB))*(svolt - eca_sl);
  const ML_CalcType I_Ca_tot_sl = (I_Ca_sl+I_cabk_sl+I_pca_sl) -  2.00000*I_ncx_sl;

  const ML_CalcType I_nak     = I_nak_junc+I_nak_sl;
  const ML_CalcType ek        =  (1.00000/FoRT)*(log((K_o/K_i)));
  const ML_CalcType I_kr      =  (v(VT_gkr))*x_kr*ptTeaP->rkr[Vi]*(svolt - ek);
  const ML_CalcType ikpconst  = (v(VT_gkp))*ptTeaP->kp_kp[Vi]*(svolt - ek);
  const ML_CalcType I_kp_junc =  Fjunc*ikpconst;
  const ML_CalcType I_kp_sl   =  Fsl*ikpconst;
  const ML_CalcType I_kp      = I_kp_junc+I_kp_sl;
  const ML_CalcType eks       =  (1.00000/FoRT)*(log(((K_o+ (v(VT_pNaK))*Na_o)/(K_i+ (v(VT_pNaK))*Na_i))));
  const ML_CalcType iksconst  = x_ks*x_ks*(svolt - eks);
  const ML_CalcType I_ks_junc =  Fjunc*(v(VT_gks_junc))*iksconst;
  const ML_CalcType I_ks_sl   =  Fsl*(v(VT_gks_sl))*iksconst;
  const ML_CalcType I_ks      = I_ks_junc+I_ks_sl;
  const ML_CalcType I_tof     =  (v(VT_GtoFast))*x_to_f*y_to_f*(svolt - ek);
  const ML_CalcType I_to      = I_tof;
  const ML_CalcType I_kur     = (v(VT_Gkur))*x_kur*y_kur*(svolt-ek);
  const ML_CalcType aki       = 1.02000/(1.00000+(exp((0.238500*((svolt - ek) - 59.2150)))));
  const ML_CalcType bki       =
    (0.491240*(exp((0.0803200*((svolt+5.47600) - ek))))+(exp((0.0617500*((svolt - ek) - 594.310)))))/
    (1.00000+(exp((-0.514300*((svolt - ek)+4.75300)))));
  const ML_CalcType kiss  = aki/(aki+bki);
  const ML_CalcType I_ki  =  (v(VT_Ikiconst))*kiss*(svolt - ek);
  const ML_CalcType ibark = ( (v(VT_pK))*(0.750000*K_i*(expVFRT) -  0.750000*K_o))/((expVFRT) - 1.00000);
  const ML_CalcType I_CaK =  ibark*dfconst*
    ( (v(VT_Fjunc_CaL))*((v(VT_fcaCaj))+(1.00000 - f_Ca_Bj))+ (v(VT_Fsl_CaL))*((v(VT_fcaCaMSL))+(1.00000 - f_Ca_Bsl)));
  const ML_CalcType I_K_tot        = ((I_to+I_kr+I_ks+I_ki) -  2.00000*I_nak)+I_CaK+I_kp+I_kur;
  const ML_CalcType I_Na_tot       = I_Na_tot_junc+I_Na_tot_sl;
  const ML_CalcType iclcaconst     = (v(VT_GClCa))*(svolt - (v(VT_ecl)));
  const ML_CalcType I_ClCa_junc    =  ((Fjunc*iclcaconst)/(1.00000+(v(VT_KdClCa))/Ca_j));
  const ML_CalcType I_ClCa_sl      =  ((Fsl*iclcaconst)/(1.00000+(v(VT_KdClCa))/Ca_sl));
  const ML_CalcType I_ClCa         = I_ClCa_junc+I_ClCa_sl;
  const ML_CalcType I_Cl_tot       = I_ClCa+ptTeaP->I_Clbk[Vi];
  const ML_CalcType I_Ca_tot       = I_Ca_tot_junc+I_Ca_tot_sl;
  const ML_CalcType I_tot          = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot-i_external;
  const ML_CalcType I_Na           = I_Na_junc+I_Na_sl;
  const ML_CalcType I_nabk         = I_nabk_junc+I_nabk_sl;
  const ML_CalcType I_Ca           = I_Ca_junc+I_Ca_sl;
  const ML_CalcType I_CaNa         = I_CaNa_junc+I_CaNa_sl;
  const ML_CalcType I_Catot        = I_Ca+I_CaK+I_CaNa;
  const ML_CalcType I_Na_tot_junc2 =  3.00000*I_ncx_junc+ 3.00000*I_nak_junc+I_CaNa_junc;
  const ML_CalcType I_ncx          = I_ncx_junc+I_ncx_sl;
  const ML_CalcType I_Na_tot_sl2   =  3.00000*I_ncx_sl+ 3.00000*I_nak_sl+I_CaNa_sl;
  const ML_CalcType I_pca          = I_pca_junc+I_pca_sl;
  const ML_CalcType I_cabk         = I_cabk_junc+I_cabk_sl;

  // calculating rates part

  int NDERIVS          = 4; // number of derivatives/accuracy
  const ML_CalcType h1 = HT/(NDERIVS)*0.5;
  const ML_CalcType h2 = h1*2.0;

  const ML_CalcType srconst1 = (J_serca - ((J_SRleak*Vmyo)/Vsr+J_SRCarel));
  const ML_CalcType srconst2 = (v(VT_kon_csqn));
  const ML_CalcType srconst3 = (v(VT_Bmax_Csqn));
  const ML_CalcType srconst4 = (v(VT_koff_csqn));
  const ML_CalcType iconst1  = ((-J_serca*Vsr)/Vmyo - J_CaB_cytosol);
  const ML_CalcType iconst2  = ((v(VT_J_ca_slmyo))/Vmyo);

  volatile ML_CalcType Ca_sr1, Ca_sr2, Ca_i1, Ca_i2, Csqn_b1, Csqn_b2;  // derivatives for each variable at beginning (1) and end
  // (2) of time step
  for (int i = 0; i < NDERIVS; i++) {
    Ca_sr1 = (srconst1 - (srconst2*Ca_sr*(srconst3 - Csqn_b) -  srconst4*Csqn_b));
    Ca_sr2 = (srconst1 - (srconst2*(Ca_sr+h1*Ca_sr1)*(srconst3 - (Csqn_b+h1*Csqn_b2)) -  srconst4*(Csqn_b+h1*Csqn_b2)));

    Ca_i1 = (iconst1+ iconst2*(Ca_sl - Ca_i));
    Ca_i2 = (iconst1+ iconst2*(Ca_sl - (Ca_i+h1*Ca_i1)));

    Csqn_b1 = srconst2*Ca_sr*(srconst3 - Csqn_b) -  srconst4*Csqn_b;
    Csqn_b2 = srconst2*(Ca_sr+h1*Ca_sr2)*(srconst3 - (Csqn_b+h1*Csqn_b1)) -  srconst4*(Csqn_b+h1*Csqn_b1);

    Ca_sr  += h2*Ca_sr2;
    Ca_i   += h2*Ca_i2;
    Csqn_b += h2*Csqn_b2;
  }

  Ry_Ro    += HT*( (koSRCa*Ca_j*Ca_j*Ry_Rr -  kom*Ry_Ro) - (kiSRCa*Ca_j*Ry_Ro -  kim*Ry_Ri));
  Ry_Rr    += HT*( (kim*RI -  kiSRCa*Ca_j*Ry_Rr) - (koSRCa*Ca_j*Ca_j*Ry_Rr -  kom*Ry_Ro));
  Ry_Ri    += HT*( (kiSRCa*Ca_j*Ry_Ro -  kim*Ry_Ri) - (kom*Ry_Ri -  koSRCa*Ca_j*Ca_j*RI));
  Tn_CL    += HT*dTn_CL;
  Tn_CHc   += HT*dTn_CHc;
  Tn_CHm   += HT*dTn_CHm;
  CaM      += HT*dCaM;
  Myo_c    += HT*dMyo_c;
  Myo_m    += HT*dMyo_m;
  SRB      += HT*dSRB;
  f_Ca_Bj  += HT*(  ((1.70000*Ca_j))*(1.00000 - f_Ca_Bj) -  0.0119000*f_Ca_Bj);
  f_Ca_Bsl += HT*(  ((1.70000*Ca_sl))*(1.00000 - f_Ca_Bsl) -  0.0119000*f_Ca_Bsl);
  SLL_j    += HT*dSLL_j;
  SLL_sl   += HT*dSLL_sl;
  SLH_j    += HT*dSLH_j;
  SLH_sl   += HT*dSLH_sl;

  // m += HT*( (ptTeaP->mss[Vi] - m)/ptTeaP->taum[Vi]);
  m = ptTeaP->mss[Vi]+(m-ptTeaP->mss[Vi])*ptTeaP->exptaum[Vi];

  // x_kr += HT*( (ptTeaP->xrss[Vi] - x_kr)/ptTeaP->tauxr[Vi]);
  x_kr = ptTeaP->xrss[Vi]+(x_kr-ptTeaP->xrss[Vi])*ptTeaP->exptauxr[Vi];

  // x_ks += HT*( (ptTeaP->xsss[Vi] - x_ks)/ptTeaP->tauxs[Vi]);
  x_ks  = ptTeaP->xsss[Vi]+(x_ks-ptTeaP->xsss[Vi])*ptTeaP->exptauxs[Vi];
  x_kur = ptTeaP->xkurss[Vi]+(x_kur-ptTeaP->xkurss[Vi])*ptTeaP->exptauxkur[Vi];
  y_kur = ptTeaP->ykurss[Vi]+(y_kur-ptTeaP->ykurss[Vi])*ptTeaP->exptauykur[Vi];

  // x_to_f += HT*( (ptTeaP->xtoss[Vi] - x_to_f)/ptTeaP->tauxtof[Vi]);
  x_to_f = ptTeaP->xtoss[Vi]+(x_to_f-ptTeaP->xtoss[Vi])*ptTeaP->exptauxtof[Vi];

  // y_to_f += HT*( (ptTeaP->ytoss[Vi] - y_to_f)/ptTeaP->tauytof[Vi]);
  y_to_f = ptTeaP->ytoss[Vi]+(y_to_f-ptTeaP->ytoss[Vi])*ptTeaP->exptauytof[Vi];

  // d += HT*( (ptTeaP->dss[Vi] - d)/taud);
  d = ptTeaP->dss[Vi]+(d-ptTeaP->dss[Vi])*ptTeaP->exptaud[Vi];

  // f += HT*( (ptTeaP->fss[Vi] - f)/ptTeaP->tauf[Vi]);
  f = ptTeaP->fss[Vi]+(f-ptTeaP->fss[Vi])*ptTeaP->exptauf[Vi];

  // h += HT*( (ptTeaP->hss[Vi] - h)/tauh);
  h = ptTeaP->hss[Vi]+(h-ptTeaP->hss[Vi])*ptTeaP->exptauh[Vi];

  // j += HT*( (ptTeaP->jss[Vi] - j)/tauj);
  j  = ptTeaP->jss[Vi]+(j-ptTeaP->jss[Vi])*ptTeaP->exptauj[Vi];
  ml = ptTeaP->mlss[Vi]+(ml-ptTeaP->mlss[Vi])*ptTeaP->exptauml[Vi];
  hl = ptTeaP->hlss[Vi]+(hl-ptTeaP->hlss[Vi])*exp(-HT/600);

  Na_Bj  += HT*(dNa_Bj_dt);
  Na_Bsl += HT*(dNa_Bsl_dt);
  Na_j   += HT*( ((-I_Na_tot_junc*Cmem)/(Vjunc*Frdy)+ ((v(VT_J_na_juncsl))/Vjunc)*(Na_sl - Na_j)) - dNa_Bj_dt);
  Na_sl  += HT*
    ( ((-I_Na_tot_sl*Cmem)/(Vsl*Frdy)+ ((v(VT_J_na_juncsl))/Vsl)*(Na_j - Na_sl)+ ((v(VT_J_na_slmyo))/Vsl)*
       (Na_i - Na_sl)) - dNa_Bsl_dt);
  Na_i  += HT*(  ((v(VT_J_na_slmyo))/Vmyo)*(Na_sl - Na_i));
  Ca_sl += HT*
    ( ((-I_Ca_tot_sl*Cmem)/(Vsl*2.00000*Frdy)+ ((v(VT_J_ca_juncsl))/Vsl)*(Ca_j - Ca_sl)+ ((v(VT_J_ca_slmyo))/Vsl)*
       (Ca_i - Ca_sl)) - J_CaB_sl);
  Ca_j += HT*
    ( (((-I_Ca_tot_junc*Cmem)/(Vjunc*2.00000*Frdy)+ ((v(VT_J_ca_juncsl))/Vjunc)*(Ca_sl - Ca_j)) - J_CaB_junction)+
      (J_SRCarel*Vsr)/Vjunc+(J_SRleak*Vmyo)/Vjunc);

  return tinc*(-I_tot);
}  // GrandiEtAlAtrium::Calc

void GrandiEtAlAtrium::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<<tArg<< ' ' << V<<' ' <<Na_j<<' ' <<Na_sl<<' ' <<K_i<<' ' <<Ca_j<<' ' <<Ca_sl<<' ' <<m<<' ' <<h<<' ' <<j<<
    ' ' <<x_kr
         <<' ' <<x_ks<<' ' <<Na_i<<' ' <<ml<<' ' <<hl<<' ' <<x_kur<<' ' <<y_kur<<' ' <<x_to_f<<' ' <<y_to_f<<' ' <<d<<
    ' ' <<f<<' ' <<f_Ca_Bj
         <<' ' <<f_Ca_Bsl<<' ' <<Ry_Rr<<' ' <<Ry_Ro<<' ' <<Ry_Ri<<' ' <<Ca_sr<<' ' <<Ca_i<<' ' <<Na_Bj<<' ' <<Na_Bsl<<
    ' ' <<Tn_CL
         <<' ' <<Tn_CHc<<' ' <<Tn_CHm<<' ' <<CaM<<' ' <<Myo_c<<' ' <<Myo_m<<' ' <<SRB<<' ' <<SLL_j<<' ' <<SLL_sl<<' ' <<
    SLH_j
         <<' ' <<SLH_sl<<' ' <<Csqn_b;
}

void GrandiEtAlAtrium::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
  const  ML_CalcType svolt       = V*1000.0;
  const int Vi                   = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);
  const ML_CalcType FoRT         = (v(VT_FoRT));
  const ML_CalcType Fjunc        = (v(VT_Fjunc));
  const ML_CalcType K_o          = (v(VT_K_o));
  const ML_CalcType Na_o         = (v(VT_Na_o));
  const ML_CalcType Ca_o         = (v(VT_Ca_o));
  const ML_CalcType Frdy         = (v(VT_Frdy));
  const ML_CalcType Fsl          = (v(VT_Fsl));
  const ML_CalcType Cmem         = (v(VT_Cmem));
  const ML_CalcType koff_na      = (v(VT_koff_na));
  const ML_CalcType kon_na       = (v(VT_kon_na));
  const ML_CalcType kon_sll      = (v(VT_kon_sll));
  const ML_CalcType koff_sll     = (v(VT_koff_sll));
  const ML_CalcType kon_slh      = (v(VT_kon_slh));
  const ML_CalcType koff_slh     = (v(VT_koff_slh));
  const ML_CalcType Bmax_myosin  = (v(VT_Bmax_myosin));
  const ML_CalcType Bmax_TnChigh = (v(VT_Bmax_TnChigh));
  const ML_CalcType Mgi          = (v(VT_Mgi));
  const ML_CalcType Vmyo         = (v(VT_Vmyo));
  const ML_CalcType Vsr          = (v(VT_Vsr));
  const ML_CalcType Vjunc        = (v(VT_Vjunc));
  const ML_CalcType Vsl          = (v(VT_Vsl));
  const ML_CalcType kom          = (v(VT_kom));
  const ML_CalcType kim          = (v(VT_kim));
  const ML_CalcType expVFRT      = exp(svolt*FoRT);
  const ML_CalcType exp2VFRT     = expVFRT*expVFRT;
  const ML_CalcType expnuVFRT    = exp(( (v(VT_nu))*svolt*FoRT));
  const ML_CalcType expnum1VFRT  = expnuVFRT/expVFRT;

  const ML_CalcType kCaSR = (v(VT_MaxSR)) - ((v(VT_MaxSR)) - (v(VT_MinSR)))/
    (1.00000+(pow(((v(VT_ec50SR))/Ca_sr), 2.50000)));
  const ML_CalcType koSRCa    = (v(VT_koCa))/kCaSR;
  const ML_CalcType kiSRCa    =  (v(VT_kiCa))*kCaSR;
  const ML_CalcType RI        = ((1.00000 - Ry_Rr) - Ry_Ro) - Ry_Ri;
  const ML_CalcType J_SRCarel =  (( (v(VT_ks))*Ry_Ro)/1.00000)*(Ca_sr - Ca_j);
  const ML_CalcType powCai    = pow((Ca_i/(v(VT_Kmf))), (v(VT_hillSRCaP)));
  const ML_CalcType powCasr   = pow((Ca_sr/(v(VT_Kmr))), (v(VT_hillSRCaP)));
  const ML_CalcType J_serca   = ( ((v(VT_powsrcap)))*(v(VT_Vmax_SRCaP))*((powCai) - (powCasr)))/
    (1.00000+(powCai)+(powCasr));
  const ML_CalcType J_SRleak =  (1.0+0.25*(v(VT_AF)))*5.34800e-06*(Ca_sr - Ca_j);

  const ML_CalcType dTn_CL  = (v(VT_kon_tncl))*Ca_i*((v(VT_Bmax_TnClow)) - Tn_CL) -  (v(VT_koff_tncl))*Tn_CL;
  const ML_CalcType dTn_CHc = (v(VT_kon_tnchca))*Ca_i*((Bmax_TnChigh - Tn_CHc) - Tn_CHm) -  (v(VT_koff_tnchca))*
    Tn_CHc;
  const ML_CalcType dTn_CHm = (v(VT_kon_tnchmg))*Mgi*((Bmax_TnChigh - Tn_CHc) - Tn_CHm) -  (v(VT_koff_tnchmg))*
    Tn_CHm;
  const ML_CalcType dCaM          = (v(VT_kon_cam))*Ca_i*((v(VT_Bmax_CaM)) - CaM) -  (v(VT_koff_cam))*CaM;
  const ML_CalcType dMyo_c        = (v(VT_kon_myoca))*Ca_i*((Bmax_myosin - Myo_c) - Myo_m) -  (v(VT_koff_myoca))*Myo_c;
  const ML_CalcType dMyo_m        = (v(VT_kon_myomg))*Mgi*((Bmax_myosin - Myo_c) - Myo_m) -  (v(VT_koff_myomg))*Myo_m;
  const ML_CalcType dSRB          = (v(VT_kon_sr))*Ca_i*((v(VT_Bmax_SR)) - SRB) -  (v(VT_koff_sr))*SRB;
  const ML_CalcType J_CaB_cytosol = dTn_CL+dTn_CHc+dTn_CHm+dCaM+dMyo_c+dMyo_m+dSRB;
  const ML_CalcType dNa_Bj_dt     =  kon_na*Na_j*((v(VT_Bmax_Naj)) - Na_Bj) -  koff_na*Na_Bj;
  const ML_CalcType dNa_Bsl_dt    =  kon_na*Na_sl*((v(VT_Bmax_Nasl)) - Na_Bsl) -  koff_na*Na_Bsl;

  const ML_CalcType ena_junc    =  (1.00000/FoRT)*(log((Na_o/Na_j)));
  const ML_CalcType I_naconst   = (v(VT_GNa))*m*m*m*h*j;
  const ML_CalcType I_Na_junc   =  Fjunc*I_naconst*(svolt - ena_junc);
  const ML_CalcType I_naLconst  = (v(VT_GNaL))*ml*ml*ml*hl;
  const ML_CalcType I_NaL_junc  = Fjunc*I_naLconst*(svolt - ena_junc);
  const ML_CalcType I_nabk_junc =  Fjunc*(v(VT_GNaB))*(svolt - ena_junc);
  ML_CalcType pownaj            = (v(VT_KmNaip))/Na_j;
  pownaj *= pownaj;
  pownaj *= pownaj;
  const ML_CalcType I_nakconst  = (v(VT_IbarNaK))*ptTeaP->fnak[Vi]*K_o/(K_o+(v(VT_KmKo)));
  const ML_CalcType I_nak_junc  = ((Fjunc*I_nakconst)/(1.00000+(pownaj)));
  const ML_CalcType dfconst     = svolt*Frdy*FoRT*d*f*((v(VT_powcal)))*0.450000;
  const ML_CalcType ibarnaconst = (v(VT_pNa))*dfconst*0.75/((expVFRT) - 1.00000);
  const ML_CalcType ibarna_j    = ibarnaconst*(Na_j*(expVFRT) - Na_o);
  const ML_CalcType I_CaNa_junc =  (v(VT_Fjunc_CaL))*ibarna_j*((1.00000 - f_Ca_Bj)+(v(VT_fcaCaj)));
  const ML_CalcType Ka_junc     = 1.00000/(1.00000+(pow(((v(VT_Kdact))/Ca_j), 2.00000)));
  const ML_CalcType powNaj      = Na_j*Na_j*Na_j;
  const ML_CalcType powNao      = Na_o*Na_o*Na_o;
  const ML_CalcType s1_junc     =  (expnuVFRT)*(powNaj)*Ca_o;
  const ML_CalcType s2_junc     =  (expnum1VFRT)*(powNao)*Ca_j;
  const ML_CalcType s3_junc     =  (v(VT_KmCai))*(powNao)*
    (1.00000+
     (pow((Na_j/(v(VT_KmNai))),
          3.00000)))+
    (pow((v(VT_KmNao)),
         3.00000))*Ca_j*(1.00000+Ca_j/(v(VT_KmCai)))+ (v(VT_KmCao))*(powNaj)+ (powNaj)*Ca_o+ (powNao)*Ca_j;
  const ML_CalcType incxconst     = (v(VT_IbarNCX))*(v(VT_powncx))/(1.00000+ (v(VT_ksat))*(expnum1VFRT));
  const ML_CalcType I_ncx_junc    = ((Fjunc*incxconst*Ka_junc*(s1_junc - s2_junc))/s3_junc);
  const ML_CalcType I_Na_tot_junc = I_Na_junc+I_nabk_junc+ 3.00000*I_ncx_junc+ 3.00000*I_nak_junc+I_CaNa_junc+
    I_NaL_junc;

  const ML_CalcType ena_sl    =  (1.00000/FoRT)*(log((Na_o/Na_sl)));
  const ML_CalcType I_Na_sl   =  Fsl*I_naconst*(svolt - ena_sl);
  const ML_CalcType I_NaL_sl  = Fsl*I_naLconst*(svolt - ena_sl);
  const ML_CalcType I_nabk_sl =  Fsl*(v(VT_GNaB))*(svolt - ena_sl);
  const ML_CalcType I_nak_sl  = ((Fsl*I_nakconst)/(1.00000+(pow(((v(VT_KmNaip))/Na_sl), 4.00000))));
  const ML_CalcType ibarna_sl = ibarnaconst*(Na_sl*(expVFRT) - Na_o);
  const ML_CalcType I_CaNa_sl =  (v(VT_Fsl_CaL))*ibarna_sl*((1.00000 - f_Ca_Bsl)+(v(VT_fcaCaMSL)));
  const ML_CalcType Ka_sl     = 1.00000/(1.00000+(pow(((v(VT_Kdact))/Ca_sl), 2.00000)));
  const ML_CalcType powNasl   = Na_sl*Na_sl*Na_sl;
  const ML_CalcType s1_sl     =  (expnuVFRT)*(powNasl)*Ca_o;
  const ML_CalcType s2_sl     =  (expnum1VFRT)*(powNao)*Ca_sl;
  const ML_CalcType s3_sl     =  (v(VT_KmCai))*(powNao)*
    (1.00000+
     (pow((Na_sl/(v(VT_KmNai))),
          3.00000)))+
    (pow((v(VT_KmNao)),
         3.00000))*Ca_sl*(1.00000+Ca_sl/(v(VT_KmCai)))+ (v(VT_KmCao))*(powNasl)+ (powNasl)*Ca_o+ (powNao)*Ca_sl;
  const ML_CalcType I_ncx_sl    = ((Fsl*incxconst*Ka_sl*(s1_sl - s2_sl))/s3_sl);
  const ML_CalcType I_Na_tot_sl = I_Na_sl+I_nabk_sl+ 3.00000*I_ncx_sl+ 3.00000*I_nak_sl+I_CaNa_sl+I_NaL_sl;

  const ML_CalcType dSLL_j  = kon_sll*Ca_j*((v(VT_Bmax_SLlowj)) - SLL_j) -  koff_sll*SLL_j;
  const ML_CalcType dSLL_sl = kon_sll*Ca_sl*((v(VT_Bmax_SLlowsl)) - SLL_sl) -  koff_sll*SLL_sl;
  const ML_CalcType dSLH_j  = kon_slh*Ca_j*((v(VT_Bmax_SLhighj)) - SLH_j) -  koff_slh*SLH_j;
  const ML_CalcType dSLH_sl = kon_slh*Ca_sl*((v(VT_Bmax_SLhighsl)) - SLH_sl) -  koff_slh*SLH_sl;

  const ML_CalcType J_CaB_junction = dSLL_j+dSLH_j;
  const ML_CalcType ibarcaconst    = (v(VT_pCa))*0.341000*4.00000/((exp2VFRT) - 1.00000)*dfconst;
  const ML_CalcType ibarca_j       = ibarcaconst*(Ca_j*(exp2VFRT) - Ca_o);
  const ML_CalcType I_Ca_junc      =  (v(VT_Fjunc_CaL))*ibarca_j*((1.00000 - f_Ca_Bj)+(v(VT_fcaCaj)));
  const ML_CalcType powCaj         = (pow(Ca_j, 1.60000));
  const ML_CalcType I_pca_junc     = (Fjunc*((v(VT_powslcap)))*(v(VT_IbarSLCaP))*powCaj)/((v(VT_powKmPCa))+powCaj);
  const ML_CalcType eca_junc       =  0.5/FoRT*log(Ca_o/Ca_j);
  const ML_CalcType I_cabk_junc    =  Fjunc*(v(VT_GCaB))*(svolt - eca_junc);
  const ML_CalcType I_Ca_tot_junc  = (I_Ca_junc+I_cabk_junc+I_pca_junc) -  2.00000*I_ncx_junc;

  const ML_CalcType J_CaB_sl    = dSLL_sl+dSLH_sl;
  const ML_CalcType ibarca_sl   = ibarcaconst*(Ca_sl*(exp2VFRT) - Ca_o);
  const ML_CalcType I_Ca_sl     =  (v(VT_Fsl_CaL))*ibarca_sl*((1.00000 - f_Ca_Bsl)+(v(VT_fcaCaMSL)));
  const ML_CalcType powCasl     = (pow(Ca_sl, 1.60000));
  const ML_CalcType I_pca_sl    = (Fsl*((v(VT_powslcap)))*(v(VT_IbarSLCaP))*powCasl)/((v(VT_powKmPCa))+powCasl);
  const ML_CalcType eca_sl      =  0.5/FoRT*log(Ca_o/Ca_sl);
  const ML_CalcType I_cabk_sl   =  Fsl*(v(VT_GCaB))*(svolt - eca_sl);
  const ML_CalcType I_Ca_tot_sl = (I_Ca_sl+I_cabk_sl+I_pca_sl) -  2.00000*I_ncx_sl;

  const ML_CalcType I_nak     = I_nak_junc+I_nak_sl;
  const ML_CalcType ek        =  (1.00000/FoRT)*(log((K_o/K_i)));
  const ML_CalcType I_kr      =  (v(VT_gkr))*x_kr*ptTeaP->rkr[Vi]*(svolt - ek);
  const ML_CalcType ikpconst  = (v(VT_gkp))*ptTeaP->kp_kp[Vi]*(svolt - ek);
  const ML_CalcType I_kp_junc =  Fjunc*ikpconst;
  const ML_CalcType I_kp_sl   =  Fsl*ikpconst;
  const ML_CalcType I_kp      = I_kp_junc+I_kp_sl;
  const ML_CalcType eks       =  (1.00000/FoRT)*(log(((K_o+ (v(VT_pNaK))*Na_o)/(K_i+ (v(VT_pNaK))*Na_i))));
  const ML_CalcType iksconst  = x_ks*x_ks*(svolt - eks);
  const ML_CalcType I_ks_junc =  Fjunc*(v(VT_gks_junc))*iksconst;
  const ML_CalcType I_ks_sl   =  Fsl*(v(VT_gks_sl))*iksconst;
  const ML_CalcType I_ks      = I_ks_junc+I_ks_sl;
  const ML_CalcType I_tof     =  (v(VT_GtoFast))*x_to_f*y_to_f*(svolt - ek);
  const ML_CalcType I_to      = I_tof;
  const ML_CalcType I_kur     = (v(VT_Gkur))*x_kur*y_kur*(svolt-ek);
  const ML_CalcType aki       = 1.02000/(1.00000+(exp((0.238500*((svolt - ek) - 59.2150)))));
  const ML_CalcType bki       =
    (0.491240*(exp((0.0803200*((svolt+5.47600) - ek))))+(exp((0.0617500*((svolt - ek) - 594.310)))))/
    (1.00000+(exp((-0.514300*((svolt - ek)+4.75300)))));
  const ML_CalcType kiss  = aki/(aki+bki);
  const ML_CalcType I_ki  =  (v(VT_Ikiconst))*kiss*(svolt - ek);
  const ML_CalcType ibark = ( (v(VT_pK))*(0.750000*K_i*(expVFRT) -  0.750000*K_o))/((expVFRT) - 1.00000);
  const ML_CalcType I_CaK =  ibark*dfconst*
    ( (v(VT_Fjunc_CaL))*((v(VT_fcaCaj))+(1.00000 - f_Ca_Bj))+ (v(VT_Fsl_CaL))*((v(VT_fcaCaMSL))+(1.00000 - f_Ca_Bsl)));
  const ML_CalcType I_K_tot        = ((I_to+I_kr+I_ks+I_ki) -  2.00000*I_nak)+I_CaK+I_kp+I_kur;
  const ML_CalcType I_Na_tot       = I_Na_tot_junc+I_Na_tot_sl;
  const ML_CalcType iclcaconst     = (v(VT_GClCa))*(svolt - (v(VT_ecl)));
  const ML_CalcType I_ClCa_junc    =  ((Fjunc*iclcaconst)/(1.00000+(v(VT_KdClCa))/Ca_j));
  const ML_CalcType I_ClCa_sl      =  ((Fsl*iclcaconst)/(1.00000+(v(VT_KdClCa))/Ca_sl));
  const ML_CalcType I_ClCa         = I_ClCa_junc+I_ClCa_sl;
  const ML_CalcType I_Cl_tot       = I_ClCa+ptTeaP->I_Clbk[Vi];
  const ML_CalcType I_Ca_tot       = I_Ca_tot_junc+I_Ca_tot_sl;
  const ML_CalcType I_tot          = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
  const ML_CalcType I_Na           = I_Na_junc+I_Na_sl;
  const ML_CalcType I_nabk         = I_nabk_junc+I_nabk_sl;
  const ML_CalcType I_Ca           = I_Ca_junc+I_Ca_sl;
  const ML_CalcType I_CaNa         = I_CaNa_junc+I_CaNa_sl;
  const ML_CalcType I_Catot        = I_Ca+I_CaK+I_CaNa;
  const ML_CalcType I_Na_tot_junc2 =  3.00000*I_ncx_junc+ 3.00000*I_nak_junc+I_CaNa_junc;
  const ML_CalcType I_ncx          = I_ncx_junc+I_ncx_sl;
  const ML_CalcType I_Na_tot_sl2   =  3.00000*I_ncx_sl+ 3.00000*I_nak_sl+I_CaNa_sl;
  const ML_CalcType I_pca          = I_pca_junc+I_pca_sl;
  const ML_CalcType I_cabk         = I_cabk_junc+I_cabk_sl;

  const ML_CalcType I_mem = I_Na_junc + I_Na_sl + I_Na + I_nabk_junc + I_nabk_sl + I_nabk + I_nak_junc + I_nak_sl +
    I_nak + I_kr + I_kp_junc + I_kp_sl + I_kp + I_ks_junc + I_ks_sl + I_ks + I_kur + I_NaL_sl + I_NaL_junc + I_tof +
    I_to + I_ki + I_ClCa_junc + I_ClCa_sl + I_ClCa + ptTeaP->I_Clbk[Vi] + I_Ca_junc + I_Ca_sl + I_Ca + I_CaK +
    I_CaNa_junc + I_CaNa_sl + I_CaNa + I_Catot + I_ncx_junc + I_ncx_sl + I_Na_tot_junc2 + I_ncx + I_pca_junc +
    I_Na_tot_sl2 + I_pca_sl + I_pca + I_cabk_junc + I_cabk_sl + I_cabk + I_Na_tot_junc + I_Na_tot_sl + I_K_tot +
    I_Ca_tot_junc + I_Ca_tot_sl + I_Na_tot + I_Cl_tot + I_Ca_tot + I_tot;
  tempstr<<' '<<I_Na_junc << ' '<<I_Na_sl << ' '<<I_Na << ' '<<I_nabk_junc << ' '<<I_nabk_sl << ' '<<I_nabk << ' '<<
    I_nak_junc << ' '<<I_nak_sl << ' '
         <<I_nak << ' '<<I_kr << ' '<<I_kp_junc << ' '<<I_kp_sl << ' '<<I_kp << ' '<<I_ks_junc << ' '<<I_ks_sl << ' '<<
    I_ks << ' '<<I_kur << ' '<<I_NaL_sl << ' '<<I_NaL_junc << ' '
         <<I_tof << ' '<<I_to << ' '<<I_ki << ' '<<I_ClCa_junc << ' '<<I_ClCa_sl << ' '<<I_ClCa << ' '<<
    ptTeaP->I_Clbk[Vi] << ' '<<I_Ca_junc << ' '<<I_Ca_sl << ' '
         <<I_Ca << ' '<<I_CaK << ' '<<I_CaNa_junc << ' '<<I_CaNa_sl << ' '<<I_CaNa << ' '<<I_Catot << ' '<<I_ncx_junc <<
    ' '<<I_ncx_sl << ' '<<I_Na_tot_junc2 << ' '
         <<I_ncx << ' '<<I_pca_junc << ' '<<I_Na_tot_sl2 << ' '<<I_pca_sl << ' '<<I_pca << ' '<<I_cabk_junc << ' '<<
    I_cabk_sl << ' '<<I_cabk << ' '<<I_Na_tot_junc << ' '
         <<I_Na_tot_sl << ' '<<I_K_tot << ' '<<I_Ca_tot_junc << ' '<<I_Ca_tot_sl << ' '<<I_Na_tot << ' '<<I_Cl_tot <<
    ' '<<I_Ca_tot << ' '<< J_SRCarel << ' '<< J_serca << ' '<< J_SRleak << ' '<< J_CaB_cytosol << ' '<<
    J_CaB_junction <<
    ' '<< J_CaB_sl << ' '<<I_tot << ' '<< I_mem << ' ';
}  // GrandiEtAlAtrium::LongPrint

void GrandiEtAlAtrium::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 40;
  const string ParaNames[numpara] =
  {"Na_j",   "Na_sl",     "K_i",       "Ca_j",       "Ca_sl",            "m",                 "h",                 "j",
   "x_kr",
   "x_ks",
   "Na_i",   "ml",
   "hl",     "x_kur",     "y_kur",
   "x_to_f", "y_to_f",    "d",         "f",          "f_Ca_Bj",          "f_Ca_Bsl",          "Ry_Rr",
   "Ry_Ro",
   "Ry_Ri",
   "Ca_sr",
   "Ca_i",   "Na_Bj",
   "Na_Bsl",
   "Tn_CL",  "Tn_CHc",    "Tn_CHm",    "CaM",        "Myo_c",            "Myo_m",             "SRB",
   "SLL_j",
   "SLL_sl",
   "SLH_j",
   "SLH_sl",
   "Csqn_b"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void GrandiEtAlAtrium::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 61;
  const string ParaNames[numpara] =
  {"I_Na_junc",      "I_Na_sl",               "I_Na",                        "I_nabk_junc",
   "I_nabk_sl",
   "I_nabk",
   "I_nak_junc",
   "I_nak_sl",       "I_nak",                 "I_kr",
   "I_kp_junc",      "I_kp_sl",               "I_kp",                        "I_ks_junc",
   "I_ks_sl",
   "I_ks",
   "I_kur",
   "I_NaL_sl",       "I_NaL_junc",            "I_tof",                       "I_to",
   "I_ki",           "I_ClCa_junc",           "I_ClCa_sl",                   "I_ClCa",
   "I_Clbk",
   "I_Ca_junc",
   "I_Ca_sl",
   "I_Ca",           "I_CaK",                 "I_CaNa_junc",
   "I_CaNa_sl",      "I_CaNa",                "I_Catot",                     "I_ncx_junc",
   "I_ncx_sl",
   "I_Na_tot_junc2", "I_ncx",
   "I_pca_junc",     "I_Na_tot_sl2",
   "I_pca_sl",       "I_pca",                 "I_cabk_junc",                 "I_cabk_sl",
   "I_cabk",
   "I_Na_tot_junc",  "I_Na_tot_sl",
   "I_K_tot",
   "I_Ca_tot_junc",  "I_Ca_tot_sl",           "I_Na_tot",                    "I_Cl_tot",
   "I_Ca_tot",
   "J_SRCarel",
   "J_serca",
   "J_SRleak",
   "J_CaB_cytosol",  "J_CaB_junction",        "J_CaB_sl",                    "I_tot",
   "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
} // GrandiEtAlAtrium::GetLongParameterNames
