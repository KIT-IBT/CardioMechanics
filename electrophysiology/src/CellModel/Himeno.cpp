/*
 * File: Himeno.cpp
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


#include <Himeno.h>


Himeno::Himeno(HimenoParameters *pp) {
  pHimP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(pHimP, NS_HimenoParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

Himeno::~Himeno() {}

#ifdef HETERO

inline bool Himeno::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool Himeno::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

/*! \fn int Himeno::GetSize(void)
 *  \brief Return the size of the memory that has to be included into the backup.
 *  \return Size of memory to be used for backup.
 */
inline int Himeno::GetSize(void) {
  return sizeof(Himeno)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(HimenoParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline ML_CalcType *Himeno::GetBase(void) {
  return &CaMCa;
}

inline unsigned char Himeno::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void Himeno::Init() {
#if KADEBUG
  cerr << "#initializing Class: Himeno ... " << endl;
#endif  // if KADEBUG
  CaMCa         = v(VT_init_CaMCa);
  SRCa          = v(VT_init_SRCa);
  TnChCa        = v(VT_init_TnChCa);
  p_O_NaT       = v(VT_init_p_O_NaT);
  p_I_2_NaT     = v(VT_init_p_I_2_NaT);
  p_I_s_NaT     = v(VT_init_p_I_s_NaT);
  p_O_NaL       = v(VT_init_p_O_NaL);
  p_I_1_NaL     = v(VT_init_p_I_1_NaL);
  p_I_2_NaL     = v(VT_init_p_I_2_NaL);
  p_I_s_NaL     = v(VT_init_p_I_s_NaL);
  chi_r_fast    = v(VT_init_chi_r_fast);
  chi_r_slow    = v(VT_init_chi_r_slow);
  para_Xs1      = v(VT_init_para_Xs1);
  para_Xs2      = v(VT_init_para_Xs2);
  i_fast        = v(VT_init_i_fast);
  i_slow        = v(VT_init_i_slow);
  P_7           = v(VT_init_P_7);
  P_8_13        = v(VT_init_P_8_13);
  P_1_6         = v(VT_init_P_1_6);
  p_E_1_NCX_blk = v(VT_init_p_E_1_NCX_blk);
  p_I_2_NCX_blk = v(VT_init_p_I_2_NCX_blk);
  p_I_1_NCX_blk = v(VT_init_p_I_1_NCX_blk);
  p_E_1_NCX_iz  = v(VT_init_p_E_1_NCX_iz);
  p_I_1_NCX_iz  = v(VT_init_p_I_1_NCX_iz);
  p_I_2_NCX_iz  = v(VT_init_p_I_2_NCX_iz);
  Y_ooo         = v(VT_init_Y_ooo);
  Y_ooc         = v(VT_init_Y_ooc);
  Y_occ         = v(VT_init_Y_occ);
  Y_coc         = v(VT_init_Y_coc);
  Y_coo         = v(VT_init_Y_coo);
  Y_cco         = v(VT_init_Y_cco);
  Y_oco         = v(VT_init_Y_oco);
  Y_co_iz       = v(VT_init_Y_co_iz);
  Y_oo_iz       = v(VT_init_Y_oo_iz);
  Y_oc_iz       = v(VT_init_Y_oc_iz);
  Y_co_blk      = v(VT_init_Y_co_blk);
  Y_oo_blk      = v(VT_init_Y_oo_blk);
  Y_oc_blk      = v(VT_init_Y_oc_blk);
  Ca_2_tot_jnc  = v(VT_init_Ca_2_tot_jnc);
  Ca_2_tot_iz   = v(VT_init_Ca_2_tot_iz);
  Ca_2_tot_blk  = v(VT_init_Ca_2_tot_blk);
  Ca_2_SRup     = v(VT_init_Ca_2_SRup);
  Ca_2_tot_SRrl = v(VT_init_Ca_2_tot_SRrl);
  Nai           = v(VT_init_Nai);
  Ki            = v(VT_init_Ki);
  TS_S          = v(VT_init_TS_S);
  TS_W          = v(VT_init_TS_W);
  TSCa_3        = v(VT_init_TSCa_3);
  TSCa_3W       = v(VT_init_TSCa_3W);
  TSCa_3S       = v(VT_init_TSCa_3S);
  hw            = v(VT_init_hw);
  hp            = v(VT_init_hp);
  Pb_spm        = v(VT_init_Pb_spm);
  a             = v(VT_init_a);
  L_bound_iz    = v(VT_init_L_bound_iz);
  H_bound_iz    = v(VT_init_H_bound_iz);
  Ca_2_iz       = 0;
  Ca_2_blk      = 0;
  L_free_jnc    = 0;
  H_free_jnc    = 0;
  Ca_2_nd_00    = 0;
  Ca_2_nd_L0    = 0;
  Ca_2_nd_0R    = 0;
  Ca_2_nd_LR    = 0;
}  // Himeno::Init

ML_CalcType Himeno::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                         int euler                                            = 2) {
  ML_CalcType Vm = V*1000;  // membrane voltage in mV
  const int   Vi = (int)(DivisionTab*(RangeTabhalf+Vm)+.5); // array position


  // TenTusscher2 legacy
  // convert from nA/cell to pA/pF with Cm: 2 mueF/cm^2=2e-2 F/m^2, S=0.2/mue m=0.2e6/m -fs
  //  i_external=i_external/(2e-2*Volume()*0.2e6*1e9);
  // hard-coded, because S is currently questionable in TenTusscher2.h
  // i_external=i_external/4.59312e-2; // convert from nA to pA/pF


  const ML_CalcType RTONF = v(VT_RTONF);
  const ML_CalcType Cao   = v(VT_Cao);
  const ML_CalcType Ko    = v(VT_Ko);
  const ML_CalcType Nao   = v(VT_Nao);

  // bulkSpace()
  const ML_CalcType dCaMCadt  = v(VT_k_on_CaM) * Ca_2_blk * (v(VT_B_tot_CaM) - CaMCa) - v(VT_k_off_CaM) * CaMCa;
  const ML_CalcType dTnChCadt = v(VT_k_on_TnCh) * Ca_2_blk * (v(VT_B_tot_TnCh) - TnChCa) - v(VT_k_off_TnCh) * TnChCa;
  const ML_CalcType dSRCadt   = v(VT_k_on_SR) * Ca_2_blk * (v(VT_B_tot_SR) - SRCa) - v(VT_k_off_SR) * SRCa;

  Ca_2_blk = Ca_2_tot_blk - (CaMCa + TnChCa + SRCa + 3 * (TSCa_3 + TSCa_3W + TSCa_3S) / 1000);

  // intermediateZone()
  ML_CalcType L_free_iz = v(VT_B_tot_L_iz) - L_bound_iz;
  ML_CalcType H_free_iz = v(VT_B_tot_H_iz) - H_bound_iz;
  for (int n = 0; n < 10; n++) {
    Ca_2_iz   = Ca_2_tot_iz / (1 + L_free_iz / v(VT_K_dL_iz) + H_free_iz / v(VT_K_dH_iz));
    L_free_iz = v(VT_B_tot_L_iz) / (1 + Ca_2_iz / v(VT_K_dL_iz));
    H_free_iz = v(VT_B_tot_H_iz) / (1 + Ca_2_iz / v(VT_K_dH_iz));
  }
  L_bound_iz = v(VT_B_tot_L_iz) - L_free_iz;
  H_bound_iz = v(VT_B_tot_H_iz) - H_free_iz;

  // junctionalSpace()
  for (int n = 0; n < 10; n++) {
    Ca_2_jnc   = Ca_2_tot_jnc / (1 + L_free_jnc / v(VT_K_dL_jnc) + H_free_jnc / v(VT_K_dH_jnc));
    L_free_jnc = v(VT_B_tot_L_jnc) / (1 + Ca_2_jnc / v(VT_K_dL_jnc));
    H_free_jnc = v(VT_B_tot_H_jnc) / (1 + Ca_2_jnc / v(VT_K_dH_jnc));
  }

  // releaseSiteOfSR()
  const ML_CalcType a_ = 1;
  const ML_CalcType b  = v(VT_B_tot_CSQN) - Ca_2_tot_SRrl + v(VT_K_d_CSQN_Ca);
  const ML_CalcType c  = -v(VT_K_d_CSQN_Ca) * Ca_2_tot_SRrl;
  Ca_2_SRrl = (-b + sqrt(pow(b, 2) - 4 * a_ * c)) / (2 * a_);

  // boundaryDiffusion()
  const ML_CalcType J_Ca_jnciz = v(VT_G_dCa_jnciz) * (Ca_2_jnc - Ca_2_iz) * v(VT_Sc_Cell);
  const ML_CalcType J_Ca_izblk = v(VT_G_dCa_izblk) * (Ca_2_iz - Ca_2_blk) * v(VT_Sc_Cell);
  const ML_CalcType J_trans_SR = v(VT_P_trans) * (Ca_2_SRup - Ca_2_SRrl) * v(VT_Sc_Cell);

  // currentCaL()
  const ML_CalcType epsilon_plus_iz  = (Ca_2_iz * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType epsilon_plus_blk = (Ca_2_blk * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType expdRTFVm        = exp(-2 * v(VT_F)/(v(VT_R)*v(VT_Tx))*Vm);
  const ML_CalcType Ca_2_iz_loc      =
    (Ca_2_iz + v(VT_f_L) * 2* v(VT_inverseRTONF) * Vm * expdRTFVm / (1 - expdRTFVm) * Cao) /
    (1 + v(VT_f_L) * 2* v(VT_inverseRTONF) * Vm / (1 - expdRTFVm));
  const ML_CalcType Ca_2_blk_loc =
    (Ca_2_blk + v(VT_f_L) * (2* v(VT_inverseRTONF) * Vm * expdRTFVm) / (1 - expdRTFVm) * Cao) /
    (1 + v(VT_f_L) * 2* v(VT_inverseRTONF) * Vm / (1 - expdRTFVm));
  const ML_CalcType epsilon_plus_iz_loc  = (Ca_2_iz_loc * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType epsilon_plus_blk_loc = (Ca_2_blk_loc * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType p_O_LCC              = Y_ooo + Y_ooc;
  const ML_CalcType E_Ca_jnc             = RTONF / 2 * log(Cao / Ca_2_nd_L0);
  const ML_CalcType E_Ca_blk             = RTONF / 2 * log(Cao / Ca_2_blk);
  const ML_CalcType E_Ca_iz              = RTONF / 2 * log(Cao / Ca_2_iz);
  const ML_CalcType E_K                  = RTONF / 1 * log(Ko / Ki);
  const ML_CalcType E_Na                 = RTONF / 1 * log(v(VT_Nao) / Nai);
  const ML_CalcType exp_VdRTF            = exp(-Vm * v(VT_F)/(v(VT_R)*v(VT_Tx)));
  const ML_CalcType GHK_Ca_LR            = 2 * Vm / RTONF * (Ca_2_nd_LR - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_L0            = 2 * Vm / RTONF * (Ca_2_nd_L0 - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_iz            = 2 * Vm / RTONF * (Ca_2_iz - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_blk           = 2 * Vm / RTONF * (Ca_2_blk - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Na               = 1 * Vm / RTONF * (Nai - Nao * exp_VdRTF) / (1 - exp_VdRTF);
  const ML_CalcType GHK_K                = 1 * Vm / RTONF * (Ki - Ko * exp_VdRTF) / (1 - exp_VdRTF);
  const ML_CalcType I_CaL_Ca_blk         = v(VT_f_CaL_blk) * v(VT_P_CaL_Ca) * GHK_Ca_blk * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_iz          = v(VT_f_CaL_iz) * v(VT_P_CaL_Ca) * GHK_Ca_iz * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_LR          = v(VT_f_CaL_jnc) * v(VT_P_CaL_Ca) * GHK_Ca_LR * Y_ooo * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_L0          = v(VT_f_CaL_jnc) * v(VT_P_CaL_Ca) * GHK_Ca_L0 * Y_ooc * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_blk         = v(VT_f_CaL_blk) * v(VT_P_CaL_Na) * GHK_Na * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_iz          = v(VT_f_CaL_iz) * v(VT_P_CaL_Na) * GHK_Na * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_jnc         = v(VT_f_CaL_jnc) * v(VT_P_CaL_Na) * GHK_Na * p_O_LCC * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_blk          = v(VT_f_CaL_blk) * v(VT_P_CaL_K) * GHK_K * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_iz           = v(VT_f_CaL_iz) * v(VT_P_CaL_K) * GHK_K * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_jnc          = v(VT_f_CaL_jnc) * v(VT_P_CaL_K) * GHK_K * p_O_LCC * v(VT_ATPfactor);
  const ML_CalcType I_CaL                = (I_CaL_Ca_LR + I_CaL_Ca_L0 + I_CaL_Na_jnc + I_CaL_K_jnc) +
    (I_CaL_Ca_iz + I_CaL_Na_iz + I_CaL_K_iz) + (I_CaL_Ca_blk + I_CaL_Na_blk + I_CaL_K_blk);
  const ML_CalcType Y_cc_iz    = 1 - (Y_co_iz + Y_oo_iz + Y_oc_iz);
  const ML_CalcType dY_co_izdt = pHimP->epsilon_minus[Vi] * Y_cc_iz + pHimP->alpha_minus[Vi] * Y_oo_iz -
    (epsilon_plus_iz + pHimP->alpha_plus[Vi]) * Y_co_iz;
  const ML_CalcType dY_oo_izdt = pHimP->alpha_plus[Vi] * Y_co_iz + pHimP->epsilon_minus[Vi] * Y_oc_iz -
    (pHimP->alpha_minus[Vi] + epsilon_plus_iz_loc) * Y_oo_iz;
  const ML_CalcType dY_oc_izdt = epsilon_plus_iz_loc * Y_oo_iz + pHimP->alpha_plus[Vi] * Y_cc_iz -
    (pHimP->epsilon_minus[Vi] + pHimP->alpha_minus[Vi]) * Y_oc_iz;

  // cerr << dY_oc_izdt << " " << epsilon_plus_iz_loc_ << " " << Y_oo_iz << " " << pHimP->alpha_plus[Vi] << " (" <<
  // alpha_plus_ << ") " << Y_cc_iz << " " << pHimP->epsilon_minus[Vi] << " " << pHimP->alpha_minus[Vi] << " " <<
  // Y_oc_iz << endl;
  const ML_CalcType Y_cc_blk    = 1 - (Y_co_blk + Y_oo_blk + Y_oc_blk);
  const ML_CalcType dY_co_blkdt = pHimP->epsilon_minus[Vi] * Y_cc_blk + pHimP->alpha_minus[Vi] * Y_oo_blk -
    (epsilon_plus_blk + pHimP->alpha_plus[Vi]) * Y_co_blk;
  const ML_CalcType dY_oo_blkdt = pHimP->alpha_plus[Vi] * Y_co_blk + pHimP->epsilon_minus[Vi] * Y_oc_blk -
    (pHimP->alpha_minus[Vi] + epsilon_plus_blk_loc) * Y_oo_blk;
  const ML_CalcType dY_oc_blkdt = epsilon_plus_blk_loc * Y_oo_blk + pHimP->alpha_plus[Vi] * Y_cc_blk -
    (pHimP->epsilon_minus[Vi] + pHimP->alpha_minus[Vi]) * Y_oc_blk;

  // currentNa()
  const ML_CalcType k_OI1      = pHimP->k_OI2[Vi];
  const ML_CalcType k_I1C      = pHimP->k_I2C[Vi];
  const ML_CalcType k_C2I1     = pHimP->k_C2I2[Vi];
  const ML_CalcType p_C_NaT    = 1.0 - p_I_2_NaT - p_I_s_NaT - p_O_NaT;
  const ML_CalcType dp_O_NaTdt = v(VT_k_I2O) * p_I_2_NaT + pHimP->f_C_Na[Vi] * pHimP->k_C2O[Vi] * p_C_NaT -
    (pHimP->k_OC[Vi] + pHimP->k_OI2[Vi]) * p_O_NaT;
  const ML_CalcType dp_I_2_NaTdt = pHimP->f_C_Na[Vi] * pHimP->k_C2I2[Vi] * p_C_NaT + pHimP->k_OI2[Vi] * p_O_NaT +
    pHimP->k_Isb[Vi] * p_I_s_NaT - (pHimP->k_I2C[Vi] + v(VT_k_I2O) + pHimP->k_Isf[Vi]) * p_I_2_NaT;
  const ML_CalcType dp_I_s_NaTdt = pHimP->k_Isf[Vi] * p_I_2_NaT + pHimP->k_Isf[Vi] * p_C_NaT - 2 * pHimP->k_Isb[Vi] *
    p_I_s_NaT;
  const ML_CalcType p_C_NaL    = 1.0 - p_I_2_NaL - p_I_s_NaL - p_I_1_NaL - p_O_NaL;
  const ML_CalcType dp_O_NaLdt = v(VT_k_I1O) * p_I_1_NaL + pHimP->f_C_Na[Vi] * pHimP->k_C2O[Vi] * p_C_NaL -
    (pHimP->k_OC[Vi] + k_OI1) * p_O_NaL;
  const ML_CalcType dp_I_1_NaLdt = k_OI1 * p_O_NaL + pHimP->f_C_Na[Vi] * k_C2I1 * p_C_NaL -
    (v(VT_k_I1O) + k_I1C + v(VT_k_I1I2)) * p_I_1_NaL;
  const ML_CalcType dp_I_2_NaLdt = pHimP->f_C_Na[Vi] * pHimP->k_C2I2[Vi] * p_C_NaL + v(VT_k_I1I2) * p_I_1_NaL +
    pHimP->k_Isb[Vi] * p_I_s_NaL - (pHimP->k_I2C[Vi] + pHimP->k_Isf[Vi]) * p_I_2_NaL;
  const ML_CalcType dp_I_s_NaLdt = pHimP->k_Isf[Vi] * p_I_2_NaL + pHimP->k_Isf[Vi] * p_C_NaL - 2 * pHimP->k_Isb[Vi] *
    p_I_s_NaL;
  const ML_CalcType I_NaT_Na = (1 - v(VT_f_LSM)) * v(VT_P_Na) * GHK_Na * p_O_NaT;
  const ML_CalcType I_NaT_K  = (1 - v(VT_f_LSM)) * v(VT_P_Na) * 0.1 * GHK_K * p_O_NaT;
  const ML_CalcType I_NaT    = I_NaT_Na + I_NaT_K;
  const ML_CalcType I_NaL_Na = v(VT_f_LSM) * v(VT_P_Na) * GHK_Na * p_O_NaL;
  const ML_CalcType I_NaL_K  = v(VT_f_LSM) * v(VT_P_Na) * 0.1 * GHK_K * p_O_NaL;
  const ML_CalcType I_NaL    = I_NaL_Na + I_NaL_K;
  const ML_CalcType I_Na     = I_NaT + I_NaL;

  // currentK1()
  const ML_CalcType alpha_Mg  = 12.0 * exp(-0.025 * (Vm - E_K));
  const ML_CalcType beta_Mg   = 28 * v(VT_Mg_2_cyt) * exp(0.025 * (Vm - E_K));
  const ML_CalcType f_O       = alpha_Mg / (alpha_Mg + beta_Mg);
  const ML_CalcType f_B       = beta_Mg / (alpha_Mg + beta_Mg);
  const ML_CalcType po_Mg     = f_O * f_O * f_O;
  const ML_CalcType po_Mg1    = 3.0 * f_O * f_O * f_B;
  const ML_CalcType po_Mg2    = 3.0 * f_O * f_B * f_B;
  const ML_CalcType alpha_SPM = 0.17 * exp(-0.07 * (Vm - E_K + 8 * v(VT_Mg_2_cyt))) /
    (1.0 + 0.01 * exp(0.12 * (Vm - E_K + 8 * v(VT_Mg_2_cyt))));
  const ML_CalcType beta_SPM = 0.28 * v(VT_SPM) * exp(0.15 * (Vm - E_K + 8 * v(VT_Mg_2_cyt))) /
    (1.0 + 0.01 * exp(0.13 * (Vm - E_K + 8 * v(VT_Mg_2_cyt))));

  const ML_CalcType dPb_spmdt = beta_SPM * po_Mg * (1 - Pb_spm) - alpha_SPM * Pb_spm;
  const ML_CalcType po_mode1  = v(VT_f_mode1) * (1 - Pb_spm) * (po_Mg + (2.0 / 3.0) * po_Mg1 + (1.0 / 3.0) * po_Mg2);
  const ML_CalcType po_mode2  = (1 - v(VT_f_mode1)) / (1.0 + v(VT_SPM) / (40.0 * exp(-(Vm - E_K) / 9.1)));
  const ML_CalcType p_O_K1    = po_mode1 + po_mode2;
  const ML_CalcType I_K1      = v(VT_G_K1) * v(VT_chi_K1) * (Vm - E_K) * p_O_K1;

  // currentKr()
  const ML_CalcType dchi_r_fastdt = (pHimP->chi_r_infinity[Vi] - chi_r_fast) / pHimP->tau_chi_r_fast[Vi];
  const ML_CalcType dchi_r_slowdt = (pHimP->chi_r_infinity[Vi] - chi_r_slow) / pHimP->tau_chi_r_slow[Vi];
  const ML_CalcType chi_r         = pHimP->A_chi_r_fast[Vi] * chi_r_fast + pHimP->A_chi_r_slow[Vi] * chi_r_slow;
  const ML_CalcType p_O_Kr        = chi_r * pHimP->R_Kr[Vi];
  const ML_CalcType I_Kr          = v(VT_G_Kr) * v(VT_chi_Kr) * (Vm - E_K) * p_O_Kr;

  // currentsKs()
  const ML_CalcType dpara_Xs1dt       = (pHimP->para_Xs1_infinity[Vi] - para_Xs1) / pHimP->tau_Xs1[Vi];
  const ML_CalcType para_Xs2_infinity = pHimP->para_Xs1_infinity[Vi];
  const ML_CalcType dpara_Xs2dt       = (para_Xs2_infinity - para_Xs2) / pHimP->tau_Xs2[Vi];
  const ML_CalcType para_RKs_blk      = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_blk, 1.4));
  const ML_CalcType para_RKs_iz       = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_iz, 1.4));
  const ML_CalcType p_O_Ks_blk        = para_Xs1 * para_Xs2 * para_RKs_blk;
  const ML_CalcType p_O_Ks_iz         = para_Xs1 * para_Xs2 * para_RKs_iz;
  const ML_CalcType I_Ks_K_blk        = v(VT_f_Ks_blk) * v(VT_P_Ks_K) * GHK_K * p_O_Ks_blk;
  const ML_CalcType I_Ks_K_iz         = v(VT_f_Ks_iz) * v(VT_P_Ks_K) * GHK_K * p_O_Ks_iz;
  const ML_CalcType I_Ks_Na_blk       = v(VT_f_Ks_blk) * v(VT_P_Ks_Na) * GHK_Na * p_O_Ks_blk;
  const ML_CalcType I_Ks_Na_iz        = v(VT_f_Ks_iz) * v(VT_P_Ks_Na) * GHK_Na * p_O_Ks_iz;
  const ML_CalcType I_Ks              = I_Ks_K_blk + I_Ks_K_iz + I_Ks_Na_blk + I_Ks_Na_iz;

  // currentKto()
  // const ML_CalcType tau_i_fast = (4.562 + 1 / (0.3933 * exp(-(Vm + 100) / 100) + 0.08004 * exp((Vm + 50) /
  // 16.59)))*(1.0-(0.95)/(1.0+exp((Vm +70.0)/5.0)));
  // const ML_CalcType tau_i_slow = (23.62 + 1 / (0.001416 * exp(-(Vm + 96.52) / 59.05) + 0.000000017808 * exp((Vm +
  // 114.1) / 8.079)))*(1.0-(0.95)/(1.0+exp((Vm +70.0)/5.0)));
  const ML_CalcType dadt      = (pHimP->a_infinity[Vi] - a) / pHimP->tau_a[Vi];
  const ML_CalcType di_fastdt = (pHimP->i_infinity[Vi] - i_fast) / pHimP->tau_i_fast[Vi];
  const ML_CalcType di_slowdt = (pHimP->i_infinity[Vi] - i_slow) / pHimP->tau_i_slow[Vi];

  // const ML_CalcType di_fastdt = (pHimP->i_infinity[Vi] - i_fast) / tau_i_fast;
  // const ML_CalcType di_slowdt = (pHimP->i_infinity[Vi] - i_slow) / tau_i_slow;
  const ML_CalcType i       = pHimP->A_i_fast[Vi] * i_fast + pHimP->A_i_slow[Vi] * i_slow;
  const ML_CalcType p_O_Kto = a * i;

  const ML_CalcType I_Kto = v(VT_G_Kto) * p_O_Kto * (Vm - E_K);


  // currentKpl()
  const ML_CalcType I_Kpl = v(VT_P_Kpl) * v(VT_chi_Kpl) * pHimP->p_O_Kpl[Vi] * GHK_K;

  // currentCab()
  const ML_CalcType I_Cab_blk = v(VT_P_Cab) * v(VT_f_Cab_blk) * GHK_Ca_blk;
  const ML_CalcType I_Cab_iz  = v(VT_P_Cab) * v(VT_f_Cab_iz) * GHK_Ca_iz;
  const ML_CalcType I_Cab     = I_Cab_iz + I_Cab_blk;

  // currentbNSC()
  const ML_CalcType I_bNSC_K  = v(VT_P_bNSC_K) * GHK_K;
  const ML_CalcType I_bNSC_Na = v(VT_P_bNSC_Na) * GHK_Na;
  const ML_CalcType I_bNSC    = I_bNSC_K + I_bNSC_Na;

  // currentlCa()
  const ML_CalcType p_O_blk       = 1.0 / (1.0 + pow(0.0012 / Ca_2_blk, 3));
  const ML_CalcType p_O_iz        = 1.0 / (1.0 + pow(0.0012 / Ca_2_iz, 3));
  const ML_CalcType I_l_Ca_Na_blk = v(VT_P_l_Ca_Na) * v(VT_f_l_Ca_blk) * GHK_Na * p_O_blk;
  const ML_CalcType I_l_Ca_Na_iz  = v(VT_P_l_Ca_Na) * v(VT_f_l_Ca_iz) * GHK_Na * p_O_iz;
  const ML_CalcType I_l_Ca_K_blk  = v(VT_P_l_Ca_K) * v(VT_f_l_Ca_blk) * GHK_K * p_O_blk;
  const ML_CalcType I_l_Ca_K_iz   = v(VT_P_l_Ca_K) * v(VT_f_l_Ca_iz) * GHK_K * p_O_iz;
  const ML_CalcType I_l_Ca        = I_l_Ca_Na_iz + I_l_Ca_K_iz + I_l_Ca_Na_blk + I_l_Ca_K_blk;

  // currentKATP()
  const ML_CalcType I_KATP = v(VT_G_KATP) * (Vm - E_K) * v(VT_p_O_KATP) * v(VT_chi_KATP);

  // currentNaK()
  const ML_CalcType Nai_bar       = Nai / pHimP->K_d_Nai[Vi];
  const ML_CalcType Ki_bar        = Ki / pHimP->K_d_Ki[Vi];
  const ML_CalcType alpha_1_plus  = v(VT_k_1_plus) * pow(Nai_bar, 3) / (pow(1 + Nai_bar, 3) + pow(1 + Ki_bar, 2) - 1);
  const ML_CalcType alpha_4_minus = v(VT_k_4_minus) * pow(Ki_bar, 2) / (pow(1 + Nai_bar, 3) + pow(1 + Ki_bar, 2) - 1);
  const ML_CalcType P_14_15       = 1 - P_1_6 - P_7 - P_8_13;
  const ML_CalcType V_step1       = alpha_1_plus * P_1_6 - v(VT_alpha_1_minus) * P_7;
  const ML_CalcType V_step2       = v(VT_alpha_2_plus) * P_7 - pHimP->alpha_2_minus[Vi] * P_8_13;
  const ML_CalcType V_step3       = pHimP->alpha_3_plus[Vi] * P_8_13 - v(VT_alpha_3_minus) * P_14_15;
  const ML_CalcType V_step4       = v(VT_alpha_4_plus) * P_14_15 - alpha_4_minus * P_1_6;
  const ML_CalcType v_cyc_NaK     = V_step2;
  const ML_CalcType dP_1_6dt      = -alpha_1_plus * P_1_6 + v(VT_alpha_1_minus) * P_7 + v(VT_alpha_4_plus) * P_14_15 -
    alpha_4_minus * P_1_6;
  const ML_CalcType dP_7dt = -v(VT_alpha_2_plus) * P_7 + pHimP->alpha_2_minus[Vi] * P_8_13 + alpha_1_plus *
    P_1_6 - v(VT_alpha_1_minus) * P_7;
  const ML_CalcType dP_8_13dt = -pHimP->alpha_3_plus[Vi] * P_8_13 + v(VT_alpha_3_minus) * P_14_15 + v(VT_alpha_2_plus) *
    P_7 - pHimP->alpha_2_minus[Vi] * P_8_13;
  const ML_CalcType I_NaK    = v(VT_Amp_NaK) * v_cyc_NaK;
  const ML_CalcType I_NaK_Na = 3 * I_NaK;
  const ML_CalcType I_NaK_K  = -2 * I_NaK;

  // currentNCX()
  const ML_CalcType f_Caina_blk  = Ca_2_blk / (Ca_2_blk + v(VT_K_m_act));
  const ML_CalcType f_Caina_iz   = Ca_2_iz / (Ca_2_iz + v(VT_K_m_act));
  const ML_CalcType q_blk_E_1_Na = 1.0 / (1.0 + pow(v(VT_K_m_Nai) / Nai, 3) * (1.0 + Ca_2_blk / v(VT_K_m_Cai)));
  const ML_CalcType q_iz_E_1_Na  = 1.0 / (1.0 + pow(v(VT_K_m_Nai) / Nai, 3) * (1.0 + Ca_2_iz / v(VT_K_m_Cai)));
  const ML_CalcType q_blk_E_1_Ca = 1.0 / (1.0 + v(VT_K_m_Cai) / Ca_2_blk * (1.0 + pow(Nai / v(VT_K_m_Nai), 3)));
  const ML_CalcType q_iz_E_1_Ca  = 1.0 / (1.0 + v(VT_K_m_Cai) / Ca_2_iz * (1.0 + pow(Nai / v(VT_K_m_Nai), 3)));
  const ML_CalcType alpha_1_blk  = q_blk_E_1_Na *
    (f_Caina_blk * v(VT_alpha_1_on) + (1 - f_Caina_blk) * v(VT_alpha_1_off));
  const ML_CalcType alpha_1_iz = q_iz_E_1_Na *
    (f_Caina_iz * v(VT_alpha_1_on) + (1 - f_Caina_iz) * v(VT_alpha_1_off));
  const ML_CalcType beta_1_blk    = f_Caina_blk * v(VT_beta_1_on) + (1 - f_Caina_blk) * v(VT_beta_1_off);
  const ML_CalcType beta_1_iz     = f_Caina_iz * v(VT_beta_1_on) + (1 - f_Caina_iz) * v(VT_beta_1_off);
  const ML_CalcType alpha_2_blk   = f_Caina_blk * v(VT_alpha_2_on) + (1 - f_Caina_blk) * v(VT_alpha_2_off);
  const ML_CalcType alpha_2_iz    = f_Caina_iz * v(VT_alpha_2_on) + (1 - f_Caina_iz) * v(VT_alpha_2_off);
  const ML_CalcType beta_2_blk    = f_Caina_blk * v(VT_beta_2_on) + (1 - f_Caina_blk) * v(VT_beta_2_off);
  const ML_CalcType beta_2_iz     = f_Caina_iz * v(VT_beta_2_on) + (1 - f_Caina_iz) * v(VT_beta_2_off);
  const ML_CalcType beta_E_blk    = pHimP->k_1[Vi] * q_blk_E_1_Na + v(VT_k_3) * q_blk_E_1_Ca;
  const ML_CalcType beta_E_iz     = pHimP->k_1[Vi] * q_iz_E_1_Na + v(VT_k_3) * q_iz_E_1_Ca;
  const ML_CalcType p_E_2_NCX_blk = 1 - p_E_1_NCX_blk - p_I_1_NCX_blk - p_I_2_NCX_blk;
  const ML_CalcType p_E_2_NCX_iz  = 1 - p_E_1_NCX_iz - p_I_1_NCX_iz - p_I_2_NCX_iz;
  const ML_CalcType v_cyc_NCX_blk = pHimP->k_1[Vi] * q_blk_E_1_Na * p_E_1_NCX_blk - pHimP->k_2[Vi] * v(VT_q_E_2_Na) *
    p_E_2_NCX_blk;
  const ML_CalcType v_cyc_NCX_iz = pHimP->k_1[Vi] * q_iz_E_1_Na * p_E_1_NCX_iz - pHimP->k_2[Vi] * v(VT_q_E_2_Na) *
    p_E_2_NCX_iz;
  const ML_CalcType dp_E_1_NCX_blkdt = p_E_2_NCX_blk * pHimP->alpha_E[Vi] + p_I_1_NCX_blk * beta_1_blk + p_I_2_NCX_blk *
    beta_2_blk - p_E_1_NCX_blk * (beta_E_blk + alpha_1_blk + alpha_2_blk);
  const ML_CalcType dp_I_1_NCX_blkdt = p_E_1_NCX_blk * alpha_1_blk - p_I_1_NCX_blk * beta_1_blk;
  const ML_CalcType dp_I_2_NCX_blkdt = p_E_1_NCX_blk * alpha_2_blk - p_I_2_NCX_blk * beta_2_blk;
  const ML_CalcType dp_E_1_NCX_izdt  = p_E_2_NCX_iz * pHimP->alpha_E[Vi] + p_I_1_NCX_iz * beta_1_iz + p_I_2_NCX_iz *
    beta_2_iz - p_E_1_NCX_iz * (beta_E_iz + alpha_1_iz + alpha_2_iz);
  const ML_CalcType dp_I_1_NCX_izdt = p_E_1_NCX_iz * alpha_1_iz - p_I_1_NCX_iz * beta_1_iz;
  const ML_CalcType dp_I_2_NCX_izdt = p_E_1_NCX_iz * alpha_2_iz - p_I_2_NCX_iz * beta_2_iz;
  const ML_CalcType I_NCX_blk       = v(VT_f_NCX_blk) * v(VT_Amp_NCX) * v_cyc_NCX_blk;
  const ML_CalcType I_NCX_iz        = v(VT_f_NCX_iz) * v(VT_Amp_NCX) * v_cyc_NCX_iz;
  const ML_CalcType I_NCX           = I_NCX_iz + I_NCX_blk;
  const ML_CalcType I_NCX_Na_blk    = 3 * I_NCX_blk;
  const ML_CalcType I_NCX_Na_iz     = 3 * I_NCX_iz;
  const ML_CalcType I_NCX_Ca_blk    = -2 * I_NCX_blk;
  const ML_CalcType I_NCX_Ca_iz     = -2 * I_NCX_iz;

  // currentPMCA()
  const ML_CalcType I_PMCA_blk = v(VT_f_PMCA_blk) * v(VT_Amp_PMCA) *
    pow(Ca_2_blk, 1.6) / (pow(v(VT_K_m), 1.6) + pow(Ca_2_blk, 1.6));
  const ML_CalcType I_PMCA_iz = v(VT_f_PMCA_iz) * v(VT_Amp_PMCA) *
    pow(Ca_2_iz, 1.6) / (pow(v(VT_K_m), 1.6) + pow(Ca_2_iz, 1.6));
  const ML_CalcType I_PMCA = I_PMCA_iz + I_PMCA_blk;

  // CaRU()
  const ML_CalcType p_O_RyR   = Y_ooo + Y_coo + Y_cco + Y_oco;
  const ML_CalcType p_O_RyR_t = p_O_RyR + v(VT_p_O_RyR_base);
  const ML_CalcType J_Ca_rel  = v(VT_P_RyR) * p_O_RyR_t * (Ca_2_SRrl - Ca_2_jnc) * v(VT_Sc_Cell);
  Ca_2_nd_00 = Ca_2_jnc;
  Ca_2_nd_L0 = (Ca_2_nd_00 + pHimP->Ca_2_nd_L02s[Vi]) / pHimP->Ca_2_nd_L0d[Vi];
  Ca_2_nd_0R = (Ca_2_nd_00 + v(VT_f_R) * Ca_2_SRrl) / (1 + v(VT_f_R));
  Ca_2_nd_LR =
    (Ca_2_nd_00 + v(VT_f_R) * Ca_2_SRrl + v(VT_f_L) * 2* v(VT_inverseRTONF) * Vm * expdRTFVm/ (1 - expdRTFVm) * Cao) /
    (1 + v(VT_f_R) + v(VT_f_L) * 2* v(VT_inverseRTONF) * Vm / (1 - expdRTFVm));
  const ML_CalcType epsilon_plus_00 = (Ca_2_nd_00 * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType epsilon_plus_L0 = (Ca_2_nd_L0 * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType epsilon_plus_0R = (Ca_2_nd_0R * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType epsilon_plus_LR = (Ca_2_nd_LR * pHimP->alpha_plus[Vi]) / v(VT_T_L_K_L);
  const ML_CalcType k_co_00         = v(VT_Q_10) * 0.4 / (1 + pow(0.025 / Ca_2_nd_00, 2.7));
  const ML_CalcType k_co_L0         = v(VT_Q_10) * 0.4 / (1 + pow(0.025 / Ca_2_nd_L0, 2.7));
  const ML_CalcType k_co_0R         = v(VT_Q_10) * 0.4 / (1 + pow(0.025 / Ca_2_nd_0R, 2.7));
  const ML_CalcType k_co_LR         = v(VT_Q_10) * 0.4 / (1 + pow(0.025 / Ca_2_nd_LR, 2.7));
  const ML_CalcType f_t_00          = k_co_00 / (k_co_00 + v(VT_k_oc));
  const ML_CalcType f_t_L0          = k_co_L0 / (k_co_L0 + v(VT_k_oc));
  const ML_CalcType k_rco_0         = v(VT_f_n) * f_t_00 * k_co_0R * (v(VT_sloc0) + Ca_2_SRrl);
  const ML_CalcType k_rco_L         = v(VT_f_n) * f_t_L0 * k_co_LR * (v(VT_sloc0) + Ca_2_SRrl);
  const ML_CalcType p_C_0           = v(VT_k_oc) / (v(VT_k_oc) + f_t_00 * (k_rco_0 / (v(VT_f_n) * f_t_00)));
  const ML_CalcType p_C_L           = v(VT_k_oc) / (v(VT_k_oc) + f_t_00 * (k_rco_L / (v(VT_f_n) * f_t_L0)));
  const ML_CalcType k_roc_0         = v(VT_k_oc) * pow(p_C_0, (v(VT_N_RyR) - 1) * 0.74);
  const ML_CalcType k_roc_L         = v(VT_k_oc) * pow(p_C_L, (v(VT_N_RyR) - 1) * 0.74);
  const ML_CalcType Y_ccc           = 1 - (Y_ooo + Y_ooc + Y_coo + Y_coc + Y_cco + Y_oco + Y_occ);
  const ML_CalcType dY_ooodt        = k_rco_L * Y_ooc + pHimP->alpha_plus[Vi] * Y_coo + pHimP->epsilon_minus[Vi] *
    Y_oco - (k_roc_L + pHimP->alpha_minus[Vi] + epsilon_plus_LR) * Y_ooo;
  const ML_CalcType dY_oocdt = pHimP->alpha_plus[Vi] * Y_coc + k_roc_L * Y_ooo + pHimP->epsilon_minus[Vi] *
    Y_occ - (pHimP->alpha_minus[Vi] + k_rco_L + epsilon_plus_L0) * Y_ooc;
  const ML_CalcType dY_coodt = k_rco_0 * Y_coc + pHimP->alpha_minus[Vi] * Y_ooo + pHimP->epsilon_minus[Vi] *
    Y_cco - (k_roc_0 + pHimP->alpha_plus[Vi] + epsilon_plus_0R) * Y_coo;
  const ML_CalcType dY_cocdt = k_roc_0 * Y_coo + pHimP->alpha_minus[Vi] * Y_ooc + pHimP->epsilon_minus[Vi] *
    Y_ccc - (k_rco_0 + pHimP->alpha_plus[Vi] + epsilon_plus_00) * Y_coc;
  const ML_CalcType dY_ccodt = k_rco_0 * Y_ccc + pHimP->alpha_minus[Vi] * Y_oco + epsilon_plus_0R * Y_coo -
    (k_roc_0 + pHimP->alpha_plus[Vi] + pHimP->epsilon_minus[Vi]) * Y_cco;
  const ML_CalcType dY_ocodt = k_rco_0 * Y_occ + pHimP->alpha_plus[Vi] * Y_cco + epsilon_plus_LR * Y_ooo -
    (k_roc_0 + pHimP->alpha_minus[Vi] + pHimP->epsilon_minus[Vi]) * Y_oco;
  const ML_CalcType dY_occdt = pHimP->alpha_plus[Vi] * Y_ccc + k_roc_0 * Y_oco + epsilon_plus_L0 * Y_ooc -
    (pHimP->alpha_minus[Vi] + k_rco_0 + pHimP->epsilon_minus[Vi]) * Y_occ;

  // SERCA()
  const ML_CalcType alpha_2 = 2540 / (1 + pow(v(VT_K_dCai) / Ca_2_blk, 1.7));
  const ML_CalcType alpha_3 = 5.35 / (1 + pow(Ca_2_SRup / v(VT_K_dCasr), 1.7));
  const ML_CalcType beta_1  = 0.1972 / (1 + pow(Ca_2_blk / v(VT_K_dCai), 1.7));
  const ML_CalcType beta_2  = 25435 * v(VT_MgADP_cyt) / (1 + pow(v(VT_K_dCasr) / Ca_2_SRup, 1.7));
  const ML_CalcType beta_3  = 149 * v(VT_Pi);
  const ML_CalcType v_cyc   = 6.86 * (v(VT_alpha_1) * alpha_2 * alpha_3 - beta_1 * beta_2 * beta_3) /
    (alpha_2 * alpha_3 + beta_1 * alpha_3 + beta_1 * beta_2 + v(VT_alpha_1) * alpha_3 + beta_2 * v(VT_alpha_1) +
     beta_2 *
     beta_3 + v(VT_alpha_1) * alpha_2 + beta_3 * beta_1 + beta_3 * alpha_2);
  const ML_CalcType J_SERCA = v(VT_Amp_SERCA) * v_cyc / (2 * v(VT_F)) * v(VT_Sc_Cell);

  // membranePotential()
  const ML_CalcType I_tot_K = (I_CaL_K_jnc + I_CaL_K_iz + I_CaL_K_blk) + I_NaT_K + I_NaL_K + I_K1 + I_Kr +
    (I_Ks_K_iz + I_Ks_K_blk) + I_Kto + I_Kpl + I_NaK_K + I_KATP + I_bNSC_K + (I_l_Ca_K_iz + I_l_Ca_K_blk);
  const ML_CalcType I_tot_Na = (I_CaL_Na_jnc + I_CaL_Na_iz + I_CaL_Na_blk) + (I_NCX_Na_iz + I_NCX_Na_blk) +
    (I_Ks_Na_iz + I_Ks_Na_blk) + I_NaT_Na + I_NaL_Na + I_NaK_Na + I_bNSC_Na + (I_l_Ca_Na_iz + I_l_Ca_Na_blk);
  const ML_CalcType I_tot_Ca_blk = I_CaL_Ca_blk + I_PMCA_blk + I_NCX_Ca_blk + I_Cab_blk;
  const ML_CalcType I_tot_Ca_iz  = I_CaL_Ca_iz + I_PMCA_iz + I_NCX_Ca_iz + I_Cab_iz;
  const ML_CalcType I_tot_Ca_jnc = I_CaL_Ca_LR + I_CaL_Ca_L0;
  const ML_CalcType I_tot_Ca     = I_tot_Ca_jnc + I_tot_Ca_iz + I_tot_Ca_blk;
  const ML_CalcType I_tot_cell   = I_tot_Na + I_tot_Ca + I_tot_K - i_external;

  // contraction()
  const ML_CalcType TS      = v(VT_TS_tot) - TSCa_3 - TSCa_3W - TSCa_3S - TS_S - TS_W;
  const ML_CalcType rate_g  = v(VT_Za) + v(VT_Yv) * (1 - exp(-v(VT_propFh) * pow(hw - v(VT_hwr), 2)));
  const ML_CalcType rate_gd = v(VT_Yd) + v(VT_Yc) *
    pow(v(VT_halfSL) - v(VT_Lc), 2) + v(VT_Yvd) * (1 - exp(-v(VT_propFh) * pow(hw - v(VT_hwr), 2)));
  const ML_CalcType dTSCa_3dt = v(VT_Yb) * TS * pow(Ca_2_blk * 1000, 3) - v(VT_Zb) * TSCa_3 + rate_g * TSCa_3W - v(
    VT_rate_f) * exp(-v(VT_convertF) * pow(v(VT_halfSL) - v(VT_eqvhalfSL), 2)) * TSCa_3;
  const ML_CalcType dTSCa_3Wdt = v(VT_rate_f) *
    exp(-v(VT_convertF) * pow(v(VT_halfSL) - v(VT_eqvhalfSL), 2)) * TSCa_3 - rate_g * TSCa_3W + v(VT_Zp) * TSCa_3S - v(
    VT_Yp) * TSCa_3W;
  const ML_CalcType dTSCa_3Sdt = v(VT_Yp) * TSCa_3W - v(VT_Zp) * TSCa_3S + v(VT_Zr) * TS_S *
    pow(Ca_2_blk * 1000, 3) - v(VT_Yr) * TSCa_3S;
  const ML_CalcType dTS_Sdt = v(VT_Yr) * TSCa_3S - v(VT_Zr) * TS_S *
    pow(Ca_2_blk * 1000, 3) + v(VT_Zq) * TS_W - v(VT_Yq) * TS_S;
  const ML_CalcType dTS_Wdt = v(VT_Yq) * TS_S - v(VT_Zq) * TS_W - rate_gd * TS_W;
  const ML_CalcType dhwdt   = -v(VT_rate_B) * (hw - v(VT_hwr));
  const ML_CalcType dhpdt   = -v(VT_rate_B) * (hp - v(VT_hpr));

  // ionConcentration()
  const ML_CalcType dCa_2_tot_jncdt = -I_tot_Ca_jnc *v(VT_C) / (v(VT_V_jnc) * 2 * v(VT_F)) + J_Ca_rel / v(VT_V_jnc) -
    J_Ca_jnciz / v(VT_V_jnc);
  const ML_CalcType dCa_2_tot_izdt = -I_tot_Ca_iz *v(VT_C) / (v(VT_V_iz) * 2 * v(VT_F)) + J_Ca_jnciz / v(VT_V_iz) -
    J_Ca_izblk / v(VT_V_iz);
  const ML_CalcType dCa_2_tot_blkdt = -I_tot_Ca_blk *v(VT_C) / (v(VT_V_blk) * 2 * v(VT_F)) - J_SERCA / v(VT_V_blk) +
    J_Ca_izblk / v(VT_V_blk);
  const ML_CalcType dCa_2_SRupdt = J_SERCA / v(VT_V_SRup) - J_trans_SR / v(VT_V_SRup);
  const ML_CalcType dCa_2_tot_SRrldt = J_trans_SR / v(VT_V_SRrl) - J_Ca_rel / v(VT_V_SRrl);
  const ML_CalcType dNaidt = -I_tot_Na *v(VT_C) / (v(VT_V_cyt) * v(VT_F));
  const ML_CalcType dKidt = -(-i_external+I_tot_K) * v(VT_C) / (v(VT_V_cyt) * v(VT_F));

  // euler()
  const ML_CalcType tincMs = tinc * 1000;
  CaMCa         += dCaMCadt * tincMs;
  TnChCa        += dTnChCadt * tincMs;
  SRCa          += dSRCadt * tincMs;
  p_O_NaT       += dp_O_NaTdt * tincMs;
  p_I_2_NaT     += dp_I_2_NaTdt * tincMs;
  p_I_s_NaT     += dp_I_s_NaTdt * tincMs;
  p_O_NaL       += dp_O_NaLdt * tincMs;
  p_I_1_NaL     += dp_I_1_NaLdt * tincMs;
  p_I_2_NaL     += dp_I_2_NaLdt * tincMs;
  p_I_s_NaL     += dp_I_s_NaLdt * tincMs;
  chi_r_fast    += dchi_r_fastdt * tincMs;
  chi_r_slow    += dchi_r_slowdt * tincMs;
  para_Xs1      += dpara_Xs1dt * tincMs;
  para_Xs2      += dpara_Xs2dt * tincMs;
  i_fast        += di_fastdt * tincMs;
  i_slow        += di_slowdt * tincMs;
  P_7           += dP_7dt * tincMs;
  P_8_13        += dP_8_13dt * tincMs;
  P_1_6         += dP_1_6dt * tincMs;
  p_E_1_NCX_blk += dp_E_1_NCX_blkdt * tincMs;
  p_I_1_NCX_blk += dp_I_1_NCX_blkdt * tincMs;
  p_I_2_NCX_blk += dp_I_2_NCX_blkdt * tincMs;
  p_E_1_NCX_iz  += dp_E_1_NCX_izdt * tincMs;
  p_I_1_NCX_iz  += dp_I_1_NCX_izdt * tincMs;
  p_I_2_NCX_iz  += dp_I_2_NCX_izdt * tincMs;
  Y_ooo         += dY_ooodt * tincMs;
  Y_ooc         += dY_oocdt * tincMs;
  Y_coo         += dY_coodt * tincMs;
  Y_coc         += dY_cocdt * tincMs;
  Y_cco         += dY_ccodt * tincMs;
  Y_oco         += dY_ocodt * tincMs;
  Y_occ         += dY_occdt * tincMs;
  Y_co_iz       += dY_co_izdt * tincMs;
  Y_oo_iz       += dY_oo_izdt * tincMs;
  Y_oc_iz       += dY_oc_izdt * tincMs;
  Y_co_blk      += dY_co_blkdt * tincMs;
  Y_oo_blk      += dY_oo_blkdt * tincMs;
  Y_oc_blk      += dY_oc_blkdt * tincMs;
  Ca_2_tot_jnc  += dCa_2_tot_jncdt * tincMs;
  Ca_2_tot_iz   += dCa_2_tot_izdt * tincMs;
  Ca_2_tot_blk  += dCa_2_tot_blkdt * tincMs;
  Ca_2_SRup     += dCa_2_SRupdt * tincMs;
  Ca_2_tot_SRrl += dCa_2_tot_SRrldt * tincMs;
  Nai           += dNaidt * tincMs;
  Ki            += dKidt * tincMs;
  TSCa_3        += dTSCa_3dt * tincMs;
  TSCa_3W       += dTSCa_3Wdt * tincMs;
  TSCa_3S       += dTSCa_3Sdt * tincMs;
  TS_S          += dTS_Sdt * tincMs;
  TS_W          += dTS_Wdt * tincMs;
  hw            += dhwdt * tincMs;
  hp            += dhpdt * tincMs;
  Pb_spm        += dPb_spmdt * tincMs;
  a             += dadt * tincMs;

  return tinc*(-I_tot_cell);
}  // Himeno::Calc

void Himeno::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  // Don't forget the blank (' ') at the end!!
  tempstr<<tArg<<' '<<V<<' '
         <<CaMCa<<' '<<TnChCa<<' '<<SRCa<<' '<< p_O_NaT << ' ' << p_I_2_NaT << ' ' << p_I_s_NaT << ' ' << p_O_NaL <<
    ' ' << p_I_1_NaL << ' ' << p_I_2_NaL << ' ' << p_I_s_NaL << ' ' << chi_r_fast << ' ' << chi_r_slow << ' ' <<
    para_Xs1 << ' ' << para_Xs2 << ' ' << i_fast << ' ' << i_slow << ' ' << P_7 << ' ' << P_8_13 << ' ' << P_1_6 <<
    ' ' <<
    p_E_1_NCX_blk << ' ' << p_I_1_NCX_blk << ' ' << p_I_2_NCX_blk << ' ' << p_E_1_NCX_iz << ' ' << p_I_1_NCX_iz <<
    ' ' <<
    p_I_2_NCX_iz << ' ' << Y_ooo << ' ' << Y_ooc << ' ' << Y_coo << ' ' << Y_coc << ' ' << Y_cco << ' ' << Y_oco <<
    ' ' <<
    Y_occ << ' ' << Y_co_iz << ' ' << Y_oo_iz << ' ' << Y_oc_iz << ' ' << Y_co_blk << ' ' << Y_oo_blk << ' ' <<
    Y_oc_blk << ' ' << Ca_2_tot_jnc << ' ' << Ca_2_tot_iz << ' ' << Ca_2_tot_blk << ' ' << Ca_2_tot_SRrl << ' ' <<
    Ca_2_jnc << ' ' << Ca_2_iz << ' ' << Ca_2_blk << ' ' << Ca_2_SRup << ' ' << Ca_2_SRrl << ' ' << Ca_2_nd_00 << ' ' <<
    Ca_2_nd_L0 << ' ' << Ca_2_nd_0R << ' ' << Ca_2_nd_LR << ' ' << Nai << ' ' << Ki << ' ' << TSCa_3 << ' ' <<
    TSCa_3W <<
    ' ' << TSCa_3S << ' ' << TS_S << ' ' << TS_W << ' ' << hw << ' ' << hp << ' ' << Pb_spm << ' ' << a << ' ' <<
    L_bound_iz << ' ' << H_bound_iz << ' ' << L_free_jnc << ' ' << H_free_jnc << ' ';
}

void Himeno::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
  ML_CalcType Vm = V*1000;  // membrane voltage in mV
  const int   Vi = (int)(DivisionTab*(RangeTabhalf+Vm)+.5); // array position

  const ML_CalcType RTONF = v(VT_RTONF);
  const ML_CalcType Cao   = v(VT_Cao);
  const ML_CalcType Ko    = v(VT_Ko);
  const ML_CalcType Nao   = v(VT_Nao);

  // boundaryDiffusion()
  const ML_CalcType J_Ca_jnciz = v(VT_G_dCa_jnciz) * (Ca_2_jnc - Ca_2_iz) * v(VT_Sc_Cell);
  const ML_CalcType J_Ca_izblk = v(VT_G_dCa_izblk) * (Ca_2_iz - Ca_2_blk) * v(VT_Sc_Cell);
  const ML_CalcType J_trans_SR = v(VT_P_trans) * (Ca_2_SRup - Ca_2_SRrl) * v(VT_Sc_Cell);

  // currentCaL()
  const ML_CalcType expdRTFVm    = exp(-2* v(VT_F)/(v(VT_R)*v(VT_Tx))*Vm);
  const ML_CalcType p_O_LCC      = Y_ooo + Y_ooc;
  const ML_CalcType E_K          = RTONF / 1 * log(Ko / Ki);
  const ML_CalcType exp_VdRTF    = exp(-Vm * v(VT_F)/(v(VT_R)*v(VT_Tx)));
  const ML_CalcType GHK_Ca_LR    = 2 * Vm / RTONF * (Ca_2_nd_LR - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_L0    = 2 * Vm / RTONF * (Ca_2_nd_L0 - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_iz    = 2 * Vm / RTONF * (Ca_2_iz - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Ca_blk   = 2 * Vm / RTONF * (Ca_2_blk - Cao * expdRTFVm) / (1 - expdRTFVm);
  const ML_CalcType GHK_Na       = 1 * Vm / RTONF * (Nai - Nao * exp_VdRTF) / (1 - exp_VdRTF);
  const ML_CalcType GHK_K        = 1 * Vm / RTONF * (Ki - Ko * exp_VdRTF) / (1 - exp_VdRTF);
  const ML_CalcType I_CaL_Ca_blk = v(VT_f_CaL_blk) * v(VT_P_CaL_Ca) * GHK_Ca_blk * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_iz  = v(VT_f_CaL_iz) * v(VT_P_CaL_Ca) * GHK_Ca_iz * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_LR  = v(VT_f_CaL_jnc) * v(VT_P_CaL_Ca) * GHK_Ca_LR * Y_ooo * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Ca_L0  = v(VT_f_CaL_jnc) * v(VT_P_CaL_Ca) * GHK_Ca_L0 * Y_ooc * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_blk = v(VT_f_CaL_blk) * v(VT_P_CaL_Na) * GHK_Na * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_iz  = v(VT_f_CaL_iz) * v(VT_P_CaL_Na) * GHK_Na * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_Na_jnc = v(VT_f_CaL_jnc) * v(VT_P_CaL_Na) * GHK_Na * p_O_LCC * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_blk  = v(VT_f_CaL_blk) * v(VT_P_CaL_K) * GHK_K * Y_oo_blk * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_iz   = v(VT_f_CaL_iz) * v(VT_P_CaL_K) * GHK_K * Y_oo_iz * v(VT_ATPfactor);
  const ML_CalcType I_CaL_K_jnc  = v(VT_f_CaL_jnc) * v(VT_P_CaL_K) * GHK_K * p_O_LCC * v(VT_ATPfactor);
  const ML_CalcType I_CaL        = (I_CaL_Ca_LR + I_CaL_Ca_L0 + I_CaL_Na_jnc + I_CaL_K_jnc) +
    (I_CaL_Ca_iz + I_CaL_Na_iz + I_CaL_K_iz) + (I_CaL_Ca_blk + I_CaL_Na_blk + I_CaL_K_blk);

  // currentNa()
  const ML_CalcType I_NaT_Na = (1 - v(VT_f_LSM)) * v(VT_P_Na) * GHK_Na * p_O_NaT;
  const ML_CalcType I_NaT_K  = (1 - v(VT_f_LSM)) * v(VT_P_Na) * 0.1 * GHK_K * p_O_NaT;
  const ML_CalcType I_NaT    = I_NaT_Na + I_NaT_K;
  const ML_CalcType I_NaL_Na = v(VT_f_LSM) * v(VT_P_Na) * GHK_Na * p_O_NaL;
  const ML_CalcType I_NaL_K  = v(VT_f_LSM) * v(VT_P_Na) * 0.1 * GHK_K * p_O_NaL;
  const ML_CalcType I_NaL    = I_NaL_Na + I_NaL_K;
  const ML_CalcType I_Na     = I_NaT + I_NaL;

  // currentK1()
  const ML_CalcType alpha_Mg = 12.0 * exp(-0.025 * (Vm - E_K));
  const ML_CalcType beta_Mg  = 28 * v(VT_Mg_2_cyt) * exp(0.025 * (Vm - E_K));
  const ML_CalcType f_O      = alpha_Mg / (alpha_Mg + beta_Mg);
  const ML_CalcType f_B      = beta_Mg / (alpha_Mg + beta_Mg);
  const ML_CalcType po_Mg    = f_O * f_O * f_O;
  const ML_CalcType po_Mg1   = 3.0 * f_O * f_O * f_B;
  const ML_CalcType po_Mg2   = 3.0 * f_O * f_B * f_B;
  const ML_CalcType po_mode1 = v(VT_f_mode1) * (1 - Pb_spm) * (po_Mg + (2.0 / 3.0) * po_Mg1 + (1.0 / 3.0) * po_Mg2);
  const ML_CalcType po_mode2 = (1 - v(VT_f_mode1)) / (1.0 + v(VT_SPM) / (40.0 * exp(-(Vm - E_K) / 9.1)));
  const ML_CalcType p_O_K1   = po_mode1 + po_mode2;
  const ML_CalcType I_K1     = v(VT_G_K1) * v(VT_chi_K1) * (Vm - E_K) * p_O_K1;

  // currentKr()
  const ML_CalcType chi_r  = pHimP->A_chi_r_fast[Vi] * chi_r_fast + pHimP->A_chi_r_slow[Vi] * chi_r_slow;
  const ML_CalcType p_O_Kr = chi_r * pHimP->R_Kr[Vi];
  const ML_CalcType I_Kr   = v(VT_G_Kr) * v(VT_chi_Kr) * (Vm - E_K) * p_O_Kr;

  // currentsKs()
  const ML_CalcType para_RKs_blk = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_blk, 1.4));
  const ML_CalcType para_RKs_iz  = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_iz, 1.4));
  const ML_CalcType p_O_Ks_blk   = para_Xs1 * para_Xs2 * para_RKs_blk;
  const ML_CalcType p_O_Ks_iz    = para_Xs1 * para_Xs2 * para_RKs_iz;
  const ML_CalcType I_Ks_K_blk   = v(VT_f_Ks_blk) * v(VT_P_Ks_K) * GHK_K * p_O_Ks_blk;
  const ML_CalcType I_Ks_K_iz    = v(VT_f_Ks_iz) * v(VT_P_Ks_K) * GHK_K * p_O_Ks_iz;
  const ML_CalcType I_Ks_Na_blk  = v(VT_f_Ks_blk) * v(VT_P_Ks_Na) * GHK_Na * p_O_Ks_blk;
  const ML_CalcType I_Ks_Na_iz   = v(VT_f_Ks_iz) * v(VT_P_Ks_Na) * GHK_Na * p_O_Ks_iz;
  const ML_CalcType I_Ks         = I_Ks_K_blk + I_Ks_K_iz + I_Ks_Na_blk + I_Ks_Na_iz;

  // currentKto()
  const ML_CalcType i       = pHimP->A_i_fast[Vi] * i_fast + pHimP->A_i_slow[Vi] * i_slow;
  const ML_CalcType p_O_Kto = a * i;
  const ML_CalcType I_Kto   = v(VT_G_Kto) * p_O_Kto * (Vm - E_K);

  // currentKpl()
  const ML_CalcType I_Kpl = v(VT_P_Kpl) * v(VT_chi_Kpl) * pHimP->p_O_Kpl[Vi] * GHK_K;

  // currentCab()
  const ML_CalcType I_Cab_blk = v(VT_P_Cab) * v(VT_f_Cab_blk) * GHK_Ca_blk;
  const ML_CalcType I_Cab_iz  = v(VT_P_Cab) * v(VT_f_Cab_iz) * GHK_Ca_iz;
  const ML_CalcType I_Cab     = I_Cab_iz + I_Cab_blk;

  // currentbNSC()
  const ML_CalcType I_bNSC_K  = v(VT_P_bNSC_K) * GHK_K;
  const ML_CalcType I_bNSC_Na = v(VT_P_bNSC_Na) * GHK_Na;
  const ML_CalcType I_bNSC    = I_bNSC_K + I_bNSC_Na;

  // currentlCa()
  const ML_CalcType p_O_blk       = 1.0 / (1.0 + pow(0.0012 / Ca_2_blk, 3));
  const ML_CalcType p_O_iz        = 1.0 / (1.0 + pow(0.0012 / Ca_2_iz, 3));
  const ML_CalcType I_l_Ca_Na_blk = v(VT_P_l_Ca_Na) * v(VT_f_l_Ca_blk) * GHK_Na * p_O_blk;
  const ML_CalcType I_l_Ca_Na_iz  = v(VT_P_l_Ca_Na) * v(VT_f_l_Ca_iz) * GHK_Na * p_O_iz;
  const ML_CalcType I_l_Ca_K_blk  = v(VT_P_l_Ca_K) * v(VT_f_l_Ca_blk) * GHK_K * p_O_blk;
  const ML_CalcType I_l_Ca_K_iz   = v(VT_P_l_Ca_K) * v(VT_f_l_Ca_iz) * GHK_K * p_O_iz;
  const ML_CalcType I_l_Ca        = I_l_Ca_Na_iz + I_l_Ca_K_iz + I_l_Ca_Na_blk + I_l_Ca_K_blk;

  // currentKATP()
  const ML_CalcType I_KATP = v(VT_G_KATP) * (Vm - E_K) * v(VT_p_O_KATP) * v(VT_chi_KATP);

  // currentNaK()
  const ML_CalcType V_step2   = v(VT_alpha_2_plus) * P_7 - pHimP->alpha_2_minus[Vi] * P_8_13;
  const ML_CalcType v_cyc_NaK = V_step2;
  const ML_CalcType I_NaK     = v(VT_Amp_NaK) * v_cyc_NaK;
  const ML_CalcType I_NaK_Na  = 3 * I_NaK;
  const ML_CalcType I_NaK_K   = -2 * I_NaK;

  // currentNCX()
  const ML_CalcType q_blk_E_1_Na  = 1.0 / (1.0 + pow(v(VT_K_m_Nai) / Nai, 3) * (1.0 + Ca_2_blk / v(VT_K_m_Cai)));
  const ML_CalcType q_iz_E_1_Na   = 1.0 / (1.0 + pow(v(VT_K_m_Nai) / Nai, 3) * (1.0 + Ca_2_iz / v(VT_K_m_Cai)));
  const ML_CalcType p_E_2_NCX_blk = 1 - p_E_1_NCX_blk - p_I_1_NCX_blk - p_I_2_NCX_blk;
  const ML_CalcType p_E_2_NCX_iz  = 1 - p_E_1_NCX_iz - p_I_1_NCX_iz - p_I_2_NCX_iz;
  const ML_CalcType v_cyc_NCX_blk = pHimP->k_1[Vi] * q_blk_E_1_Na * p_E_1_NCX_blk - pHimP->k_2[Vi] * v(VT_q_E_2_Na) *
    p_E_2_NCX_blk;
  const ML_CalcType v_cyc_NCX_iz = pHimP->k_1[Vi] * q_iz_E_1_Na * p_E_1_NCX_iz - pHimP->k_2[Vi] * v(VT_q_E_2_Na) *
    p_E_2_NCX_iz;
  const ML_CalcType I_NCX_blk    = v(VT_f_NCX_blk) * v(VT_Amp_NCX) * v_cyc_NCX_blk;
  const ML_CalcType I_NCX_iz     = v(VT_f_NCX_iz) * v(VT_Amp_NCX) * v_cyc_NCX_iz;
  const ML_CalcType I_NCX        = I_NCX_iz + I_NCX_blk;
  const ML_CalcType I_NCX_Na_blk = 3 * I_NCX_blk;
  const ML_CalcType I_NCX_Na_iz  = 3 * I_NCX_iz;
  const ML_CalcType I_NCX_Ca_blk = -2 * I_NCX_blk;
  const ML_CalcType I_NCX_Ca_iz  = -2 * I_NCX_iz;

  // currentPMCA()
  const ML_CalcType I_PMCA_blk = v(VT_f_PMCA_blk) * v(VT_Amp_PMCA) *
    pow(Ca_2_blk, 1.6) / (pow(v(VT_K_m), 1.6) + pow(Ca_2_blk, 1.6));
  const ML_CalcType I_PMCA_iz = v(VT_f_PMCA_iz) * v(VT_Amp_PMCA) *
    pow(Ca_2_iz, 1.6) / (pow(v(VT_K_m), 1.6) + pow(Ca_2_iz, 1.6));
  const ML_CalcType I_PMCA = I_PMCA_iz + I_PMCA_blk;

  // SERCA()
  const ML_CalcType alpha_2 = 2540 / (1 + pow(v(VT_K_dCai) / Ca_2_blk, 1.7));
  const ML_CalcType alpha_3 = 5.35 / (1 + pow(Ca_2_SRup / v(VT_K_dCasr), 1.7));
  const ML_CalcType beta_1  = 0.1972 / (1 + pow(Ca_2_blk / v(VT_K_dCai), 1.7));
  const ML_CalcType beta_2  = 25435 * v(VT_MgADP_cyt) / (1 + pow(v(VT_K_dCasr) / Ca_2_SRup, 1.7));
  const ML_CalcType beta_3  = 149 * v(VT_Pi);
  const ML_CalcType v_cyc   = 6.86 * (v(VT_alpha_1) * alpha_2 * alpha_3 - beta_1 * beta_2 * beta_3) /
    (alpha_2 * alpha_3 + beta_1 * alpha_3 + beta_1 * beta_2 + v(VT_alpha_1) * alpha_3 + beta_2 * v(VT_alpha_1) +
     beta_2 *
     beta_3 + v(VT_alpha_1) * alpha_2 + beta_3 * beta_1 + beta_3 * alpha_2);
  const ML_CalcType J_SERCA = v(VT_Amp_SERCA) * v_cyc / (2 * v(VT_F)) * v(VT_Sc_Cell);

  // membranePotential()
  const ML_CalcType I_tot_K = (I_CaL_K_jnc + I_CaL_K_iz + I_CaL_K_blk) + I_NaT_K + I_NaL_K + I_K1 + I_Kr +
    (I_Ks_K_iz + I_Ks_K_blk) + I_Kto + I_Kpl + I_NaK_K + I_KATP + I_bNSC_K + (I_l_Ca_K_iz + I_l_Ca_K_blk);
  const ML_CalcType I_tot_Na = (I_CaL_Na_jnc + I_CaL_Na_iz + I_CaL_Na_blk) + (I_NCX_Na_iz + I_NCX_Na_blk) +
    (I_Ks_Na_iz + I_Ks_Na_blk) + I_NaT_Na + I_NaL_Na + I_NaK_Na + I_bNSC_Na + (I_l_Ca_Na_iz + I_l_Ca_Na_blk);
  const ML_CalcType I_tot_Ca_blk = I_CaL_Ca_blk + I_PMCA_blk + I_NCX_Ca_blk + I_Cab_blk;
  const ML_CalcType I_tot_Ca_iz  = I_CaL_Ca_iz + I_PMCA_iz + I_NCX_Ca_iz + I_Cab_iz;
  const ML_CalcType I_tot_Ca_jnc = I_CaL_Ca_LR + I_CaL_Ca_L0;
  const ML_CalcType I_tot_Ca     = I_tot_Ca_jnc + I_tot_Ca_iz + I_tot_Ca_blk;
  const ML_CalcType I_tot_cell   = I_tot_Na + I_tot_Ca + I_tot_K;


  tempstr << I_tot_cell << " " << I_tot_Ca << " " << I_tot_Na << " " << I_tot_K << " " << I_tot_Ca_blk << " " <<
    I_tot_Ca_iz << " " << I_tot_Ca_jnc <<" "<< J_SERCA << " " << I_PMCA << " " <<I_PMCA_iz << " "<< I_PMCA_blk<< " "<<
    I_NCX_Ca_iz<<" "<<I_NCX_Ca_blk << " " <<I_NCX_Na_iz<< " " <<I_NCX_Na_blk << " " << I_NCX << " " <<I_NCX_iz<< " " <<
    I_NCX_blk << " " << I_NaK_K<< " " << I_NaK_Na << " "<< I_NaK << " " << I_KATP << " " << I_l_Ca << " " <<
    I_l_Ca_K_iz << " " << I_l_Ca_K_blk << " " << I_l_Ca_Na_iz << " " << I_l_Ca_Na_blk << " " << I_bNSC<< " " <<
    I_bNSC_Na << " "<< I_bNSC_K << " " << I_Cab << " " << I_Cab_iz << " " << I_Cab_blk << " " << I_Kpl << " " <<
    I_Kto <<
    " " << I_Ks << " " << I_Ks_Na_iz << " " << I_Ks_Na_blk << " " << I_Ks_K_iz<< " " <<I_Ks_K_blk << " " << I_Kr <<
    " " <<
    I_K1 << " " << I_Na << " " << I_NaL << " " << I_NaL_K << " " << I_NaL_Na << " " << I_NaT << " " << I_NaT_K << " " <<
    I_NaT_Na<< " " << I_CaL << " " << I_CaL_K_jnc << " " << I_CaL_K_iz << " " << I_CaL_K_blk << " " << I_CaL_Na_jnc <<
    " " << I_CaL_Na_iz << " " << I_CaL_Na_blk << " " << I_CaL_Ca_L0 << " " << I_CaL_Ca_LR << " " <<I_CaL_Ca_iz << " " <<
    I_CaL_Ca_blk <<' ';
}  // Himeno::LongPrint

void Himeno::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {"CaMCa",         "TnChCa",         "SRCa",                   "p_O_NaT",                   "p_I_2_NaT",
   "p_I_s_NaT",
   "p_O_NaL",
   "p_I_1_NaL",     "p_I_2_NaL",      "p_I_s_NaL",
   "chi_r_fast",    "chi_r_slow",     "para_Xs1",               "para_Xs2",                  "i_fast",
   "i_slow",
   "P_7",
   "P_8_13",        "P_1_6",          "p_E_1_NCX_blk",
   "p_I_1_NCX_blk", "p_I_2_NCX_blk",  "p_E_1_NCX_iz",           "p_I_1_NCX_iz",              "p_I_2_NCX_iz",
   "Y_ooo",
   "Y_ooc",
   "Y_coo",
   "Y_coc",         "Y_cco",          "Y_oco",                  "Y_occ",                     "Y_co_iz",
   "Y_oo_iz",
   "Y_oc_iz",
   "Y_co_blk",      "Y_oo_blk",       "Y_oc_blk",
   "Ca_2_tot_jnc",  "Ca_2_tot_iz",    "Ca_2_tot_blk",           "Ca_2_tot_SRrl",             "Ca_2_jnc",
   "Ca_2_iz",
   "Ca_2_blk",
   "Ca_2_SRup",
   "Ca_2_SRrl",     "Ca_2_nd_00",     "Ca_2_nd_L0",             "Ca_2_nd_0R",                "Ca_2_nd_LR",
   "Nai",
   "Ki",
   "TSCa_3",        "TSCa_3W",        "TSCa_3S",
   "TS_S",          "TS_W",           "hw",                     "hp",                        "Pb_spm",
   "a",
   "L_bound_iz",
   "H_bound_iz",    "L_free_jnc",     "H_free_jnc"};

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
} // Himeno::GetParameterNames

void Himeno::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {"I_tot_cell",    "I_tot_Ca",        "I_tot_Na",           "I_tot_K",             "I_tot_Ca_blk",
   "I_tot_Ca_iz",
   "I_tot_Ca_jnc",
   "J_SERCA",       "I_PMCA",
   "I_PMCA_iz",     "I_PMCA_blk",      "I_NCX_Ca_iz",        "I_NCX_Ca_blk",        "I_NCX_Na_iz",
   "I_NCX_Na_blk",
   "I_NCX",
   "I_NCX_iz",
   "I_NCX_blk",     "I_NaK_K",         "I_NaK_Na",           "I_NaK",               "I_KATP",
   "I_l_Ca",
   "I_l_Ca_K_iz",
   "I_l_Ca_K_blk",  "I_l_Ca_Na_iz",
   "I_l_Ca_Na_blk", "I_bNSC",          "I_bNSC_Na",          "I_bNSC_K",            "I_Cab",
   "I_Cab_iz",
   "I_Cab_blk",
   "I_Kpl",         "I_Kto",           "I_Ks",
   "I_Ks_Na_iz",    "I_Ks_Na_blk",     "I_Ks_K_iz",          "I_Ks_K_blk",          "I_Kr",                     "I_K1",
   "I_Na",
   "I_NaL",         "I_NaL_K",         "I_NaL_Na",
   "I_NaT",         "I_NaT_K",         "I_NaT_Na",           "I_CaL",               "I_CaL_K_jnc",
   "I_CaL_K_iz",
   "I_CaL_K_blk",
   "I_CaL_Na_jnc",  "I_CaL_Na_iz",
   "I_CaL_Na_blk",  "I_CaL_Ca_L0",     "I_CaL_Ca_LR",        "I_CaL_Ca_iz",         "I_CaL_Ca_blk"};
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
} // Himeno::GetLongParameterNames
