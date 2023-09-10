/**@file BeelerReuterParameters.cpp
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

#include <BeelerReuterParameters.h>

BeelerReuterParameters::BeelerReuterParameters(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "BeelerReuterParameters::BeelerReuterParameters "<< initFile << endl;
#endif  // if KADEBUG

  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

void BeelerReuterParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"BeelerReuterParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}

void BeelerReuterParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "BeelerReuterParameters:init " << initFile << endl;
#endif  // if KADEBUG

  // Initialization of parameter names
  P[VT_C_m].name       = "C_m";
  P[VT_g_Na].name      = "g_Na";
  P[VT_g_NaC].name     = "g_NaC";
  P[VT_E_Na].name      = "E_Na";
  P[VT_g_s].name       = "g_s";
  P[VT_Vol].name       = "Vol";
  P[VT_Amp].name       = "Amp";
  P[VT_dC_m].name      = "dC_m";
  P[VT_Init_Ca_i].name = "init_Ca_i";
  P[VT_Init_m].name    = "init_m";
  P[VT_Init_h].name    = "init_h";
  P[VT_Init_j].name    = "init_j";
  P[VT_Init_d].name    = "init_d";
  P[VT_Init_f].name    = "init_f";
  P[VT_Init_x1].name   = "init_x1";
  P[VT_Init_Vm].name   = "init_Vm";
  P[VT_RC_ax1C1].name  = "RC_ax1C1";
  P[VT_RC_ax1C2].name  = "RC_ax1C2";
  P[VT_RC_ax1C3].name  = "RC_ax1C3";
  P[VT_RC_ax1C4].name  = "RC_ax1C4";
  P[VT_RC_ax1C5].name  = "RC_ax1C5";
  P[VT_RC_ax1C6].name  = "RC_ax1C6";
  P[VT_RC_ax1C7].name  = "RC_ax1C7";
  P[VT_RC_bx1C1].name  = "RC_bx1C1";
  P[VT_RC_bx1C2].name  = "RC_bx1C2";
  P[VT_RC_bx1C3].name  = "RC_bx1C3";
  P[VT_RC_bx1C4].name  = "RC_bx1C4";
  P[VT_RC_bx1C5].name  = "RC_bx1C5";
  P[VT_RC_bx1C6].name  = "RC_bx1C6";
  P[VT_RC_bx1C7].name  = "RC_bx1C7";
  P[VT_RC_amC1].name   = "RC_amC1";
  P[VT_RC_amC2].name   = "RC_amC2";
  P[VT_RC_amC3].name   = "RC_amC3";
  P[VT_RC_amC4].name   = "RC_amC4";
  P[VT_RC_amC5].name   = "RC_amC5";
  P[VT_RC_amC6].name   = "RC_amC6";
  P[VT_RC_amC7].name   = "RC_amC7";
  P[VT_RC_bmC1].name   = "RC_bmC1";
  P[VT_RC_bmC2].name   = "RC_bmC2";
  P[VT_RC_bmC3].name   = "RC_bmC3";
  P[VT_RC_bmC4].name   = "RC_bmC4";
  P[VT_RC_bmC5].name   = "RC_bmC5";
  P[VT_RC_bmC6].name   = "RC_bmC6";
  P[VT_RC_bmC7].name   = "RC_bmC7";
  P[VT_RC_ahC1].name   = "RC_ahC1";
  P[VT_RC_ahC2].name   = "RC_ahC2";
  P[VT_RC_ahC3].name   = "RC_ahC3";
  P[VT_RC_ahC4].name   = "RC_ahC4";
  P[VT_RC_ahC5].name   = "RC_ahC5";
  P[VT_RC_ahC6].name   = "RC_ahC6";
  P[VT_RC_ahC7].name   = "RC_ahC7";
  P[VT_RC_bhC1].name   = "RC_bhC1";
  P[VT_RC_bhC2].name   = "RC_bhC2";
  P[VT_RC_bhC3].name   = "RC_bhC3";
  P[VT_RC_bhC4].name   = "RC_bhC4";
  P[VT_RC_bhC5].name   = "RC_bhC5";
  P[VT_RC_bhC6].name   = "RC_bhC6";
  P[VT_RC_bhC7].name   = "RC_bhC7";
  P[VT_RC_ajC1].name   = "RC_ajC1";
  P[VT_RC_ajC2].name   = "RC_ajC2";
  P[VT_RC_ajC3].name   = "RC_ajC3";
  P[VT_RC_ajC4].name   = "RC_ajC4";
  P[VT_RC_ajC5].name   = "RC_ajC5";
  P[VT_RC_ajC6].name   = "RC_ajC6";
  P[VT_RC_ajC7].name   = "RC_ajC7";
  P[VT_RC_bjC1].name   = "RC_bjC1";
  P[VT_RC_bjC2].name   = "RC_bjC2";
  P[VT_RC_bjC3].name   = "RC_bjC3";
  P[VT_RC_bjC4].name   = "RC_bjC4";
  P[VT_RC_bjC5].name   = "RC_bjC5";
  P[VT_RC_bjC6].name   = "RC_bjC6";
  P[VT_RC_bjC7].name   = "RC_bjC7";
  P[VT_RC_adC1].name   = "RC_adC1";
  P[VT_RC_adC2].name   = "RC_adC2";
  P[VT_RC_adC3].name   = "RC_adC3";
  P[VT_RC_adC4].name   = "RC_adC4";
  P[VT_RC_adC5].name   = "RC_adC5";
  P[VT_RC_adC6].name   = "RC_adC6";
  P[VT_RC_adC7].name   = "RC_adC7";
  P[VT_RC_bdC1].name   = "RC_bdC1";
  P[VT_RC_bdC2].name   = "RC_bdC2";
  P[VT_RC_bdC3].name   = "RC_bdC3";
  P[VT_RC_bdC4].name   = "RC_bdC4";
  P[VT_RC_bdC5].name   = "RC_bdC5";
  P[VT_RC_bdC6].name   = "RC_bdC6";
  P[VT_RC_bdC7].name   = "RC_bdC7";
  P[VT_RC_afC1].name   = "RC_afC1";
  P[VT_RC_afC2].name   = "RC_afC2";
  P[VT_RC_afC3].name   = "RC_afC3";
  P[VT_RC_afC4].name   = "RC_afC4";
  P[VT_RC_afC5].name   = "RC_afC5";
  P[VT_RC_afC6].name   = "RC_afC6";
  P[VT_RC_afC7].name   = "RC_afC7";
  P[VT_RC_bfC1].name   = "RC_bfC1";
  P[VT_RC_bfC2].name   = "RC_bfC2";
  P[VT_RC_bfC3].name   = "RC_bfC3";
  P[VT_RC_bfC4].name   = "RC_bfC4";
  P[VT_RC_bfC5].name   = "RC_bfC5";
  P[VT_RC_bfC6].name   = "RC_bfC6";
  P[VT_RC_bfC7].name   = "RC_bfC7";
  P[VT_RC_ix1C1].name  = "RC_ix1C1";
  P[VT_RC_ix1C2].name  = "RC_ix1C2";
  P[VT_RC_ix1C3].name  = "RC_ix1C3";
  P[VT_RC_ix1C4].name  = "RC_ix1C4";
  P[VT_RC_ix1C5].name  = "RC_ix1C5";
  P[VT_RC_ix1C6].name  = "RC_ix1C6";
  P[VT_RC_ix1C7].name  = "RC_ix1C7";

  P[VT_dC_m].readFromFile = false;

  ParameterLoader EPL(initFile, EMT_BeelerReuter);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  Calculate();
  InitTable(tinc);
}  // BeelerReuterParameters::Init

void BeelerReuterParameters::Calculate() {
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
#if KADEBUG
  cerr << "BeelerReuterParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_dC_m].value = .001/(P[VT_C_m].value);
}

void BeelerReuterParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "BeelerReuterParameters - InitTable()" << endl;
#endif  // if KADEBUG

  // double *p;
  // double oldValue;

  double a, b;

  tinc *= -1000.0;  // sec -> ms., negation for exptau calculation

  for (double V = -RangeTabhalf+0.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);

    a =
      (P[VT_RC_ax1C1].value*exp(P[VT_RC_ax1C2].value*(V+(P[VT_RC_ax1C3].value)))+(P[VT_RC_ax1C4].value)*
       (V+(P[VT_RC_ax1C5].value)))/(exp(P[VT_RC_ax1C6].value*(V+(P[VT_RC_ax1C3].value)))+(P[VT_RC_ax1C7].value));
    b = a +
      (P[VT_RC_bx1C1].value*exp(P[VT_RC_bx1C2].value*(V+(P[VT_RC_bx1C3].value)))+(P[VT_RC_bx1C4].value)*
       (V+(P[VT_RC_bx1C5].value)))/(exp(P[VT_RC_bx1C6].value*(V+(P[VT_RC_bx1C3].value)))+(P[VT_RC_bx1C7].value));
    x1_inf[Vi]    = a / b;
    exptau_x1[Vi] = exp(tinc * b);

    a =
      (P[VT_RC_amC1].value*exp(P[VT_RC_amC2].value*(V+(P[VT_RC_amC3].value)))+(P[VT_RC_amC4].value)*
       (V+(P[VT_RC_amC5].value)))/(exp(P[VT_RC_amC6].value*(V+(P[VT_RC_amC3].value)))+(P[VT_RC_amC7].value));
    b = a +
      (P[VT_RC_bmC1].value*exp(P[VT_RC_bmC2].value*(V+(P[VT_RC_bmC3].value)))+(P[VT_RC_bmC4].value)*
       (V+(P[VT_RC_bmC5].value)))/(exp(P[VT_RC_bmC6].value*(V+(P[VT_RC_bmC3].value)))+(P[VT_RC_bmC7].value));
    m_inf[Vi]    = a / b;
    exptau_m[Vi] = exp(tinc * b);

    a =
      (P[VT_RC_ahC1].value*exp(P[VT_RC_ahC2].value*(V+(P[VT_RC_ahC3].value)))+(P[VT_RC_ahC4].value)*
       (V+(P[VT_RC_ahC5].value)))/(exp(P[VT_RC_ahC6].value*(V+(P[VT_RC_ahC3].value)))+(P[VT_RC_ahC7].value));
    b = a +
      (P[VT_RC_bhC1].value*exp(P[VT_RC_bhC2].value*(V+(P[VT_RC_bhC3].value)))+(P[VT_RC_bhC4].value)*
       (V+(P[VT_RC_bhC5].value)))/(exp(P[VT_RC_bhC6].value*(V+(P[VT_RC_bhC3].value)))+(P[VT_RC_bhC7].value));
    h_inf[Vi]    = a / b;
    exptau_h[Vi] = exp(tinc * b);

    a =
      (P[VT_RC_ajC1].value*exp(P[VT_RC_ajC2].value*(V+(P[VT_RC_ajC3].value)))+(P[VT_RC_ajC4].value)*
       (V+(P[VT_RC_ajC5].value)))/(exp(P[VT_RC_ajC6].value*(V+(P[VT_RC_ajC3].value)))+(P[VT_RC_ajC7].value));
    b = a +
      (P[VT_RC_bjC1].value*exp(P[VT_RC_bjC2].value*(V+(P[VT_RC_bjC3].value)))+(P[VT_RC_bjC4].value)*
       (V+(P[VT_RC_bjC5].value)))/(exp(P[VT_RC_bjC6].value*(V+(P[VT_RC_bjC3].value)))+(P[VT_RC_bjC7].value));
    j_inf[Vi]    = a / b;
    exptau_j[Vi] = exp(tinc * b);

    a =
      (P[VT_RC_adC1].value*exp(P[VT_RC_adC2].value*(V+(P[VT_RC_adC3].value)))+(P[VT_RC_adC4].value)*
       (V+(P[VT_RC_adC5].value)))/(exp(P[VT_RC_adC6].value*(V+(P[VT_RC_adC3].value)))+(P[VT_RC_adC7].value));
    b = a +
      (P[VT_RC_bdC1].value*exp(P[VT_RC_bdC2].value*(V+(P[VT_RC_bdC3].value)))+(P[VT_RC_bdC4].value)*
       (V+(P[VT_RC_bdC5].value)))/(exp(P[VT_RC_bdC6].value*(V+(P[VT_RC_bdC3].value)))+(P[VT_RC_bdC7].value));
    d_inf[Vi]    = a / b;
    exptau_d[Vi] = exp(tinc * b);

    a =
      (P[VT_RC_afC1].value*exp(P[VT_RC_afC2].value*(V+(P[VT_RC_afC3].value)))+(P[VT_RC_afC4].value)*
       (V+(P[VT_RC_afC5].value)))/(exp(P[VT_RC_afC6].value*(V+(P[VT_RC_afC3].value)))+(P[VT_RC_afC7].value));
    b = a +
      (P[VT_RC_bfC1].value*exp(P[VT_RC_bfC2].value*(V+(P[VT_RC_bfC3].value)))+(P[VT_RC_bfC4].value)*
       (V+(P[VT_RC_bfC5].value)))/(exp(P[VT_RC_bfC6].value*(V+(P[VT_RC_bfC3].value)))+(P[VT_RC_bfC7].value));
    f_inf[Vi]    = a / b;
    exptau_f[Vi] = exp(tinc * b);

    /*for (int y=0; y<12; y++)
       {
        switch(y)
          {
          case 0: {oldValue = 0; p=&a_x1[Vi];break;}
          case 1: {oldValue = *p; p=&b_x1[Vi];break;}
          case 2: {oldValue = 0; p=&a_m[Vi];break;}
          case 3: {oldValue = *p; p=&b_m[Vi];break;}
          case 4: {oldValue = 0; p=&a_h[Vi];break;}
          case 5: {oldValue = *p; p=&b_h[Vi];break;}
          case 6: {oldValue = 0; p=&a_j[Vi];break;}
          case 7: {oldValue = *p; p=&b_j[Vi];break;}
          case 8: {oldValue = 0; p=&a_d[Vi];break;}
          case 9: {oldValue = *p; p=&b_d[Vi];break;}
          case 10: {oldValue = 0; p=&a_f[Vi];break;}
          case 11: {oldValue = *p; p=&b_f[Vi];break;}
          }
     * p=oldValue+(v(VT_RC[y]C1)*exp(v(VT_RC[y]C2)*(V+(v(VT_RC[y]C3))))+(v(VT_RC[y]C4))*(V+(v(VT_RC[y]C5))))/(exp(v(VT_RC[y]C6)*(V+(v(VT_RC[y]C3))))+(v(VT_RC[y]C7)));
       }*/

    i_x_1[Vi] = P[VT_RC_ix1C1].value*(exp(P[VT_RC_ix1C2].value*(V+(P[VT_RC_ix1C3].value)))-(P[VT_RC_ix1C4].value))/exp(
      P[VT_RC_ix1C5].value*(V+(P[VT_RC_ix1C6].value)));
    i_K_1[Vi] = 0.35*
      (4.0*(exp(.04*(V+85.0))-1.0)/(exp(.08*(V+53.0))+exp(.04*(V+53.0)))+.2*(V+23.0)/(1.0-exp(-.04*(V+23.0))));

    /*cout << Vi << ' ' << a_x1[Vi] << ' ' << b_x1[Vi] << ' ' << a_m[Vi] << ' ' << b_m[Vi] << ' ' << a_h[Vi] << ' '
       <<b_h[Vi]
       << ' ' << a_j[Vi] << ' ' << b_j[Vi] << ' ' << a_d[Vi] << ' ' << b_d[Vi] << ' ' << a_f[Vi] << ' ' << b_f[Vi]
       << ' ' << i_K_1_[Vi] << ' ' << i_x_1_[Vi] << endl;*/

    // p=NULL;
    // delete p;
  }
}  // BeelerReuterParameters::InitTable
