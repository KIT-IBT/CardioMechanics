/*
 * File: ForceModelTest.cpp
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





#include <CellModel.h>
#include <dwArray.h>
#include <kaVersion.h>
#include <time.h>


PrintParameterModus PrintParameterMode = PrintParameterModeOff;

int main(int argc, char **argv) {
  kaVersion EMTVersion(1, 1, 0, 2023, 9, 1, argc, argv);

  EMTVersion.addOption("tbegin", "<time [s]>", OT_optional, "0", "start time of the simulation");
  EMTVersion.addOption("tend", "<time [s]>", OT_optional, "1.5", "end time of the simulation");
  EMTVersion.addOption("tinc", "<time [s]>", OT_optional, "1e-5", "time increment for the calculation");
  EMTVersion.addOption("toutbegin", "<time [s]>", OT_optional, "0", "start time of the output");
  EMTVersion.addOption("toutend", "<time [s]>", OT_optional, "99999", "end time of the output");
  EMTVersion.addOption("toutinc", "<value>", OT_optional, "10", "factor of the output increment related to tinc");
  EMTVersion.addOption("outfile", "<file>", OT_optional, "", "file in which the output is printed");
  EMTVersion.addOption("outputtime", "<file>", OT_optional, "",
                       "format: lines containing only one double value setting the times where an output is produced");
  EMTVersion.addOption("outputtimeasset", "", OT_optional, "false", "");
  EMTVersion.addOption("save", "<file>", OT_optional, "", "save the current cell state in a backup file");
  EMTVersion.addOption("load", "<file>", OT_optional, "", "start the calculation with a backup file");
  EMTVersion.addOption("wait", "<file>", OT_optional, "", "");
  EMTVersion.addOption("verbose", "", OT_optional, "false", "output of additional informations");
  EMTVersion.addOption("quiet", "", OT_optional, "false", "");
  EMTVersion.addOption("progress", "", OT_optional, "false", "show progress on stderr");
#ifdef USE_EMT
  EMTVersion.addOption("amplitude", "<value>", OT_optional, "cellmodel specific",
                       "amplitude of the stimulation current");
  EMTVersion.addOption("thresholdvm", "<value>", OT_optional, "cellmodel specific", "");
  EMTVersion.addOption("text", "<value>", OT_optional, "1", "time interval between the stimulations");
  EMTVersion.addOption("textlen", "<value>", OT_optional, "cellmodel specific", "length of the stimulation");
  EMTVersion.addOption("textbegin", "<value>", OT_optional, "0.1", "start time of the first stimulation");
  EMTVersion.addOption("textend", "<value>", OT_optional, "cellmodel specific", "end time of the first stimulation");
  EMTVersion.addOption("textdiv", "<value>", OT_optional, "1", "");
  EMTVersion.addOption("clamp", "<file>", OT_optional, "", "use a clampfile for clamping the transmembran voltage");
#endif  // ifdef USE_EMT
  EMTVersion.addOption("stretch", "<value [m/m]>", OT_optional, "1", "");
  EMTVersion.addOption("sl", "<value mue m>", OT_optional, "2", "");
  EMTVersion.addOption("velocity", "<value [stretch/s]>", OT_optional, "0", "");
  EMTVersion.addOption("stretchbegin", "<time [s]>", OT_optional, "0", "");
  EMTVersion.addOption("stretchend", "<time [s]>", OT_optional, "1", "");
  EMTVersion.addOption("restretchtime", "<time [s]>", OT_optional, "cellmodel specific", "");
  EMTVersion.addOption("bp|-bpi|-sp|-lp|-lpi", "", OT_optional, "StatusPrintNormal", "");
  EMTVersion.addOption("pp", "", OT_optional, "off", "");
  EMTVersion.addOption("euler|-rk2|-rk4", "", OT_optional, "euler", "");
  EMTVersion.addOption("h", "", OT_optional, "", "print help");
  EMTVersion.addOption("multitest", "<amount>", OT_optional, "0", "");
#ifdef USE_EMT
  EMTVersion.addOption("CaO", "<value>", OT_optional, "cellmodel specific", "");
  EMTVersion.addOption("opt", "<value>", OT_optional, "false", "");
  EMTVersion.addOption("freq", "<value>", OT_optional, "cellmodel specific", "");
  EMTVersion.addOption("monophasic", "", OT_optional, "true", "");
  EMTVersion.addOption("biphasic", "", OT_optional, "false", "");
  EMTVersion.addOption("R", "<cell#> <cell#> <resistance [Ohm]>", OT_optional, "", "");
  EMTVersion.addOption("evfile", "<file>", OT_optional, "",
                       "evfile with the name of the applied cellmodel and its parameters");
#endif  // ifdef USE_EMT
#ifdef USE_FMT
  EMTVersion.addOption("fvfile", "<file>", OT_optional, "", "");
# ifndef USE_EMT
  EMTVersion.addOption("tCabegin", "<value>", OT_optional, "", "");
  EMTVersion.addOption("CaRice", "", OT_optional, "", "");
  EMTVersion.addOption("CaInc", "", OT_optional, "", "");
  EMTVersion.addOption("CaPanerai", "", OT_optional, "", "");
  EMTVersion.addOption("CaCoppini", "", OT_optional, "", "");
  EMTVersion.addOption("CaConst", "<value>", OT_optional, "", "");
# else  // ifndef USE_EMT
  EMTVersion.addOption("feedback", "", OT_optional, "false", "");
# endif  // ifndef USE_EMT
#endif  // ifdef USE_FMT


  if (argc < 2) {
    EMTVersion.printHelp();
#ifdef USE_EMT
    for (int emt = EMT_Dummy+1; emt != EMT_Last; emt++)
      cerr << "\t[-" << FindElphyModelDescriptor((ElphyModelType)emt) << "]"  << endl;
#endif  // ifdef USE_EMT
#ifdef USE_FMT
    for (int fmt = FMT_Dummy+1; fmt != FMT_Last; fmt++)
      cerr << "\t[-" << FindForceModelDescriptor((ForceModelType)fmt) << "]"  << endl;
#endif  // ifdef USE_FMT
    exit(-1);
  }

  int tout = 0, toutinc = 10;
  double t, tbegin = 0.0, tend = 1.5, tinc = .00001, toutbegin = 0., toutend = 99999;
  bool   verbose = false, quiet = false, progress = false;

  /*     added for outputtime support     */
  char *outputtime;
  bool  setoutputtime = false;
  dwArray<double> times;
  int index            = 0;
  bool outputtimeasset = false;

  /*     added for outputtime support     */
  StatusPrint spm = StatusPrintNormal;

  vector<CellParameters> cp;
#ifdef USE_EMT
  vector<CellCoupling> cc;
#endif  // ifdef USE_EMT

  CellParameters dummyCP;
  dummyCP.cellnr = 0;
  cp.push_back(dummyCP);

  CellParameters *actCP = &cp.back();

  try {
    for (int ac = 1; ac < argc; ac++) {
      const char *parg = argv[ac];
      if (!strcasecmp("-h", parg) || !strcasecmp("-help", parg)) {
        EMTVersion.printHelp();
      } else if (!strcasecmp("-verbose", parg) || !strcasecmp("-v", parg)) {
        verbose = true;
      } else if (!strcasecmp("-quiet", parg)) {
        quiet = true;
      } else if (!strcasecmp("-progress", parg)) {
        progress = true;
      } else if (!strcasecmp("-lp", parg)) {
        spm = StatusPrintLong;
      } else if (!strcasecmp("-lpi", parg)) {
        spm = StatusPrintLongI;
      } else if (!strcasecmp("-sp", parg)) {
        spm = StatusPrintShort;
      } else if (!strcasecmp("-bp", parg)) {
        spm = StatusPrintBrief;
      } else if (!strcasecmp("-bpi", parg)) {
        spm = StatusPrintBriefI;
      } else if (!strcasecmp("-tbegin", parg) && (ac+1 < argc) ) {
        tbegin = atof(argv[++ac]);
      } else if (!strcasecmp("-tend", parg) && (ac+1 < argc) ) {
        tend = atof(argv[++ac]);
      } else if (!strcasecmp("-tinc", parg) && (ac+1 < argc) ) {
        tinc = atof(argv[++ac]);
      } else if (!strcasecmp("-toutbegin", parg) && (ac+1 < argc) ) {
        toutbegin = atof(argv[++ac]);
      } else if (!strcasecmp("-toutend", parg) && (ac+1 < argc) ) {
        toutend = atof(argv[++ac]);
      } else if (!strcasecmp("-toutinc", parg) && (ac+1 < argc) ) {
        toutinc = (int)atof(argv[++ac]);
      } else if (!strcasecmp("-multitest", parg) && (ac+1 < argc) ) {
        int mt = atoi(argv[++ac]);
        for (int num = 0; num < mt; num++)
          cp.push_back(*actCP);
      } else if ((*parg == '-') && isdigit(parg[1])) {
        int cellnr = atoi(parg+1);
        vector<CellParameters>::iterator iCP, bCP = cp.begin(), eCP = cp.end();
        for (iCP = bCP; iCP != eCP; iCP++)
          if ((*iCP).cellnr == cellnr)
            break;
        if (iCP != eCP) {
          actCP = &(*iCP);
        } else {
          cp.push_back(dummyCP);
          actCP         = &cp.back();
          actCP->cellnr = cellnr;
        }
      }

#if defined(USE_FMT) && defined(USE_EMT)
      else if (!strcasecmp("-feedback", parg)) {
        actCP->fb = true;
      }
#endif  // if defined(USE_FMT) && defined(USE_EMT)
#ifdef USE_EMT
      else if (!strcasecmp("-text", parg) && (ac+1 < argc)) {
        actCP->text = atof(argv[++ac]);
      } else if (!strcasecmp("-textlen", parg) && (ac+1 < argc) ) {
        actCP->textlen = atof(argv[++ac]);
      } else if (!strcasecmp("-textbegin", parg) && (ac+1 < argc) ) {
        actCP->textbegin = atof(argv[++ac]);
      } else if (!strcasecmp("-textend", parg) && (ac+1 < argc) ) {
        actCP->textend = atof(argv[++ac]);
      } else if (!strcasecmp("-textdiv", parg) && (ac+1 < argc) ) {
        actCP->textdiv = atof(argv[++ac]);
      } else if (!strcasecmp("-evfile", parg) && (ac+1 < argc) ) {
        actCP->initEVFile = argv[++ac];
      } else if ((*parg == '-') && (FindElphyModelType(parg+1) != EMT_Last) ) {
        actCP->emt = FindElphyModelType(parg+1);
      } else if (!strcasecmp("-CaO", parg) && (ac+1 < argc) ) {
        actCP->Ca_o = atof(argv[++ac]);
      } else if (!strcasecmp("-opt", parg)) {
        actCP->opt = true;
      } else if (!strcasecmp("-freq", parg) && (ac+1 < argc) ) {
        actCP->freq = atof(argv[++ac]);
      } else if (!strcasecmp("-monophasic", parg)) {
        actCP->monophasic = true;
      } else if (!strcasecmp("-biphasic", parg)) {
        actCP->monophasic = false;
      } else if (!strcasecmp("-amplitude", parg) && (ac+1 < argc) ) {
        actCP->amplitude = atof(argv[++ac]);
      } else if (!strcasecmp("-thresholdvm", parg) && (ac+1 < argc) ) {
        actCP->thresVm = atof(argv[++ac]);
      } else if (!strcasecmp("-clamp", parg) && (ac+1 < argc) ) {
        actCP->readClampFile(argv[++ac]);
      } else if (!strcasecmp("-R", parg) && (ac+3 < argc) ) {
        CellCoupling dummyCC;
        dummyCC.cellnr1 = atoi(argv[++ac]);
        dummyCC.cellnr2 = atoi(argv[++ac]);
        dummyCC.R       = atof(argv[++ac]);
        cc.push_back(dummyCC);
      }
#endif  // ifdef USE_EMT
#ifdef USE_FMT
      else if (!strcasecmp("-fvfile", parg) && (ac+1 < argc)) {
        actCP->initFVFile = argv[++ac];
      } else if ((*parg == '-') && (FindForceModelType(parg+1) != FMT_Last) ) {
        actCP->fmt = FindForceModelType(parg+1);
      }
# ifndef USE_EMT
      else if (!strcasecmp("-CaRice", parg)) {
        actCP->cf = CF_Rice;
      } else if (!strcasecmp("-CaInc", parg)) {
        actCP->cf = CF_Inc;
      } else if (!strcasecmp("-CaPanerai", parg)) {
        actCP->cf = CF_Panerai;
      } else if (!strcasecmp("-CaCoppini", parg)) {
        actCP->cf = CF_Coppini;
      } else if (!strcasecmp("-CaConst", parg) && (ac+1 < argc) ) {
        actCP->cf      = CF_Const;
        actCP->CaConst = atof(argv[++ac]);
      } else if (!strcasecmp("-tCabegin", parg) && (ac+1 < argc) ) {
        actCP->tCabegin = atof(argv[++ac]);
      }
# endif  // ifndef USE_EMT
#endif  // ifdef USE_FMT
      else if (!strcasecmp("-euler", parg)) {
        actCP->cs = CS_Euler;
      } else if (!strcasecmp("-rk2", parg)) {
        actCP->cs = CS_RK2;
      } else if (!strcasecmp("-rk4", parg)) {
        actCP->cs = CS_RK4;
      } else if (!strcasecmp("-stretchAmplitude", parg) && (ac+1 < argc) ) {
        actCP->stretchAmplitude = atof(argv[++ac]);
        actCP->sf               = SF_RECT;
      } else if (!strcasecmp("-stretch", parg) && (ac+1 < argc)) {
        actCP->stretchdefault = atof(argv[++ac]);
      } else if (!strcasecmp("-sl", parg) && (ac+1 < argc) ) {
        double sl = atof(argv[++ac]);
        actCP->stretchdefault = sl/2.0;
      } else if (!strcasecmp("-stretchbegin", parg) && (ac+1 < argc)) {
        actCP->stretchbegin = atof(argv[++ac]);
      } else if (!strcasecmp("-stretchend", parg) && (ac+1 < argc) ) {
        actCP->stretchend = atof(argv[++ac]);
      } else if (!strcasecmp("-restretchtime", parg) && (ac+1 < argc) ) {
        actCP->restretch = atof(argv[++ac]);
      } else if (!strcasecmp("-velocity", parg) && (ac+1 < argc) ) {
        actCP->velocitydefault = atof(argv[++ac]);
      } else if (!strcasecmp("-pp", parg)) {
        PrintParameterMode = PrintParameterModeOn;
      } else if (!strcasecmp("-outfile", parg) && (ac+1 < argc) ) {
        actCP->outfile = argv[++ac];
      }

      /*     added for outputtime support     */
      else if (!strcasecmp("-outputtime", parg) && (ac+1 < argc)) {
        outputtime    = argv[++ac];
        setoutputtime = true;
      } else if (!strcasecmp("-outputtimeasset", parg)) {
        outputtimeasset = true;
      }

      /*     added for outputtime support     */
      else if (!strcasecmp("-load", parg) && (ac+1 < argc)) {
        actCP->loadStatus = argv[++ac];
      } else if (!strcasecmp("-save", parg) && (ac+1 < argc) ) {
        actCP->saveStatus = argv[++ac];
      } else {
        throw kaBaseException("Unknown or incomplete option %s", parg);
      }
    }

    /*     added for outputtime support     */
    if (setoutputtime == true) {
      double zwischen = 0;
      ifstream in(outputtime);

      if (!in.is_open()) {
        throw kaBaseException("File %s could not be opened!", outputtime);
      } else {
        char instr[1024];
        do {
          in.getline(instr, sizeof(instr));
          int check = sscanf(instr, "%lf", &zwischen);
          if (check == 1)
            times.addEntry(zwischen);
          else if (check == 0)
            throw kaBaseException("Please check the format of your input-file %s!\n", outputtime);
        } while (!in.eof());
      }
    }

    /*     added for outputtime support     */

    vector<CellParameters>::iterator iCP, bCP = cp.begin(), eCP = cp.end();
    for (iCP = bCP; iCP != eCP; iCP++) {
      (*iCP).init(verbose, tinc);
      if ((*iCP).loadStatus)
        (*iCP).loadStatusFile();
      if (!quiet)
        (*iCP).outputStatusLine(spm);
    }

#ifdef USE_EMT
    vector<CellCoupling>::iterator iCC, bCC = cc.begin(), eCC = cc.end();
    for (iCC = bCC; iCC != eCC; iCC++) {
      const double sf = 1e9;  // Vm [mV], Cm [mue F] => R [Ohm] with scaling by 1e-9
      (*iCC).S = ((*iCC).R ? sf/(*iCC).R : 0.);
      for (iCP = bCP; iCP != eCP; iCP++)
        if ((*iCP).cellnr == (*iCC).cellnr1)
          break;
      if (iCP == eCP)
        throw kaBaseException("Cell %d cannot be found for coupling to cell %d", (*iCC).cellnr1, (*iCC).cellnr2);
      (*iCC).cell1 = &(*iCP);
      for (iCP = bCP; iCP != eCP; iCP++)
        if ((*iCP).cellnr == (*iCC).cellnr2)
          break;
      if (iCP == eCP)
        throw kaBaseException("Cell %d cannot be found for coupling to cell %d", (*iCC).cellnr2, (*iCC).cellnr1);
      (*iCC).cell2 = &(*iCP);
    }

    for (iCP = bCP; iCP != eCP; iCP++)
      if (PrintParameterMode == PrintParameterModeOn)
        (*iCP).ppem->PrintParameters();

#endif  // ifdef USE_EMT


    int pcnt = 0;
    for (t = tbegin; t < tend; t += tinc) {
      if (progress) {
        if (pcnt++ == 0)
          fprintf(stderr, "\rProgress: % 6.2f%%", (t-tbegin)/(tend-tbegin)*100);
        else if (pcnt > 10) pcnt = 0;
      }

#ifdef USE_EMT
      for (iCP = bCP; iCP != eCP; iCP++)
        (*iCP).i_internal = 0.;

      for (iCC = bCC; iCC != eCC; iCC++) {
        CalcType i_internal = ((*iCC).cell2->Vm-(*iCC).cell1->Vm)*(*iCC).S;
        (*iCC).cell1->i_internal += i_internal;
        (*iCC).cell2->i_internal += -i_internal;
      }
#endif  // ifdef USE_EMT

      for (iCP = bCP; iCP != eCP; iCP++)
        (*iCP).calculate(t, tinc);

      /*     added for outputtime support     */
      if (setoutputtime == true) {
        if ((index <= times.uBound()) && (t >= times[index])) {
          for (iCP = bCP; iCP != eCP; iCP++)
            if (!quiet) {
              if (outputtimeasset == true)
                (*iCP).output(spm, times[index]);
              else
                (*iCP).output(spm, t);
            }
          index++;
        }
      }

      /*     added for outputtime support     */
      else {
        if ((t >= toutbegin) && (t < toutend) ) {
          if (!tout) {
            for (iCP = bCP; iCP != eCP; iCP++)
              if (!quiet)
                (*iCP).output(spm, t);
          }
          if (++tout >= toutinc)
            tout = 0;
        }
      }
    }

    if (progress) {
      fprintf(stderr, "\n");
    }


    for (iCP = bCP; iCP != eCP; iCP++)
      if ((*iCP).saveStatus)
        (*iCP).saveStatusFile();
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }

  return 0;
}  // main
