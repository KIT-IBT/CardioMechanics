/*
 * File: CellModel.h
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


#include <kaBasicIO.h>
using namespace nskaGlobal;

// needed for kaWrite/kaRead

#include <CellModelLayer.h>
#include <vector>

#include <limits>
#undef HUGE
#define HUGE std::numeric_limits<ML_CalcType>::max()

typedef ML_CalcType CalcType;

enum StatusPrint {
  StatusPrintBrief, StatusPrintBriefI, StatusPrintShort, StatusPrintLong, StatusPrintLongI, StatusPrintNormal
};
enum StretchFunction { SF_VELOCITY, SF_RECT, SF_NONE };
enum ForceCoupling { TroponinCoupling, CalciumCoupling, CalciumFeedBackCoupling };
enum CaFunction { CF_Default, CF_Rice, CF_Panerai, CF_Const, CF_Inc, CF_Coppini };
enum CellSolver { CS_Euler, CS_RK2, CS_RK4 };

//! Class handles single step in clamping protocol
class Clamp {
 public:
  double t;
  int type;
  double val;
};

//! Class handles parameters for a single cell model
class CellParameters {
 public:
  int cellnr;

#ifdef USE_EMT
  CalcType Vm;
  ElphyModelType emt;
  string initEVFile;
  vbElphyModel<CalcType> *pem;
  vbElphyParameters<CalcType> *ppem;

  CalcType Ca_o;
  CalcType i_internal;
  double text, textlen, textbegin, textend, textdiv, freq;
  double amplitude, thresVm;
  bool opt, monophasic;
  double currentEndTime;
  vector<Clamp> vClamp;
  int currentClampIndex;
  typedef CalcType (vbElphyModel<CalcType>::*csMethod)(double, CalcType, CalcType, CalcType, int);

  csMethod csFunction;

#endif  // ifdef USE_EMT

#ifdef USE_FMT
  CalcType Force;
  ForceModelType fmt;
  string initFVFile;
  vbForceModel<CalcType> *pfm;
  vbForceParameters<CalcType> *ppfm;

# ifdef USE_EMT
  ForceCoupling fc;
  bool fb;
# else  // ifdef USE_EMT
  CaFunction cf;
  double tCabegin;
  CalcType CaConst;
  vector<double> extCaTransient;
# endif  // ifdef USE_EMT
#endif  // ifdef USE_FMT
  double stretch;
  double stretchdefault, stretchbegin, stretchend, restretch, stretchAmplitude;
  double velocity, velocitydefault;
  StretchFunction sf;

  CalcType Cai;
  char *loadStatus;
  char *saveStatus;

  char *outfile;
  FILE *fpo;

  CellSolver cs;

  CellParameters() {
    cellnr = 0;
#ifdef USE_EMT
    Vm                = 0.0;
    emt               = EMT_Dummy;
    initEVFile        = "";
    pem               = NULL;
    ppem              = NULL;
    Ca_o              = HUGE;
    i_internal        = 0.;
    text              = 1.0;
    textbegin         = 0.1;
    textlen           = HUGE;
    textend           = HUGE;
    textdiv           = 1;
    freq              = HUGE;
    amplitude         = HUGE;
    thresVm           = HUGE;
    opt               = false;
    monophasic        = true;
    currentClampIndex = 0;
    currentEndTime    = 0;
#endif  // ifdef USE_EMT

#ifdef USE_FMT
    Force      = 0.;
    fmt        = FMT_Dummy;
    initFVFile = "";
    pfm        = NULL;
    ppfm       = NULL;
# ifdef USE_EMT
    fb = false;
# else  // ifdef USE_EMT
    cf       = CF_Default;
    tCabegin = 0.;
    CaConst  = 0.;
# endif  // ifdef USE_EMT
#endif  // ifdef USE_FMT

    stretch          = stretchdefault = 1.;
    stretchbegin     = 0.;
    stretchend       = 0.;
    velocity         = velocitydefault = 0.;
    restretch        = HUGE;
    stretchAmplitude = 1.;

    Cai        = 0.;
    sf         = SF_NONE;
    loadStatus = NULL;
    saveStatus = NULL;

    outfile = NULL;
    fpo     = NULL;

    cs = CS_Euler;
  }

  virtual ~CellParameters() {
    if (fpo)
      fclose(fpo);

#ifdef USE_EMT
    delete ppem;
    delete pem;
#endif  // ifdef USE_EMT
#ifdef USE_FMT
    delete pfm;
    delete ppfm;
#endif  // ifdef USE_FMT
  }

#ifdef USE_EMT

  // Read clamp configuration from file clampFile into vector

  void readClampFile(char *clampFile) {
    FILE *fp = fopen(clampFile, "r");

    if (!fp)
      throw kaBaseException("Can't open clampfile %s", clampFile);
    for (int numClamps = 0; !feof(fp); numClamps++) {
      Clamp tmp;
      char  buf[256];
      int   rc = fscanf(fp, "%s %lf %lf ", buf, &tmp.t, &tmp.val);
      if (rc < 3)
        break;
      if (!strcasecmp(buf, "Vm")) {
        tmp.type = 1;
      } else if (!strcasecmp(buf, "Cai")) {
        tmp.type = 2;
      } else if (!strcasecmp(buf, "Cai+")) {
        tmp.type = 3; 
        currentClampIndex = 1;
      } else if (!strcasecmp(buf, "stretch")) {
        tmp.type = 4;
      } else {
        throw kaBaseException("Unknown clamp type %s", buf);
      }

      // cerr << "# "<<numClamps<<" "<<tmp.type<<" "<<tmp.t<<" "<<tmp.val<<endl;
      vClamp.push_back(tmp);
    }
    if (vClamp.size() > 0) {
      currentEndTime = vClamp[0].t;

      // cerr<<"t=0, setting to "<<vClamp[currentClampIndex].val<<" until "<<currentEndTime<<endl;
    }
    fclose(fp);
  }  // readClampFile

#endif  // ifdef USE_EMT

  void init(bool verbose, ML_CalcType tinc) {
#ifdef USE_EMT
    if (initEVFile.size())
      emt = FindElphyModelFileType(initEVFile.c_str());
    else
      initEVFile = FindElphyModelFileName(emt);

    if (verbose)
      cerr << "ElphyModelType " << emt << " used in " << initEVFile << endl;

    initElphyParameters<CalcType>(&ppem, initEVFile.c_str(), emt, tinc);
    initElphyModel<CalcType>(&pem, ppem, emt, opt);

    if (amplitude == HUGE)
      amplitude = pem->GetAmplitude();
    if (textlen == HUGE) {
      textlen = pem->GetStimTime();
      if (verbose)
        cout << "Setting StimTime to " << textlen << endl;
    }
    Vm = pem->GetVm();

    /*!< Function to set the new tinc value given by parameter at command line */
    pem->SetTinc(tinc);

    if (Ca_o != HUGE)
      pem->SetCao(Ca_o);
#endif  // USE_EMT

#ifdef USE_FMT
    if (initFVFile.size())
      fmt = FindForceModelFileType(initFVFile.c_str());

    if (!initFVFile.size())
      initFVFile = FindForceModelFileName(fmt);
    if (verbose)
      cerr << "ForceModelType " << fmt << " used in " << initFVFile << endl;

    initForceParameters(&ppfm, initFVFile.c_str(), fmt);
    initForceModel(&pfm, ppfm, fmt);

# ifdef USE_EMT // and USE_FMT
    // Inquiry if CellModel gives TCa and if ForceModel can use it
    if (pem->OutIsTCa() && pfm->InIsTCa())
      fc = TroponinCoupling;
    else
      fc = (fb ? CalciumFeedBackCoupling : CalciumCoupling);
# else // USE_EMT
    if (cf == CF_Coppini) {       // human CaT used in Land et al. 2017 model (ref. 9 there, array from matlab code from Niederer lab web page; R. Coppini, C. Ferrantini, L. Yao, P. Fan, M. Del Lungo, F. Stillitano, L. Sartiani, B. Tosi, S. Suffredini, C. Tesi, M. Yacoub, I. Olivotto, L. Belardinelli, C. Poggesi, E. Cerbai, A. Mugelli, Late sodium current inhibition reverses electromechan- ical dysfunction in human hypertrophic cardiomyopathy, Circulation 127 (5) (February 2013) 575â584.
      double tempCaTransient[] =
      {0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.170, 0.175, 0.181, 0.188, 0.196, 0.205, 0.215, 0.225,
       0.235, 0.247, 0.259, 0.270, 0.280, 0.289, 0.297, 0.304, 0.312, 0.319, 0.328, 0.336, 0.345, 0.353, 0.359, 0.364,
       0.368, 0.373, 0.377, 0.382, 0.388, 0.394, 0.401, 0.408, 0.413, 0.416, 0.420, 0.423, 0.427, 0.432, 0.438, 0.444,
       0.450, 0.455, 0.460, 0.463, 0.465, 0.468, 0.470, 0.475, 0.480, 0.485, 0.491, 0.496, 0.499, 0.502, 0.504, 0.504,
       0.506, 0.510, 0.513, 0.519, 0.524, 0.529, 0.531, 0.532, 0.532, 0.532, 0.533, 0.535, 0.538, 0.543, 0.547, 0.551,
       0.554, 0.554, 0.554, 0.554, 0.553, 0.555, 0.559, 0.561, 0.565, 0.568, 0.570, 0.569, 0.568, 0.568, 0.566, 0.567,
       0.570, 0.574, 0.577, 0.581, 0.581, 0.581, 0.579, 0.578, 0.578, 0.578, 0.579, 0.582, 0.586, 0.588, 0.588, 0.587,
       0.585, 0.583, 0.581, 0.580, 0.583, 0.587, 0.590, 0.593, 0.593, 0.592, 0.590, 0.586, 0.584, 0.585, 0.586, 0.590,
       0.592, 0.594, 0.593, 0.593, 0.590, 0.588, 0.586, 0.586, 0.587, 0.589, 0.592, 0.594, 0.593, 0.591, 0.589, 0.586,
       0.584, 0.583, 0.585, 0.587, 0.589, 0.591, 0.589, 0.588, 0.585, 0.582, 0.580, 0.579, 0.580, 0.582, 0.585, 0.587,
       0.586, 0.584, 0.581, 0.579, 0.578, 0.577, 0.578, 0.580, 0.581, 0.583, 0.581, 0.579, 0.576, 0.572, 0.570, 0.569,
       0.572, 0.574, 0.576, 0.578, 0.577, 0.574, 0.572, 0.569, 0.567, 0.566, 0.567, 0.569, 0.570, 0.570, 0.569, 0.567,
       0.563, 0.560, 0.558, 0.557, 0.558, 0.559, 0.561, 0.562, 0.563, 0.560, 0.557, 0.554, 0.551, 0.550, 0.550, 0.552,
       0.553, 0.555, 0.553, 0.550, 0.546, 0.544, 0.541, 0.541, 0.542, 0.543, 0.545, 0.545, 0.543, 0.541, 0.537, 0.534,
       0.531, 0.531, 0.530, 0.532, 0.534, 0.535, 0.533, 0.530, 0.527, 0.523, 0.520, 0.518, 0.518, 0.519, 0.521, 0.521,
       0.520, 0.518, 0.515, 0.512, 0.509, 0.508, 0.508, 0.508, 0.510, 0.510, 0.508, 0.505, 0.502, 0.499, 0.497, 0.496,
       0.496, 0.496, 0.496, 0.496, 0.496, 0.493, 0.490, 0.486, 0.484, 0.483, 0.483, 0.483, 0.484, 0.483, 0.481, 0.479,
       0.476, 0.472, 0.469, 0.467, 0.467, 0.467, 0.468, 0.469, 0.466, 0.464, 0.460, 0.457, 0.454, 0.452, 0.453, 0.453,
       0.454, 0.454, 0.453, 0.450, 0.447, 0.444, 0.441, 0.439, 0.439, 0.440, 0.441, 0.441, 0.440, 0.437, 0.433, 0.430,
       0.427, 0.426, 0.426, 0.426, 0.426, 0.426, 0.424, 0.421, 0.418, 0.415, 0.413, 0.411, 0.410, 0.411, 0.412, 0.411,
       0.409, 0.407, 0.404, 0.401, 0.399, 0.398, 0.397, 0.397, 0.397, 0.397, 0.396, 0.393, 0.390, 0.387, 0.384, 0.383,
       0.383, 0.383, 0.384, 0.384, 0.382, 0.380, 0.377, 0.374, 0.372, 0.370, 0.369, 0.370, 0.370, 0.370, 0.368, 0.366,
       0.364, 0.361, 0.358, 0.357, 0.357, 0.357, 0.357, 0.356, 0.355, 0.353, 0.350, 0.348, 0.345, 0.344, 0.343, 0.344,
       0.344, 0.344, 0.343, 0.341, 0.338, 0.336, 0.333, 0.332, 0.332, 0.332, 0.332, 0.332, 0.332, 0.330, 0.327, 0.325,
       0.323, 0.322, 0.321, 0.321, 0.321, 0.321, 0.320, 0.318, 0.316, 0.313, 0.311, 0.310, 0.310, 0.310, 0.311, 0.311,
       0.310, 0.308, 0.306, 0.303, 0.301, 0.301, 0.300, 0.300, 0.301, 0.301, 0.300, 0.298, 0.296, 0.294, 0.292, 0.291,
       0.291, 0.291, 0.292, 0.292, 0.291, 0.289, 0.287, 0.285, 0.284, 0.283, 0.282, 0.283, 0.283, 0.282, 0.281, 0.280,
       0.279, 0.277, 0.275, 0.274, 0.274, 0.274, 0.275, 0.274, 0.273, 0.272, 0.270, 0.269, 0.268, 0.267, 0.267, 0.267,
       0.267, 0.267, 0.266, 0.265, 0.263, 0.261, 0.260, 0.260, 0.260, 0.260, 0.260, 0.260, 0.259, 0.258, 0.257, 0.255,
       0.253, 0.253, 0.253, 0.253, 0.254, 0.254, 0.253, 0.252, 0.250, 0.248, 0.247, 0.246, 0.246, 0.247, 0.247, 0.247,
       0.247, 0.246, 0.244, 0.243, 0.241, 0.241, 0.241, 0.242, 0.242, 0.242, 0.242, 0.240, 0.239, 0.237, 0.236, 0.236,
       0.236, 0.236, 0.237, 0.237, 0.236, 0.236, 0.234, 0.233, 0.232, 0.231, 0.232, 0.232, 0.233, 0.232, 0.232, 0.231,
       0.229, 0.228, 0.227, 0.226, 0.226, 0.227, 0.228, 0.228, 0.227, 0.226, 0.225, 0.224, 0.223, 0.222, 0.222, 0.223,
       0.223, 0.223, 0.223, 0.222, 0.221, 0.220, 0.219, 0.218, 0.219, 0.219, 0.220, 0.220, 0.220, 0.219, 0.218, 0.217,
       0.216, 0.215, 0.215, 0.216, 0.216, 0.216, 0.216, 0.215, 0.214, 0.213, 0.212, 0.211, 0.212, 0.212, 0.213, 0.213,
       0.213, 0.212, 0.211, 0.209, 0.209, 0.209, 0.209, 0.209, 0.210, 0.210, 0.209, 0.208, 0.207, 0.206, 0.205, 0.205,
       0.205, 0.206, 0.206, 0.206, 0.206, 0.206, 0.205, 0.204, 0.203, 0.203, 0.203, 0.204, 0.204, 0.204, 0.204, 0.203,
       0.202, 0.201, 0.200, 0.200, 0.200, 0.201, 0.201, 0.201, 0.201, 0.200, 0.199, 0.199, 0.198, 0.197, 0.198, 0.198,
       0.199, 0.199, 0.199, 0.198, 0.197, 0.196, 0.195, 0.195, 0.196, 0.196, 0.197, 0.197, 0.197, 0.196, 0.195, 0.194,
       0.194, 0.193, 0.194, 0.194, 0.195, 0.195, 0.195, 0.194, 0.193, 0.192, 0.191, 0.191, 0.191, 0.192, 0.192, 0.193,
       0.192, 0.192, 0.191, 0.190, 0.190, 0.189, 0.190, 0.190, 0.191, 0.191, 0.191, 0.190, 0.189, 0.189, 0.188, 0.188,
       0.188, 0.188, 0.189, 0.189, 0.189, 0.188, 0.187, 0.187, 0.186, 0.186, 0.186, 0.187, 0.187, 0.188, 0.187, 0.187,
       0.186, 0.185, 0.185, 0.184, 0.185, 0.185, 0.186, 0.186, 0.186, 0.185, 0.184, 0.184, 0.183, 0.183, 0.183, 0.184,
       0.185, 0.185, 0.185, 0.184, 0.183, 0.182, 0.182, 0.181, 0.182, 0.182, 0.183, 0.184, 0.184, 0.183, 0.182, 0.181,
       0.181, 0.181, 0.181, 0.181, 0.182, 0.182, 0.182, 0.182, 0.181, 0.180, 0.179, 0.179, 0.180, 0.180, 0.181, 0.181,
       0.181, 0.180, 0.180, 0.179, 0.178, 0.178, 0.178, 0.179, 0.180, 0.180, 0.180, 0.179, 0.179, 0.178, 0.177, 0.177,
       0.177, 0.178, 0.179, 0.179, 0.179, 0.178, 0.178, 0.177, 0.176, 0.176, 0.176, 0.177, 0.178, 0.178, 0.178, 0.177,
       0.176, 0.175, 0.175, 0.175, 0.175, 0.176, 0.176, 0.177, 0.177, 0.177, 0.176, 0.175, 0.175, 0.175, 0.175, 0.176,
       0.176, 0.177, 0.176, 0.176, 0.175, 0.174, 0.174, 0.174, 0.174, 0.174, 0.175, 0.175, 0.175, 0.175, 0.174, 0.173,
       0.173, 0.173, 0.173, 0.174, 0.174, 0.174, 0.174, 0.174, 0.173, 0.172, 0.172, 0.172, 0.172, 0.172, 0.173, 0.173,
       0.173, 0.173, 0.172, 0.171, 0.171, 0.171, 0.171, 0.172, 0.172, 0.173, 0.173, 0.172, 0.171, 0.171, 0.170, 0.170,
       0.170, 0.171, 0.172, 0.172, 0.172, 0.172, 0.171, 0.170, 0.170, 0.169, 0.170, 0.171, 0.171, 0.171, 0.171, 0.171,
       0.170, 0.169, 0.169, 0.169, 0.169, 0.170, 0.170, 0.171, 0.171, 0.170, 0.170, 0.169, 0.168, 0.168, 0.168, 0.169,
       0.170, 0.170, 0.170, 0.170, 0.169, 0.168, 0.168, 0.167, 0.168, 0.169, 0.169, 0.170, 0.170, 0.169, 0.168, 0.168,
       0.167, 0.167, 0.168, 0.168, 0.169, 0.169, 0.169, 0.169, 0.168, 0.167, 0.167, 0.167, 0.167, 0.168, 0.168, 0.169,
       0.169, 0.168, 0.167, 0.167, 0.166, 0.167, 0.167, 0.168, 0.168, 0.168, 0.168, 0.168, 0.167, 0.166, 0.166, 0.166,
       0.166, 0.166, 0.167, 0.168, 0.167, 0.167, 0.166, 0.166, 0.165, 0.165, 0.166, 0.167, 0.167, 0.168, 0.168, 0.167,
       0.166, 0.166, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.167, 0.167, 0.166, 0.165, 0.165, 0.165, 0.166, 0.166,
       0.167, 0.167, 0.167, 0.166, 0.166, 0.165, 0.165, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.166, 0.166, 0.165,
       0.165, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166,
       0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.165, 0.166, 0.167, 0.167, 0.166, 0.165, 0.165, 0.164, 0.164,
       0.165, 0.165, 0.166, 0.166, 0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166, 0.166, 0.166,
       0.165, 0.165, 0.164, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166};
      unsigned dataArraySize = sizeof(tempCaTransient) / sizeof(double);
      extCaTransient.insert(extCaTransient.end(), &tempCaTransient[0], &tempCaTransient[dataArraySize]);
    }
# endif  // USE_EMT
#endif  // USE_FMT

    stretch = stretchdefault;

    if (outfile && strcmp(outfile, "NULL")) {
      fpo = fopen(outfile, "w");
      if (!fpo)
        throw kaBaseException("Opening file %s for output", outfile);
    }

#ifdef USE_EMT
    switch (cs) {
      case CS_Euler:
        csFunction = &vbElphyModel<CalcType>::Calc;
        break;
      case CS_RK2:
        csFunction = &vbElphyModel<CalcType>::CalcRungeKutta2;
        break;
      case CS_RK4:
        csFunction = &vbElphyModel<CalcType>::CalcRungeKutta4;
        break;
    }
#endif  // ifdef USE_EMT
  }  // init

  void calculate(const double t, const double tinc) {
    switch (sf) {
      case SF_VELOCITY: {
        if ((t >= stretchbegin) && (t <= stretchend))
          velocity = velocitydefault;
        else if ((t-stretchend >= restretch) && (t-2.*stretchend+stretchbegin <= restretch))
          velocity = velocitydefault;
        else
          velocity = 0;
        stretch += tinc*velocity;
      }
      break;

      case SF_RECT: {
        if ((t >= stretchbegin) && (t <= stretchend))
          stretch = stretchAmplitude;
        else
          stretch = stretchdefault;
      }
      break;

      default:
        break;
    }

#ifdef USE_EMT
    if (currentClampIndex < vClamp.size()) {
      if (t > currentEndTime) {
        currentClampIndex++;
        if (vClamp[currentClampIndex].type == 3)
          currentClampIndex++;
        currentEndTime += vClamp[currentClampIndex].t;

        // if (currentClampIndex<vClamp.size())
        //    cerr<<"t="<<t<<", setting to "<<vClamp[currentClampIndex].val<<" until "<<currentEndTime<<endl;
      }
    }
    if (currentClampIndex < vClamp.size()) {
      if (vClamp[currentClampIndex].type == 1) {
        Vm = vClamp[currentClampIndex].val;
      } else if (vClamp[currentClampIndex].type == 2) {
        pem->SetCai(vClamp[currentClampIndex].val);
      } else if (vClamp[currentClampIndex].type == 3) {
        Vm = vClamp[currentClampIndex-1].val;
        pem->SetCai(vClamp[currentClampIndex].val);
      } else if (vClamp[currentClampIndex].type == 4) {
        stretch = vClamp[currentClampIndex].val;
      }
    }
    CalcType i_external;
    if ((t >= textbegin) && (t < textend)) {
      i_external = (Vm >= thresVm ? 0. : amplitude);
      if (freq != HUGE) {
        i_external = sin(((t-textbegin)*2*M_PI*freq)+.000001)*i_external;
        if (monophasic && (i_external < 0.0) )
          i_external = 0.0;
      }
      if (t >= textbegin+textlen) {
        text      /= textdiv;
        textbegin += text;
      }
    } else {
      i_external = 0.0;
    }
    Vm += (pem->*csFunction)(tinc, Vm, i_external+i_internal, stretch, 1);
    if (std::isnan(Vm))
      throw kaBaseException("Vm is NaN");

# ifdef USE_FMT
    switch (fc) {
      case TroponinCoupling:
        Force = pfm->CalcTrop(tinc, stretch, velocity, pem->GetTCa()*1000.0);
        break;

      case CalciumCoupling: {
        CalcType Calcium = pem->GetCai()*1000.0;
        Force = pfm->Calc(tinc, stretch, velocity, Calcium);
      }
      break;

      case CalciumFeedBackCoupling: {
        CalcType Calcium = pem->GetCai()*1000.0;
        Force = pfm->Calc(tinc, stretch, velocity, Calcium);
        pem->SetCai(Calcium*.001);
      }
      break;
    }
# endif  // ifdef USE_FMT

#else  // ifdef USE_EMT
    switch (cf) {
      case CF_Const:
        Cai = CaConst;  // mikroM
        break;
      case CF_Default:
        Cai = 0.001;   // mikroM
        break;
      case CF_Rice:
        Cai = 0.1;     // mikroM
        break;
      case CF_Inc:
        Cai = 0.0;     // mikroM
        break;
      case CF_Coppini:
        Cai = 0.166; // mikroM
        break;
      default:
        Cai = 0.0;
        break;
    }

    if (t >= tCabegin) {
      const double tb = t-tCabegin;
      switch (cf) {
        case CF_Default: {
          const double tau_Ca = 0.06;  // s
          const double Ca_max = 3;  // mikroM
          Cai += (Ca_max-Cai)*tb/tau_Ca*exp(1-tb/tau_Ca);
          break;
        }

        case CF_Inc:
          Cai += 100*tb;  // Cai will be linear from Cai=0 at t=tbegin and Cai=100 at t+=1s
          break;

        case CF_Rice: {
          const double tau1   = 0.025; // s
          const double tau2   = .1;
          const double Ca_max = 1;  // mikroM Rice99 H1744
          Cai += 1.67*Ca_max*(1-exp(-tb/tau1))*exp(1-tb/tau2);
          break;
        }

        case CF_Coppini: {
          unsigned caIndex = round(tb / 0.001);
          caIndex = std::min(caIndex, (unsigned)1000);
          Cai = extCaTransient[caIndex];
          break;
        }

        case CF_Panerai: {
          const double Ca_max = 0.45;  // mikroM
          Cai = Ca_max*(1-exp(-200*tb*tb))*exp((tb < 0.3 ? 0 : -5)*(tb-0.3)*(tb-0.3));
          break;
        }

        case CF_Const:
          break;
      }  // switch
    }

    Force = pfm->Calc(tinc, stretch, velocity, Cai);
#endif  // ifdef USE_EMT
  }  // calculate

  void outputStatusLine(const StatusPrint spm) {
    vector<string> paraNames;

    paraNames.push_back("t");
#ifdef USE_EMT
    paraNames.push_back("V_m");
#endif  // ifdef USE_EMT
    switch (spm) {
      case StatusPrintBrief:
#ifdef USE_FMT
        paraNames.push_back("Tension");
#endif  // ifdef USE_FMT
        break;

      case StatusPrintBriefI:
#ifdef USE_EMT
        paraNames.push_back("I_internal");
#endif  // ifdef USE_EMT
        break;

      case StatusPrintShort:
        paraNames.push_back("Ca_i");
#ifdef USE_FMT
        paraNames.push_back("Tension");
#endif  // ifdef USE_FMT
        break;

      case StatusPrintNormal:
#ifdef USE_EMT
        pem->GetParameterNames(paraNames);
# ifdef USE_FMT

        //        paraNames.push_back("Tension");
        pfm->GetParameterNames(paraNames);
# endif  // ifdef USE_FMT
#else  // ifdef USE_EMT
        paraNames.push_back("Ca_i");
        paraNames.push_back("Tension");
        pfm->GetParameterNames(paraNames);
#endif  // ifdef USE_EMT
        break;

      case StatusPrintLong:
#ifdef USE_EMT
        pem->GetLongParameterNames(paraNames);
# ifdef USE_FMT
        paraNames.push_back("Tension");
        pfm->GetParameterNames(paraNames);
# endif  // ifdef USE_FMT
#else  // ifdef USE_EMT
        paraNames.push_back("Ca_i");
        paraNames.push_back("Tension");
        pfm->GetParameterNames(paraNames);
#endif  // ifdef USE_EMT
        break;

      case StatusPrintLongI:
#ifdef USE_EMT
        pem->GetLongParameterNames(paraNames);
        paraNames.push_back("I_internal");
#endif  // ifdef USE_EMT
        break;

      default:
        throw kaBaseException("Undefined print option");
        break;
    }  // switch

    ostringstream tostr(ostringstream::out);

    // tostr<<"#  "; lh326 deletet # on output file to import into other programms more easily
    tostr<<"";
    int length = paraNames.size();
    for (int i = 0; i < length; i++)
      tostr<<i+1<<":"<<paraNames[i]<<" ";
    tostr<<endl;
    tostr << '\0';
    output(tostr.str().c_str());
  }  // outputStatusLine

  inline void output(const char *buf) {
    if (fpo)
      fprintf(fpo, "%s", buf);
    else if (!outfile)
      cout<<buf;
  }

  void output(const StatusPrint spm, const double t) {
    ostringstream tostr(ostringstream::out);

    tostr.setf(ios_base::scientific);

    switch (spm) {
      case StatusPrintBrief:
        tostr << t <<' ';
#ifdef USE_EMT
        tostr << Vm << ' ';
#endif  // ifdef USE_EMT
#ifdef USE_FMT
        tostr << Force;
#endif  // ifdef USE_FMT
        break;

      case StatusPrintBriefI:
        tostr<<t;
#ifdef USE_EMT
        tostr<<' '<<Vm;
        tostr<<' '<<i_internal;
#endif  // ifdef USE_EMT
        break;

      case StatusPrintShort:
        tostr<<t<<' ';
#ifdef USE_EMT
        tostr<<Vm<<' ';
        Cai = pem->GetCai()*1000.0;
#endif  // ifdef USE_EMT
        tostr<<Cai << ' ';
#ifdef USE_FMT
        tostr << Force << ' ';
#endif  // ifdef USE_FMT
        break;

      case StatusPrintNormal:
#ifdef USE_EMT
        pem->Print(tostr, t, Vm);
# ifdef USE_FMT

        //        tostr << ' ' << Force << ' ';
        pfm->Print(tostr);
# endif  // ifdef USE_FMT
#else  // ifdef USE_EMT
        pfm->Print(tostr, t, Cai, Force);
#endif  // ifdef USE_EMT
        break;

      case StatusPrintLong:
#ifdef USE_EMT
        pem->LongPrint(tostr, t, Vm);
# ifdef USE_FMT

        // tostr << ' ' << Force << ' '; tg319: produced an extra space before the force values
        tostr << Force << ' ';
        pfm->Print(tostr);
# endif  // ifdef USE_FMT
#else  // ifdef USE_EMT
        pfm->Print(tostr, t, Cai, Force);
#endif  // ifdef USE_EMT
        break;

      case StatusPrintLongI:
#ifdef USE_EMT
        pem->LongPrint(tostr, t, Vm);
        tostr<<' '<<i_internal;
#endif  // ifdef USE_EMT
        break;

      default:
        throw kaBaseException("Undefined print option");
        break;
    }  // switch
    tostr << endl;
    tostr << '\0';

    output(tostr.str().c_str());
  }  // output

  void loadStatusFile() {
    FILE *fp = fopen(loadStatus, "rb");

    if (!fp)
      throw kaBaseException("Reading file %s", loadStatus);
#ifdef USE_EMT
    kaRead(&Vm, sizeof(CalcType), 1, fp);
    pem->ReadStatus(fp);
#endif  // ifdef USE_EMT
#ifdef USE_FMT
    pfm->ReadStatus(fp);
#endif  // ifdef USE_FMT
    fclose(fp);
  }

  void saveStatusFile() {
    FILE *fp = fopen(saveStatus, "wb");

    if (!fp)
      throw kaBaseException("Writing file %s", saveStatus);
#ifdef USE_EMT
    kaWrite(&Vm, sizeof(CalcType), 1, fp);
    pem->WriteStatus(fp);
#endif  // ifdef USE_EMT
#ifdef USE_FMT
    pfm->WriteStatus(fp);
#endif  // ifdef USE_FMT
    fclose(fp);
  }
};  // class CellParameters

#ifdef USE_EMT

//! Class handles coupling resistance between two cells
class CellCoupling {
 public:
  int cellnr1, cellnr2;
  double R, S;
  CellParameters *cell1, *cell2;
};
#endif  // ifdef USE_EMT
