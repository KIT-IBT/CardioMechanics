/*! \file APAnalysis.cpp
   \brief Functions to analyse files with numerical data ordered
   in rows (lines) and columns separated by spaces or tabs.

   \author gs, IBT - Universität Karlsruhe (TH), fs, CVRTI - University
   of Utah
 */

#include <kaExceptions.h>
#include <ExpFit.h>

#include <cmath>
#include <limits>

#undef HUGE
#define HUGE std::numeric_limits<double>::max()

char formatx[256];  //!< Format string for parsing of lines with scanf
char formaty[1024];
double range0 = -HUGE, range1 = HUGE;

//! Set formal string to get values in 1. and specified column
void setFormat(int col) {
  strcpy(formatx, "%lf");

  formaty[0] = 0;
  for (int i = 1; i < col; i++)
    strcat(formaty, "%*s ");
  strcat(formaty, "%lf");
}

//! Read and check line, skip comments
inline bool getLine(ifstream &in, double &a, double &b) {
  char instr[4096];
  bool ok;

  do {
    do {
      in.getline(instr, sizeof(instr));
      if (in.eof())
        return false;
    } while (instr[0] == '#');

    int rc = sscanf(instr, formatx, &a);
    if (!rc)
      a = NAN;
    rc = sscanf(instr, formaty, &b);
    ok = (rc == 1);
  } while (ok && (a < range0 || a > range1));

  return ok;
}

inline void getFirstLine(char *infile, double &v, double &t) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);

  if (!getLine(in, t, v))
    t = v = -HUGE;
}

inline void getLastLine(char *infile, double &v, double &t) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  double tv, tt;

  if (!getLine(in, tt, tv)) {
    v = t = -HUGE;
  } else {
    do {
      t = tt;
      v = tv;
    } while (getLine(in, tt, tv));
  }
}

inline void getAtLine(char *infile, double &v, double t) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  double tv, tt;
  v = -HUGE;

  while (getLine(in, tt, tv))
    if (tt == t) {  // ?
      v = tv;
      break;
    }
}

inline void getMinValue(char *infile, double &minvalue, double &tmv) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  double v, t;

  if (!getLine(in, tmv, minvalue)) {
    minvalue = tmv = -HUGE;
  } else {
    while (getLine(in, t, v))
      if (v < minvalue) {
        minvalue = v;
        tmv      = t;
      }
  }
}

inline void getMaxValue(char *infile, double &maxvalue, double &tmv) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  double v, t;

  if (!getLine(in, tmv, maxvalue)) {
    maxvalue = tmv = -HUGE;
  } else {
    while (getLine(in, t, v))
      if (v > maxvalue) {
        maxvalue = v;
        tmv      = t;
      }
  }
}

inline void getMinMaxValue(char *infile, double &minvalue, double &maxvalue, double &mintime, double &maxtime) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  double v, t;

  if (!getLine(in, mintime, minvalue)) {
    minvalue = maxvalue = mintime = maxtime = -HUGE;
  } else {
    maxvalue = minvalue;
    maxtime  = mintime;
    while (getLine(in, t, v)) {
      if (v < minvalue) {
        minvalue = v;
        mintime  = t;
      }
      if (v > maxvalue) {
        maxvalue = v;
        maxtime  = t;
      }
    }
  }
}

inline double getTime(char *infile, double threshold, bool up = true) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);

  double lowervalue  = 0.0;
  double highervalue = 0.0;
  double timelow     = 0.0;
  double timehigh    = 0.0;

  if (!getLine(in, timelow, lowervalue))
    return -HUGE;

  if (up) {
    if (lowervalue > threshold)
      return timelow;
  }

  while (getLine(in, timehigh, highervalue)) {
    if (up && (highervalue > threshold) && (lowervalue <= threshold))
      break;
    if (!up && (highervalue < threshold) && (lowervalue >= threshold))
      break;
    timelow    = timehigh;
    lowervalue = highervalue;
  }
  return in.eof() ? -HUGE : timelow+(threshold-lowervalue)/(highervalue-lowervalue)*(timehigh-timelow);
}  // getTime

inline int getCnt(char *infile, double threshold, bool up = true) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);

  double value1 = 0.0;
  double value2 = 0.0;
  double t;
  int cnt = 0;

  if (!getLine(in, t, value2))
    return 0;

  if (up) {
    if (value2 > threshold)
      cnt++;
  }

  while (getLine(in, t, value1)) {
    if (up && (value1 > threshold) && (value2 <= threshold))
      cnt++;
    if (!up && (value2 < threshold) && (value2 >= threshold))
      cnt++;
    value2 = value1;
  }
  return cnt;
}  // getCnt

inline void getFitSub(char *infile, double *p, int funcNum) {
  ifstream in(infile);

  if (!in.is_open())
    throw kaBaseException("Can't open file %s\n", infile);
  const int maxline = 10000;
  double   *x = new double[maxline], *y = new double[maxline];

  int num;
  for (num = 1; getLine(in, x[num], y[num]); num++) {
    //     printf("%lf %lf\n", x[num], y[num]);
    if (num >= maxline)
      throw kaBaseException("Too many lines in file %s\n", infile);
  }
  num--;

  //  cout<<"Number of read x-y pairs: "<<num<<endl;

  const int maxdim = 10;

  double **alpha = new double *[maxdim];
  for (int i = 0; i < maxdim; i++)
    alpha[i] = new double[maxdim];

  double **covar = new double *[maxdim];
  for (int i = 0; i < maxdim; i++)
    covar[i] = new double[maxdim];

  int ma, ia[maxdim];
  double a[maxdim];

  double avgx = 0, avgy = 0, x1 = x[1], miny = HUGE, maxy = -HUGE;
  for (int i = 1; i <= num; i++) {
    //    if (funcNum==1 || funcNum==0) x[i]-=x1;
    if ((funcNum == 2) || (funcNum == 3)) {
      if (x[i] < 0)
        throw kaBaseException("getFitSub: negative x-values are not allowed");
      if (x[i] == 0)
        x[i] = 1e-15;
    }
    avgx += x[i];
    avgy += y[i];
    if (y[i] > maxy)
      maxy = y[i];
    if (y[i] < miny)
      miny = y[i];
  }
  avgx /= num;
  avgy /= num;

  switch (funcNum) {
    case 0:
      ma    = 3;
      ia[1] = 1; ia[2] = 1; ia[3] = 1;
      a[1]  = maxy; a[2] = (x[num-1]-x[1])/2; a[3] = -(x[num-1]-x[1])/(y[num-1]-y[1]);

      //    a[1]=-200; a[2]=0; a[3]=0.03;
      //    cerr << a[1] << '\t' << a[2] << '\t' << a[3] << endl;
      break;
    case 1:
      ma    = 4;
      ia[1] = ia[2] = 1; ia[3] = ia[4] = 0;
      a[1]  = 1.; a[2] = 1.; a[3] = x[1]; a[4] = 0.;
      break;
    case 2:
      ma    = 3;
      ia[1] = ia[2] = ia[3] = 1;
      a[1]  = 1.; a[2] = avgx; a[3] = maxy;

      //    printf(" %lf %lf %lf\n", a[1], a[2], a[3]);
      break;
    case 3:
      ma    = 4;
      ia[1] = ia[2] = ia[3] = 1;
      ia[4] = 1;
      a[1]  = 1.; a[2] = avgx; a[3] = maxy; a[4] = 0;

      //    printf(" %lf %lf %lf\n", a[1], a[2], a[3]);
      break;
    default:
      throw kaBaseException("getFitSub: no usable function number.");
  }  // switch

  double *sig = new double[maxline];
  for (int i = 0; i < num+1; i++)
    sig[i] = 1/sqrt((double)num);

  double alambda = -1., chisq = HUGE;
  for (int cnt = 0; cnt < 1000; cnt++) {
    mrqmin(x, y, sig, num, a, ia, ma, covar, alpha, &chisq, funcNum, &alambda);

    //    printf("%d %le %le\t%lf %lf %lf %lf\n", cnt, chisq, alambda, a[1], a[2], a[3], a[4]);
    if ((chisq < 1e-7) || (alambda > 1e10))
      break;
  }
  alambda = 0;
  mrqmin(x, y, sig, num, a, ia, ma, covar, alpha, &chisq, funcNum, &alambda);
  for (int i = 0; i < ma; i++)
    p[i] = a[i+1];

  for (int i = 0; i < maxdim; i++)
    delete alpha[i];
  delete[] alpha;
  for (int i = 0; i < maxdim; i++)
    delete covar[i];

  delete[] covar;
  delete[] sig;
  delete[] y;
  delete[] x;
}  // getFitSub

int main(int argc, char *argv[]) {
  if (argc <= 2) {
    cerr << argv[0] << " <file1>" << endl;
    cerr << "\t[-col <number>][-range <t0> <t1>]" << endl;
    cerr << "\t[-vel <file2> <distance>]" << endl;
    cerr << "\t[-max][-min][-ext][-amplitude]" << endl;
    cerr << "\t[-mean][-integral][-stddev][-stderr]" << endl;
    cerr << "\t[-maxdvdt][-maxdvdt2][-mindvdt][-mindvdt2]" << endl;
    cerr << "\t[-dynintegrate][-dynmin][-dynmax]" << endl;
    cerr << "\t[-timemax][-timemin][-timeext][-timemaxdvdt][-timemaxdvdt2][-timemindvdt][-timemindvdt2]" << endl;
    cerr << "\t[-apd <percentage>][-apd0 <percentage>]" << endl;
    cerr << "\t[-timeupstroke <val>][-timedownstroke <val>]" << endl;
    cerr << "\t[-cntupstroke <val>][-cntdownstroke <val>]" << endl;
    cerr << "\t[-last][-first]" << endl;
    cerr << "\t[-at <time> [-at <time>]...]" << endl;
    cerr << "\t[-atValueUp <targetValue> [-atValueUp <targetValue>]...]" << endl;
    cerr << "\t[-atValueDown <targetValue> [-atValueDown <targetValue>]...]" << endl;
    cerr << "\t[-boltzmann|-hill|-hillOffset|-exp]" << endl;
    cerr << "\t[-tau]" << endl;
    cerr << "\t[-normMin|-normMax|-normMinMax]" << endl;
    cerr << "\t[-mult <file2 or value>][-add <file2 or value>]" << endl;
    exit(-1);
  }
  try {
    char *infile1 = argv[1];
    char *infile2;
    double distance, percentage, threshold, number;
    vector<double> attime;
    vector<double> atValue;
    enum atValueDirectionEnum { Up, Down };
    vector<int> atValueDirection;

    int col = 2;

    enum analysisFunction {
      analysisVelocity,
      analysisMin,
      analysisMax,
      analysisExt,
      analysisAmplitude,
      analysisMean,
      analysisStdDev,
      analysisStdErr,
      analysisIntegral,
      analysisTimeMin,
      analysisTimeExt,
      analysisTimeMax,
      analysisMaxDvDt,
      analysisMaxDvDt2,
      analysisMinDvDt,
      analysisMinDvDt2,
      analysisTimeMaxDvDt,
      analysisTimeMaxDvDt2,
      analysisTimeMinDvDt,
      analysisTimeMinDvDt2,
      analysisAPD,
      analysisAPD0,
      analysisDynintegrate,
      analysisDynMin,
      analysisDynMax,
      analysisTimeUpstroke,
      analysisTimeDownstroke,
      analysisCntUpstroke,
      analysisCntDownstroke,
      analysisLast,
      analysisFirst,
      analysisAt,
      analysisBoltzmann,
      analysisHill,
      analysisHillOffset,
      analysisExp,
      analysisTau,
      analysisNormMin,
      analysisNormMax,
      analysisNormMinMax,
      analysisMult,
      analysisAdd,
      analysisAtValue,
    }
    aF;
    aF = analysisMaxDvDt;

    for (int i = 2; i < argc; i++)
      if (!strcasecmp(argv[i], "-vel") && (i < argc-2)) {
        aF       = analysisVelocity;
        infile2  = argv[++i];
        distance = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-maxdvdt")) {
        aF = analysisMaxDvDt;
      } else if (!strcasecmp(argv[i], "-maxdvdt2")) {
        aF = analysisMaxDvDt2;
      } else if (!strcasecmp(argv[i], "-timemaxdvdt")) {
        aF = analysisTimeMaxDvDt;
      } else if (!strcasecmp(argv[i], "-timemaxdvdt2")) {
        aF = analysisTimeMaxDvDt2;
      } else if (!strcasecmp(argv[i], "-mindvdt")) {
        aF = analysisMinDvDt;
      } else if (!strcasecmp(argv[i], "-mindvdt2")) {
        aF = analysisMinDvDt2;
      } else if (!strcasecmp(argv[i], "-timemindvdt")) {
        aF = analysisTimeMinDvDt;
      } else if (!strcasecmp(argv[i], "-timemindvdt2")) {
        aF = analysisTimeMinDvDt2;
      } else if (!strcasecmp(argv[i], "-max")) {
        aF = analysisMax;
      } else if (!strcasecmp(argv[i], "-min")) {
        aF = analysisMin;
      } else if (!strcasecmp(argv[i], "-mean")) {
        aF = analysisMean;
      } else if (!strcasecmp(argv[i], "-stddev")) {
        aF = analysisStdDev;
      } else if (!strcasecmp(argv[i], "-stderr")) {
        aF = analysisStdErr;
      } else if (!strcasecmp(argv[i], "-ext")) {
        aF = analysisExt;
      } else if (!strcasecmp(argv[i], "-amplitude")) {
        aF = analysisAmplitude;
      } else if (!strcasecmp(argv[i], "-integral")) {
        aF = analysisIntegral;
      } else if (!strcasecmp(argv[i], "-timemax")) {
        aF = analysisTimeMax;
      } else if (!strcasecmp(argv[i], "-timemin")) {
        aF = analysisTimeMin;
      } else if (!strcasecmp(argv[i], "-timeext")) {
        aF = analysisTimeExt;
      } else if (!strcasecmp(argv[i], "-dynintegrate")) {
        aF = analysisDynintegrate;
      } else if (!strcasecmp(argv[i], "-dynmin")) {
        aF = analysisDynMin;
      } else if (!strcasecmp(argv[i], "-dynmax")) {
        aF = analysisDynMax;
      } else if (!strcasecmp(argv[i], "-apd") && (i < argc-1) ) {
        aF         = analysisAPD;
        percentage = atof(argv[++i])/100.;
      } else if (!strcasecmp(argv[i], "-apd0") && (i < argc-1)) {
        aF         = analysisAPD0;
        percentage = atof(argv[++i])/100.;
      } else if (!strcasecmp(argv[i], "-timeupstroke") && (i < argc-1)) {
        aF        = analysisTimeUpstroke;
        threshold = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-timedownstroke") && (i < argc-1)) {
        aF        = analysisTimeDownstroke;
        threshold = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-cntupstroke") && (i < argc-1)) {
        aF        = analysisCntUpstroke;
        threshold = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-cntdownstroke") && (i < argc-1)) {
        aF        = analysisCntDownstroke;
        threshold = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-col") && (i < argc-1)) {
        col = atoi(argv[++i]);
        if (col < 1)
          throw kaBaseException("Column number must be larger 0");
      } else if (!strcasecmp(argv[i], "-range") && (i < argc-2)) {
        range0 = atof(argv[++i]);
        range1 = atof(argv[++i]);
      } else if (!strcasecmp(argv[i], "-first")) {
        aF = analysisFirst;
      } else if (!strcasecmp(argv[i], "-last")) {
        aF = analysisLast;
      } else if (!strcasecmp(argv[i], "-at") && (i < argc-1) ) {
        aF = analysisAt;
        attime.push_back(atof(argv[++i]));
      } else if (!strcasecmp(argv[i], "-atValueUp") && (i < argc-1)) {
        aF = analysisAtValue;
        atValue.push_back(atof(argv[++i]));
        atValueDirection.push_back(Up);
      } else if (!strcasecmp(argv[i], "-atValueDown") && (i < argc-1)) {
        aF = analysisAtValue;
        atValue.push_back(atof(argv[++i]));
        atValueDirection.push_back(Down);
      } else if (!strcasecmp(argv[i], "-boltzmann")) {
        aF = analysisBoltzmann;
      } else if (!strcasecmp(argv[i], "-exp")) {
        aF = analysisExp;
      } else if (!strcasecmp(argv[i], "-hill")) {
        aF = analysisHill;
      } else if (!strcasecmp(argv[i], "-hillOffset")) {
        aF = analysisHillOffset;
      } else if (!strcasecmp(argv[i], "-tau")) {
        aF = analysisTau;
      } else if (!strcasecmp(argv[i], "-normMin")) {
        aF = analysisNormMin;
      } else if (!strcasecmp(argv[i], "-normMax")) {
        aF = analysisNormMax;
      } else if (!strcasecmp(argv[i], "-normMinMax")) {
        aF = analysisNormMinMax;
      } else if (!strcasecmp(argv[i], "-add") && (i < argc-1) ) {
        aF = analysisAdd;
        char *isNotNumber = NULL;
        number = strtod(argv[++i], &isNotNumber);
        if (argv[i] == isNotNumber) {
          infile2 = argv[i];
          number  = HUGE;
        }
      } else if (!strcasecmp(argv[i], "-mult") && (i < argc-1) ) {
        aF = analysisMult;
        char *isNotNumber = NULL;
        number = strtod(argv[++i], &isNotNumber);
        if (argv[i] == isNotNumber) {
          infile2 = argv[i];
          number  = HUGE;
        }
      } else {
        throw kaBaseException("Unknown option %s", argv[i]);
      }

    cout.precision(10);

    setFormat(col);
    switch (aF) {
      case analysisMaxDvDt:
      case analysisTimeMaxDvDt: {
        double lowervalue, highervalue, timelow, timehigh;
        double maxdvdtvalue = -HUGE, timemaxdvdtvalue = -HUGE;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, timelow, lowervalue);
        while (getLine(in, timehigh, highervalue)) {
          double dV = highervalue-lowervalue;
          double dt = 1.0/(timehigh-timelow);
          if (dV*dt > maxdvdtvalue) {
            maxdvdtvalue     = dV*dt;
            timemaxdvdtvalue = (timehigh+timelow)/2.0;
          }
          timelow    = timehigh;
          lowervalue = highervalue;
        }
        cout << (aF == analysisMaxDvDt ? maxdvdtvalue : timemaxdvdtvalue) << endl;
        break;
      }

      case analysisMaxDvDt2:
      case analysisTimeMaxDvDt2: {
        double v1, v2, v3, t1, t2, t3;
        double maxdvdtvalue = -HUGE, timemaxdvdtvalue = -HUGE;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, t1, v1);
        getLine(in, t2, v2);
        while (getLine(in, t3, v3)) {
          double dV = v1-2.*v2+v3;
          double dt = 1./(t3-t1);
          if (dV*dt > maxdvdtvalue) {
            maxdvdtvalue     = dV*dt;
            timemaxdvdtvalue = t2;
          }
          t1 = t2; t2 = t3;
          v1 = v2; v2 = v3;
        }
        cout << (aF == analysisMaxDvDt2 ? maxdvdtvalue : timemaxdvdtvalue) << endl;
        break;
      }

      case analysisMinDvDt:
      case analysisTimeMinDvDt: {
        double lowervalue, highervalue, timelow, timehigh;
        double mindvdtvalue = HUGE, timemindvdtvalue = HUGE;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, timelow, lowervalue);
        while (getLine(in, timehigh, highervalue)) {
          double dV = highervalue-lowervalue;
          double dt = 1.0/(timehigh-timelow);
          if (dV*dt < mindvdtvalue) {
            mindvdtvalue     = dV*dt;
            timemindvdtvalue = (timehigh+timelow)/2.0;
          }
          timelow    = timehigh;
          lowervalue = highervalue;
        }
        cout << (aF == analysisMinDvDt ? mindvdtvalue : timemindvdtvalue) << endl;
        break;
      }

      case analysisMinDvDt2:
      case analysisTimeMinDvDt2: {
        double v1, v2, v3, t1, t2, t3;
        double mindvdtvalue = HUGE, timemindvdtvalue = HUGE;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, t1, v1);
        getLine(in, t2, v2);
        while (getLine(in, t3, v3)) {
          double dV = v1-2.*v2+v3;
          double dt = 1./(t3-t1);
          if (dV*dt < mindvdtvalue) {
            mindvdtvalue     = dV*dt;
            timemindvdtvalue = t2;
          }
          t1 = t2; t2 = t3;
          v1 = v2; v2 = v3;
        }
        cout << (aF == analysisMinDvDt2 ? mindvdtvalue : timemindvdtvalue) << endl;
        break;
      }
      case analysisVelocity: {
        double zerotime1 = getTime(infile1, -0.03);
        double zerotime2 = getTime(infile2, -0.03);
        cout << distance/(zerotime2-zerotime1) << endl;
        break;
      }

      case analysisMax: {
        double value, time;
        getMaxValue(infile1, value, time);
        cout << value << endl;
        break;
      }
      case analysisMin: {
        double value, time;
        getMinValue(infile1, value, time);
        cout << value << endl;
        break;
      }
      case analysisMean: {
        double t, v, sum = 0;
        int cnt = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v)) {
          sum += v;
          cnt++;
        }
        cout << sum/cnt << endl;
        break;
      }
      case analysisStdDev: {
        double t, v, sum = 0, sum2 = 0;
        int cnt = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v)) {
          sum  += v;
          sum2 += v*v;
          cnt++;
        }
        cout << sqrt(cnt*sum2-sum*sum)/cnt << endl;
        break;
      }
      case analysisStdErr: {
        double t, v, sum = 0, sum2 = 0;
        int cnt = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v)) {
          sum  += v;
          sum2 += v*v;
          cnt++;
        }
        cout << sqrt(cnt*sum2-sum*sum)/cnt/sqrt(cnt) << endl;
        break;
      }
      case analysisExt: {
        double maxvalue, minvalue, dummy;
        getMinMaxValue(infile1, minvalue, maxvalue, dummy, dummy);
        double extvalue;
        if (fabs(minvalue) > fabs(maxvalue))
          extvalue = minvalue;
        else
          extvalue = maxvalue;

        cout << extvalue << endl;
        break;
      }
      case analysisAmplitude: {
        double maxvalue, minvalue, dummy;
        getMinMaxValue(infile1, minvalue, maxvalue, dummy, dummy);
        cout << maxvalue-minvalue << endl;
        break;
      }
      case analysisIntegral: {
        double t, v, sum = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v))
          sum += v;
        cout << ' ' << sum << endl;
        break;
      }
      case analysisTimeMax: {
        double maxvalue, maxtime;
        getMaxValue(infile1, maxvalue, maxtime);
        cout << maxtime << endl;
        break;
      }
      case analysisTimeMin: {
        double minvalue, mintime;
        getMinValue(infile1, minvalue, mintime);
        cout << mintime << endl;
        break;
      }
      case analysisTimeExt: {
        double maxvalue, minvalue, maxtime, mintime;
        getMinMaxValue(infile1, minvalue, maxvalue, mintime, maxtime);
        double exttime;
        if (fabs(minvalue) > fabs(maxvalue))
          exttime = mintime;
        else
          exttime = maxtime;

        cout << exttime << endl;
        break;
      }
      case analysisAPD: {
        double maxvalue, minvalue, dummy;
        getMinMaxValue(infile1, minvalue, maxvalue, dummy, dummy);
        double starttime = getTime(infile1, maxvalue-percentage*(maxvalue-minvalue));
        double endtime   = getTime(infile1, maxvalue-percentage*(maxvalue-minvalue), false);
        cout << endtime-starttime << endl;
        break;
      }
      case analysisAPD0: {
        double maxvalue, minvalue, dummy;
        getMinMaxValue(infile1, minvalue, maxvalue, dummy, dummy);
        double endtime = getTime(infile1, maxvalue-percentage*(maxvalue-minvalue), false);
        cout << endtime << endl;
        break;
      }
      case analysisTimeUpstroke: {
        double ltime = getTime(infile1, threshold);
        cout << ltime << endl;
        break;
      }
      case analysisTimeDownstroke: {
        double stime = getTime(infile1, threshold, false);
        cout << stime << endl;
        break;
      }
      case analysisCntUpstroke: {
        int cnt = getCnt(infile1, threshold);
        cout << cnt << endl;
        break;
      }
      case analysisCntDownstroke: {
        int cnt = getCnt(infile1, threshold, false);
        cout << cnt << endl;
        break;
      }
      case analysisDynintegrate: {
        double t, v, sum = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v)) {
          sum += v;
          cout << t << ' ' << sum << endl;
        }
        break;
      }
      case analysisDynMin: {
        double t, v, min = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, t, min);
        cout << t << ' ' << min << endl;
        while (getLine(in, t, v)) {
          if (v < min)
            min = v;
          cout << t << ' ' << min << endl;
        }
        break;
      }
      case analysisDynMax: {
        double t, v, max = 0;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        getLine(in, t, max);
        cout << t << ' ' << max << endl;
        while (getLine(in, t, v)) {
          if (v > max)
            max = v;
          cout << t << ' ' << max << endl;
        }
        break;
      }
      case analysisFirst: {
        double t, v;
        getFirstLine(infile1, v, t);
        cout << v << endl;
        break;
      }
      case analysisLast: {
        double t, v;
        getLastLine(infile1, v, t);
        cout << v << endl;
        break;
      }
      case analysisAt: {
        double v;
        for (int i = 0; i < attime.size(); i++) {
          getAtLine(infile1, v, attime[i]);
          cout << v << ' ';
        }
        cout << endl;
        break;
      }
      case analysisAtValue: {
        for (int i = 0; i < atValue.size(); i++) {
          double time = getTime(infile1, atValue[i], (atValueDirection[i] == Up ? true : false));
          cout << time << ' ';
        }
        cout << endl;
        break;
      }
      case analysisBoltzmann: {
        double p[3];
        getFitSub(infile1, p, 0);
        cout << p[0] << ' ' << p[1] << ' ' << p[2] << endl;
        break;
      }
      case analysisExp: {
        double p[4];
        getFitSub(infile1, p, 1);
        cout << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << p[3] << endl;
        break;
      }
      case analysisHill: {
        double p[3];
        getFitSub(infile1, p, 2);
        cout << p[0] << ' ' << p[1] << ' ' << p[2] << endl;
        break;
      }
      case analysisHillOffset: {
        double p[4];
        getFitSub(infile1, p, 3);
        cout << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << p[3] << endl;
        break;
      }
      case analysisTau: {
        double maxvalue, minvalue, maxtime, mintime;
        getMinMaxValue(infile1, minvalue, maxvalue, mintime, maxtime);

        //         cerr << minvalue <<  '\t' << maxvalue << endl;
        //         cerr << mintime <<  '\t' << maxtime << endl;
        double extvalue, exttime;
        if (fabs(minvalue) > fabs(maxvalue)) {
          extvalue = minvalue;
          range0   = exttime = mintime;
        } else {
          extvalue = maxvalue;
          range0   = exttime = maxtime;
        }
        double tau = getTime(infile1, extvalue/exp(1.0), extvalue < 0.);

        cout << tau-exttime << endl;
        break;
      }
      case analysisNormMin: {
        double minvalue, mintime;
        getMinValue(infile1, minvalue, mintime);
        if (!minvalue)
          minvalue = 1.;

        double t, v;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v))
          cout << t << ' ' << (v-minvalue)/minvalue << endl;
        break;
      }
      case analysisNormMax: {
        double maxvalue, maxtime;
        getMaxValue(infile1, maxvalue, maxtime);
        if (!maxvalue)
          maxvalue = 1.;

        double t, v;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v))
          cout << t << ' ' << v/maxvalue << endl;
        break;
      }
      case analysisNormMinMax: {
        double maxvalue, minvalue, dummy;
        getMinMaxValue(infile1, minvalue, maxvalue, dummy, dummy);
        double amplitude = maxvalue-minvalue;
        if (!amplitude)
          amplitude = 1.;

        double t, v;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        while (getLine(in, t, v))
          cout << t << ' ' << (v-minvalue)/amplitude << endl;
        break;
      }
      case analysisAdd: {
        double t1, v1, t2, v2;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        if (number == HUGE) {
          ifstream in2(infile2);
          if (!in2.is_open())
            throw kaBaseException("Can't open file %s\n", infile2);
          while (getLine(in, t1, v1)) {
            getLine(in2, t2, v2);
            cout << t1 << ' ' << v1+v2 << endl;
          }
        } else {
          while (getLine(in, t1, v1))
            cout << t1 << ' ' << v1+number << endl;
        }
        break;
      }
      case analysisMult: {
        double t1, v1, t2, v2;
        ifstream in(infile1);
        if (!in.is_open())
          throw kaBaseException("Can't open file %s\n", infile1);
        if (number == HUGE) {
          ifstream in2(infile2);
          if (!in2.is_open())
            throw kaBaseException("Can't open file %s\n", infile2);
          while (getLine(in, t1, v1)) {
            getLine(in2, t2, v2);
            cout << t1 << ' ' << v1*v2 << endl;
          }
        } else {
          while (getLine(in, t1, v1))
            cout << t1 << ' ' << v1*number << endl;
        }
        break;
      }
      default:
        throw kaBaseException("Unknown type of analysis");
    }  // switch
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
