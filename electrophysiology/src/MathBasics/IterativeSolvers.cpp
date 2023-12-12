/*
 * File: IterativeSolvers.cpp
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


#include <IterativeSolvers.h>

#include <sys/types.h>
#include <sys/times.h>

IterativeSolvers::IterativeSolvers() {
  tolerance = 1e-15;
  max_iter  = 1000;
  alfa      = 0.00005;
}

IterativeSolvers::~IterativeSolvers() {}

void IterativeSolvers::SetTolerance(const double &t) {tolerance = t;}

void IterativeSolvers::SetMaxIter(const int &i) {max_iter = i;}

void IterativeSolvers::SetAlfa(const double &a) {alfa = a;}

int IterativeSolvers::CG(const CMatrix &A, const std::valarray<double> &b, std::valarray<double> &x, int pm) {
  struct tms tmsbuf;
  double usrtime;

  if (max_iter == 0) {
    info_stream(2)<<"max_iter == 0\n";
    return 1;
  }

  int retcode;  // return code
  int MAXN = A.MAXN;
  CMatrix *CH;  // for Cholesky factorisation
  int (*chol)(CMatrix &);

  //  CholFactSolve * pCFS;


  info_stream(2)<<"CG: MAXN="<<MAXN<<'\n';

  if (pm == -1) {
    info_stream(2)<<"CG:no preconditioning\n";
  } else {
    if (pm == 0) {
      info_stream(2)<<"CG:diagonal preconditioning\n";
    } else {
      if (pm == 1) {
        info_stream(2)<<"CG: JPICC (improved incomplete Cholesky preconditioning)\n";
        chol = JPICC;
      } else {
        info_stream(2)<<"CG: ISTDIC (standard no-fill incomplete Cholesky preconditioning)\n";
        chol = ISTDIC;
      }
      CH = new CMatrix(A);
      double factor = 1.0;
      double a      = 0.0;
      int retcode;

      times(&tmsbuf);
      usrtime = -(double)tmsbuf.tms_utime/(double)CLK_TCK;

      while (retcode = chol(*CH)) {
        info_stream(2)<<"CG:Cholesky failed for alfa="<<a<<", return code:"<<retcode<<'\n';
        factor *= 2.0;
        a       = alfa*factor;

        // restore CH and try again
        *CH        = A;
        CH->CDIAG += a;
      }
      info_stream(2)<<"CG: Cholesky OK, alfa="<<a<<'\n';

      times(&tmsbuf);
      usrtime += (double)tmsbuf.tms_utime/(double)CLK_TCK;
      cout<<"Time for incomplete Cholesky [sec]:"<<setprecision(2)<<fixed<<usrtime<<'\n';
    }
  }

  valarray<double> res(0.0, MAXN);
  valarray<double> p(0.0, MAXN);
  valarray<double> q(0.0, MAXN);
  valarray<double> z(0.0, MAXN);

  double ro0  = 0;
  double ro1  = 0;
  double beta = 0;

  A.SymmMatrixVectorProduct(x, p);  // p is used to save memory

  for (int i = 0; i < MAXN; i++) {
    res[i] = b[i]-p[i];
    p[i]   = 0.0;
  }

  double norm_b = 0.0;
  for (int i = 0; i < MAXN; i++) norm_b += b[i]*b[i];
  norm_b = sqrt(norm_b);

  double rlimit = tolerance*norm_b;

  if (norm_b == 0.0) {
    info_stream(2) << " CG: Norm of B is zero! Solution is set to zero.";
    x = 0;
    return 0;
  }


  bool NextIteration = true;
  int  Iteration     = 1;

  while (NextIteration) {
    if (pm == -1)
      z = res;
    else if (pm == 0)
      z = res/A.CDIAG;
    else
      CholFactSolve(*CH, res, z); // (*pCFS)(res, z);

    ro0 = ro1;
    ro1 = 0;

    // for(int i=0; i<MAXN; i++)  ro1=ro1+res[i]*z[i];
    ro1 = inner_product(&res[0], &res[MAXN], &z[0], 0.0);

    if (Iteration == 1) {
      copy(&z[0], &z[MAXN], &p[0]);  // for(int i=0; i<MAXN; i++) p[i]=z[i];
    } else {
      beta = ro1/ro0;
      Saypx(&p[0], &p[MAXN], &z[0], beta);

      // for(int i=0; i<MAXN; i++) p[i]=z[i]+beta*p[i];
    }


    A.SymmMatrixVectorProduct(p, q);
    double denom = 0;

    // for(int i=0; i<MAXN; i++) denom+=p[i]*q[i];
    denom = inner_product(&p[0], &p[MAXN], &q[0], 0.0);
    double aa = ro1/denom;

    Saxpy(&x[0], &x[MAXN], &p[0], aa);
    Saxpy(&res[0], &res[MAXN], &q[0], -aa);

    double norm_r = 0;
    norm_r = inner_product(&res[0], &res[MAXN], &res[0], 0.0);
    norm_r = sqrt(norm_r);


    info_stream(2)<<setprecision(8)<<scientific<<" CG: Iteration = "<< Iteration <<" ||res|| = "<<norm_r;
    info_stream(2)<<setprecision(8)<<scientific<<" rlimit = "<<rlimit<<'\n';

    if (norm_r < rlimit) {NextIteration = false; retcode = 0;}
    if (Iteration == max_iter) {
      info_stream(2) << "CG: Max iteration number reached. Exiting...\n";
      NextIteration = false;
      retcode       = 1;
    }
    Iteration++;
  }

  info_stream(2)<<"CG: completed \n";
  if (CH)
    delete CH;
  return retcode;
}  // IterativeSolvers::CG

// int IterativeSolvers::BiCG(CMatrix &A, const std::valarray<double> &b, std::valarray<double> &x, char pm){
//  int retcode;
//  int MAXN=A.MAXN;
//  CMatrix CH(0,0); //for Cholesky factorisation
//  int (*chol)(CMatrix &);
//  info_stream(2)<<"BiCGSTAB: MAXN="<<MAXN<<'\n';
//   if(pm==0){
//     info_stream(2)<<"BiCGSTAB:diagonal preconditioning\n";
//   }
//   else{
//     if(pm==1 ){
//       info_stream(2)<<"BiCGSTAB: JPICC (improved incomplete Cholesky preconditioning)\n";
//       chol=JPICC;
//     }
//     else{
//        info_stream(2)<<"BiCGSTAB: ISTDIC (standard no-fill incomplete Cholesky preconditioning)\n";
//        chol=ISTDIC;
//     }
//     CH = A;
//     double factor=1.0;
//     double a=0;
//     while(chol(*CH)!=0){
//       info_stream(2)<<"CG:Cholesky failed for alfa="<<a<<'\n';
//       factor*=2.0;
//       a=alfa*factor;
//       //restore CH and try again
//       CH=A;
//       CH->CDIAG+=a;
//     }
//     info_stream(2)<<"BiCGSTAB: Cholesky OK, alfa="<<a<<'\n';
//   }

//   valarray<double> res(0.0, MAXN);
//   valarray<double> res1(0.0, MAXN);
//   valarray<double> p(0.0, MAXN);
//   valarray<double> p1(0.0, MAXN);

//   double ro0=0;
//   double ro1=0;
//   double beta=0;

//   SymmMatrixVectorProduct(A, x, p);//p is used to save memory
//   for(int i=0; i<MAXN; i++){
//   res1[i]=res[i]=b[i]-p[i];
//   p[i]=0.0;
//   }

//   double norm_b=0.0;
//   for(int i=0; i<MAXN; i++)  norm_b+=b[i]*b[i];
//   norm_b=sqrt(norm_b);

//   double norm_A=FrobeniusNorm(A);

//   bool NextIteration=true;
//   int Iteration=1;

//   while(NextIteration){
//     valarray<double> z(0.0, MAXN);
//     valarray<double> z1(0.0, MAXN);

//     if(pm==0) {z=res/A.CDIAG; z1=res1/A.CDIAG;}
//     else {CholFactSolve(*CH, res, z); CholFactSolve(CH, res1, z1);}

//     ro0=ro1;
//     ro1=0;
//     for(int i=0; i<MAXN; i++)  ro1=ro1+res1[i]*z[i];
//     if(ro1==0) throw Error("BiCG failed\n");
//     if(Iteration == 1) {p=z; p1=z1;}
//     else{
//       beta=ro1/ro0;
//       p=z+beta*p;
//       p1=z1+beta*p1;
//     }
//     valarray<double> q(0.0, MAXN);
//     valarray<double> q1(0.0, MAXN);
//     SymmMatrixVectorProduct(A,p,q);
//     SymmMatrixVectorProduct(A,p1,q1);
//     double denom=0;
//     for(int i=0; i<MAXN; i++) denom+=p1[i]*q[i];
//     double aa=ro1/denom;
//     x+=aa*p;
//     res-=aa*q;
//     res1-=aa*q1;
//     double norm_r=0;
//     double norm_x=0;
//     for(int i=0; i<MAXN; i++){
//       norm_x+=x[i]*x[i];
//       norm_r+=res[i]*res[i];
//     }
//      norm_x=sqrt(norm_x);
//     norm_r=sqrt(norm_r);
//     info_stream(2)<<"norm_x= "<<norm_x<<'\n';
//     info_stream(2)<<"norm_r= "<<norm_r<<'\n';
//     info_stream(2)<<"tolerance= "<<tolerance<<'\n';

//     double rlimit=tolerance*(norm_A*norm_x+norm_b);
//     info_stream(2)<<" BiCG: Iteration = "<< Iteration <<" ||res|| = "<<norm_r;
//     info_stream(2)<<" rlimit = "<<rlimit<<'\n';
//     Iteration++;

//     if(norm_r<rlimit){ NextIteration=false; retcode=0;}
//     if(Iteration==max_iter){
//       info_stream(2) << "BiCG: Max iteration number reached. Exiting...\n";
//       NextIteration=false;
//       retcode=1;
//     }
//   }
//  info_stream(2)<<"BiCG: completed \n";
//   return retcode;
// };
