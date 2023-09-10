/*!\file RandTemp.h
   \brief Class for generating random numbers

   \author os - IBT, Universitaet Karlsruhe (TH)

   os May 9, 2001

   NOTICE:
   this class includes methods for generating
   uniformly distributed pseudo-random values (numerical recipes, p.280)
   and
   normaly distributed pseudo-random values (ACM Algorithm 712, J.L.Leva)

   The <double> version of the template should be tested.
   Although I've already used it with double and had no complains. -os


 */


#ifndef _RandTemp_h
#define _RandTemp_h

#include <cmath>

#include <kaMachineOS.h>
#include <ctime>

/*!\class RandTemp -
   \brief template class for random number generation

   Description:
   The class includes methods for generating
   distributed pseudo-random values (numerical recipes, p.280)

   distributed pseudo-random values ( numerical recipes p. 289)
 */
template<class T> class RandTemp {
 public:
  RandTemp();
  ~RandTemp() {delete[] iv;}

  /*!Call init with a negative integer to (re-)initialize the generator
     DO NOT call init between successive deviates in a sequence
     Normally you don't need to call init, since the generator
     is initialized automatically with system time
     during object construction
   */
  inline void init(int a) {idum = a;}

  //! Uniformly distributed random number generator
  inline T randu(void);

  //! Normally distributed random number generator. Standard deviation 1.0
  inline T randn(void);

 private:
  RandTemp(const RandTemp &) {}

  void operator=(const RandTemp &) {}

  int idum;  // seed
  int iy;

  static const int IA;
  static const int IM;
  static const T AM;
  static const int IQ;
  static const int IR;
  static const int NTAB;
  static const int NDIV;
  static T EPS;
  static T RNMX;

  int *iv;
};  // class RandTemp

// ----static constant data member definitions
// WARNING: double version on this template may (or must:)) not
// work properly


template<class T>
const int RandTemp<T>::IA = 16807;

template<class T>
const int RandTemp<T>::IM = 2147483647;

template<class T>
const T RandTemp<T>::AM = (1.0/RandTemp<T>::IM);

template<class T>
const int RandTemp<T>::IQ =  127773;

template<class T>
const int RandTemp<T>::IR = 2836;

template<class T>
const int RandTemp<T>::NTAB = 32;

template<class T>
const int RandTemp<T>::NDIV = (1+(RandTemp<T>::IM-1)/RandTemp<T>::NTAB);

template<class T>
T RandTemp<T>::EPS;

template<class T>
T RandTemp<T>::RNMX;


template<class T> RandTemp<T>::RandTemp() {
  iy = 0;
  time_t a;
  init(-time(&a));  // seed the generator with the system time

  assert( ((int)a) >= 0);  // check if  system time >= 0

  // get machine precision
  EPS = 1.2e-7;
  while ((T)(1.0 + EPS) > 1.0) EPS *= 0.5;

  RNMX = 1.0 - EPS;

  iv = new int[NTAB];
}

template<class T>
inline T RandTemp<T>::randu() {
  int j;
  int k;
  T   temp;

  if ((idum <= 0) || !iy) {
    if (-(idum) < 1)
      idum = 1;
    else
      idum = -(idum);
    for (j = NTAB+7; j >= 0; j--) {
      k    = (idum)/IQ;
      idum = IA*(idum-k*IQ)-IR*k;
      if (idum < 0)
        idum += IM;
      if (j < NTAB)
        iv[j] = idum;
    }
    iy = iv[0];
  }
  k    = (idum)/IQ;
  idum = IA*(idum-k*IQ)-IR*k;
  if (idum < 0)
    idum += IM;
  j     = iy/NDIV;
  iy    = iv[j];
  iv[j] = idum;
  if ((temp = AM*iy) > RNMX)
    return RNMX;
  else
    return temp;
}  // >::randu

template<class T>
inline T RandTemp<T>::randn() {
  static int iset = 0;
  static T   gset;
  T fac, rsq, v1, v2;

  if  (iset == 0) {
    do {
      v1  = 2.0*randu()-1.0;
      v2  = 2.0*randu()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac  = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

#endif  // ifndef _RandTemp_h
