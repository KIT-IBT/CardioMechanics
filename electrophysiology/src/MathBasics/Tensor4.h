/*
**      Name
**              Tensor4.h
**
**      Usages
**              eg FEMechanicScalar.h
**
**      Remarks
**              manage asymmetric and symmetric tensors of 4nd order
**
**      History
**              9.2.99 -fs creation
**              7.4.01 -fs asymmetric tensors, deformation gradient tensors
**              8.8.01 -fs Tensor4
**
**  created at IBT - Universität Karlsruhe
*/

#ifndef __TENSOR4_H
#define __TENSOR4_H

#include <Tensor2.h>

template<class T> class SymmetricTensor4;

/*!template class Tensor4
 * tensor of 4. order represented by a 3 x 3 x 3 x 3 matrix
 * compatible to SymmetricTensor4 */

template<class T>
class Tensor4 {
 public:
  Tensor2<T> a[3][3];

  inline Tensor4() {}

  inline Tensor4(const kaMatrixN<T, 6> &rm) {
    *this = SymmetricTensor4<T>(rm);
  }

  inline Tensor4(const SymmetricTensor4<T> &st) {
    *this = st;
  }

  inline void Null() {
    int i, j;

    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        a[i][j].Null();
  }

  inline Tensor2<T> *operator[](unsigned int i) {return a[i];}

  inline const Tensor2<T> *operator[](unsigned int i) const {return a[i];}

  inline const Tensor2<T> & Get(unsigned int i, unsigned int j) const {return a[i][j];}

  inline Tensor4 & operator=(const SymmetricTensor4<T> &st) {
    a[0][0] = st.a;
    a[0][1] = st.b;
    a[0][2] = st.c;

    a[1][0] = st.b;
    a[1][1] = st.d;
    a[1][2] = st.e;

    a[2][0] = st.c;
    a[2][1] = st.e;
    a[2][2] = st.f;

    return *this;
  }

  inline const Tensor4 operator-(const Tensor4 &S) const {
    Tensor4 tmp;
    int i, j;

    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        tmp[i][j] = a[i][j]-S[i][j];

    return tmp;
  }

  inline void CoordinateTransform2(const kaMatrixN<T, 3> R) {
    Tensor4<T> tmp;
    int i, j, k, l, m, n;
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++)
          for (l = 0; l < 3; l++) {
            const T *Rl = R[l];
            T t0        = 0;
            for (m = 0; m < 3; m++) {
              T t1 = 0;
              for (n = 0; n < 3; n++) {
                const Tensor2<T> &amn = a[m][n];
                T t2                  = R[k][0]*(Rl[0]*amn[0][0]+Rl[1]*amn[0][1]+Rl[2]*amn[0][2]);
                t2 += R[k][1]*(Rl[0]*amn[1][0]+Rl[1]*amn[1][1]+Rl[2]*amn[1][2]);
                t2 += R[k][2]*(Rl[0]*amn[2][0]+Rl[1]*amn[2][1]+Rl[2]*amn[2][2]);
                t1 += R[j][n]*t2;
              }
              t0 += R[i][m]*t1;
            }
            tmp[i][j][k][l] = t0;
          }
    *this = tmp;
  }

  inline void CoordinateTransform(const kaMatrixN<T, 3> R) {
    int i, j;
    T   tmp[3][3][9];

    const T *R0 = R[0], *R1 = R[1], *R2 = R[2];
    const T  R00 = R0[0], R01 = R0[1], R02 = R0[2];
    const T  R10 = R1[0], R11 = R1[1], R12 = R1[2];
    const T  R20 = R2[0], R21 = R2[1], R22 = R2[2];
    T *p = tmp[0][0];

    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) {
        const T *aijk = a[i][j].Pointer();

        *p++  = R00*aijk[0]+R01*aijk[1]+R02*aijk[2];
        *p++  = R10*aijk[0]+R11*aijk[1]+R12*aijk[2];
        *p++  = R20*aijk[0]+R21*aijk[1]+R22*aijk[2];
        aijk += 3;
        *p++  = R00*aijk[0]+R01*aijk[1]+R02*aijk[2];
        *p++  = R10*aijk[0]+R11*aijk[1]+R12*aijk[2];
        *p++  = R20*aijk[0]+R21*aijk[1]+R22*aijk[2];
        aijk += 3;
        *p++  = R00*aijk[0]+R01*aijk[1]+R02*aijk[2];
        *p++  = R10*aijk[0]+R11*aijk[1]+R12*aijk[2];
        *p++  = R20*aijk[0]+R21*aijk[1]+R22*aijk[2];
      }
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) {
        p = a[i][j].Pointer();
        T *pTmp = tmp[i][j];

        *p++ = R00*pTmp[0]+R01*pTmp[3]+R02*pTmp[6];
        *p++ = R00*pTmp[1]+R01*pTmp[4]+R02*pTmp[7];
        *p++ = R00*pTmp[2]+R01*pTmp[5]+R02*pTmp[8];

        *p++ = R10*pTmp[0]+R11*pTmp[3]+R12*pTmp[6];
        *p++ = R10*pTmp[1]+R11*pTmp[4]+R12*pTmp[7];
        *p++ = R10*pTmp[2]+R11*pTmp[5]+R12*pTmp[8];

        *p++ = R20*pTmp[0]+R21*pTmp[3]+R22*pTmp[6];
        *p++ = R20*pTmp[1]+R21*pTmp[4]+R22*pTmp[7];
        *p   = R20*pTmp[2]+R21*pTmp[5]+R22*pTmp[8];
      }
    p = tmp[0][0];
    for (i = 0; i < 3; i++) {
      const T *ai0c = a[i][0].Pointer(), *ai1c = a[i][1].Pointer(), *ai2c = a[i][2].Pointer();
      for (j = 0; j < 3; j++) {
        const T *Rj = R[j];
        const T  Rj0 = Rj[0], Rj1 = Rj[1], Rj2 = Rj[2];
        const T *ai0 = ai0c, *ai1 = ai1c, *ai2 = ai2c;

        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;

        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;

        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
        ai0++; ai1++; ai2++;
        *p++ = Rj0 **ai0+Rj1 **ai1+Rj2 **ai2;
      }
    }
    for (i = 0; i < 3; i++) {
      const T *Ri = R[i];
      const T  Ri0 = Ri[0], Ri1 = Ri[1], Ri2 = Ri[2];
      for (j = 0; j < 3; j++) {
        T *p = a[i][j].Pointer();
        T *pTmp0 = tmp[0][j], *pTmp1 = tmp[1][j], *pTmp2 = tmp[2][j];

        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;

        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;

        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
        pTmp0++; pTmp1++; pTmp2++;
        *p++ = Ri0 **pTmp0+Ri1 **pTmp1+Ri2 **pTmp2;
      }
    }
  }  // CoordinateTransform
};  // class Tensor4


template<class T>
inline ostream & operator<<(ostream &os, const Tensor4<T> &t) {
  int i, j;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      os << t[i][j] << '\n';

  return os;
}

/*!template class SymmetricTensor4
 * tensor of 4. order represented by a 6 x 6 dimensional vector
 * compatible to Tensor4 */

template<class T>
class SymmetricTensor4 : public SymmetricTensor2<SymmetricTensor2<T>> {
 public:
  inline SymmetricTensor4() {}

  inline SymmetricTensor4(const Tensor4<T> &rm) {
    SymmetricTensor2<SymmetricTensor2<T>>::a = rm.Get(0);
    SymmetricTensor2<SymmetricTensor2<T>>::b = rm.Get(1);
    SymmetricTensor2<SymmetricTensor2<T>>::c = rm.Get(2);
    SymmetricTensor2<SymmetricTensor2<T>>::d = rm.Get(4);
    SymmetricTensor2<SymmetricTensor2<T>>::e = rm.Get(5);
    SymmetricTensor2<SymmetricTensor2<T>>::f = rm.Get(8);
  }

  inline SymmetricTensor4(const kaMatrixN<T, 6> &rm) {
    const T *pm = rm.Pointer();

    SymmetricTensor2<SymmetricTensor2<T>>::a.a = pm[0];
    SymmetricTensor2<SymmetricTensor2<T>>::a.b = pm[3];
    SymmetricTensor2<SymmetricTensor2<T>>::a.c = pm[4];
    SymmetricTensor2<SymmetricTensor2<T>>::a.d = pm[1];
    SymmetricTensor2<SymmetricTensor2<T>>::a.e = pm[5];
    SymmetricTensor2<SymmetricTensor2<T>>::a.f = pm[2];

    SymmetricTensor2<SymmetricTensor2<T>>::d.a = pm[6];
    SymmetricTensor2<SymmetricTensor2<T>>::d.b = pm[9];
    SymmetricTensor2<SymmetricTensor2<T>>::d.c = pm[10];
    SymmetricTensor2<SymmetricTensor2<T>>::d.d = pm[7];
    SymmetricTensor2<SymmetricTensor2<T>>::d.e = pm[11];
    SymmetricTensor2<SymmetricTensor2<T>>::d.f = pm[8];

    SymmetricTensor2<SymmetricTensor2<T>>::f.a = pm[12];
    SymmetricTensor2<SymmetricTensor2<T>>::f.b = pm[15];
    SymmetricTensor2<SymmetricTensor2<T>>::f.c = pm[16];
    SymmetricTensor2<SymmetricTensor2<T>>::f.d = pm[13];
    SymmetricTensor2<SymmetricTensor2<T>>::f.e = pm[17];
    SymmetricTensor2<SymmetricTensor2<T>>::f.f = pm[14];

    SymmetricTensor2<SymmetricTensor2<T>>::b.a = pm[18];
    SymmetricTensor2<SymmetricTensor2<T>>::b.b = pm[21];
    SymmetricTensor2<SymmetricTensor2<T>>::b.c = pm[22];
    SymmetricTensor2<SymmetricTensor2<T>>::b.d = pm[19];
    SymmetricTensor2<SymmetricTensor2<T>>::b.e = pm[23];
    SymmetricTensor2<SymmetricTensor2<T>>::b.f = pm[20];

    SymmetricTensor2<SymmetricTensor2<T>>::c.a = pm[24];
    SymmetricTensor2<SymmetricTensor2<T>>::c.b = pm[27];
    SymmetricTensor2<SymmetricTensor2<T>>::c.c = pm[28];
    SymmetricTensor2<SymmetricTensor2<T>>::c.d = pm[25];
    SymmetricTensor2<SymmetricTensor2<T>>::c.e = pm[29];
    SymmetricTensor2<SymmetricTensor2<T>>::c.f = pm[26];

    SymmetricTensor2<SymmetricTensor2<T>>::e.a = pm[30];
    SymmetricTensor2<SymmetricTensor2<T>>::e.b = pm[33];
    SymmetricTensor2<SymmetricTensor2<T>>::e.c = pm[34];
    SymmetricTensor2<SymmetricTensor2<T>>::e.d = pm[31];
    SymmetricTensor2<SymmetricTensor2<T>>::e.e = pm[35];
    SymmetricTensor2<SymmetricTensor2<T>>::e.f = pm[32];
  }
};  // class SymmetricTensor4


template<class T>
inline ostream & operator<<(ostream &os, const SymmetricTensor4<T> &t) {
  os << t.a << t.b << t.c << t.d << t.e << t.f;
  return os;
}

#endif  // ifndef __TENSOR4_H
