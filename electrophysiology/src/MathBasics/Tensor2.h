/*
**      Name
**              Tensor2.h
**
**      Usages
**              eg FDYZ.cpp
**
**      Remarks
**              manage asymmetric and symmetric tensors of 2nd order
**              compatible to class kaMatrixN
**
**      History
**              9.2.99 -fs creation
**              7.4.01 -fs asymmetric tensors, deformation gradient tensors
**
**  created at IBT - Universität Karlsruhe
*/

#ifndef TENSOR2_H
#define TENSOR2_H

#include <kaMatrixN.h>
#include <kaTabularizedFunction.h>
#include <Jacobi.h>
#include <kaPoint.h>

template<class T> class SymmetricTensor2;

/*!template class Tensor2
 * tensor of 2. order represented by a 3 x 3 matrix
 * compatible to kaMatrixN (base class) and SymmetricTensor2 */

template<class T>
class Tensor2 : public kaMatrixN<T, 3>  {
 public:
  /*!default constructor */
  inline Tensor2() {}

  /*!constructor for transfer with kaMatrixN */
  inline Tensor2(const kaMatrixN<T, 3> &m) : kaMatrixN<T, 3>(m) {}

  /*!constructor for transfer with SymmetricTensor2 */
  inline Tensor2(const SymmetricTensor2<T> &st2) {
    T *ps = this->pElements;

    *ps++ = st2.a;
    *ps++ = st2.b;
    *ps++ = st2.c;

    *ps++ = st2.b;
    *ps++ = st2.d;
    *ps++ = st2.e;

    *ps++ = st2.c;
    *ps++ = st2.e;
    *ps   = st2.f;
  }

  /*!assignment operator for transfer with kaMatrixN */
  inline Tensor2 & operator=(const kaMatrixN<T, 3> &m) {
    T *ps        = this->pElements;
    const T *pss = m.Pointer();

    *ps++ = *pss++;
    *ps++ = *pss++;
    *ps++ = *pss++;

    *ps++ = *pss++;
    *ps++ = *pss++;
    *ps++ = *pss++;

    *ps++ = *pss++;
    *ps++ = *pss++;
    *ps   = *pss;

    return *this;
  }

  /*!assignment operator for transfer with SymmetricTensor2 */
  inline Tensor2 & operator=(const SymmetricTensor2<T> &st2) {
    T *ps = this->pElements;

    *ps++ = st2.a;
    *ps++ = st2.b;
    *ps++ = st2.c;

    *ps++ = st2.b;
    *ps++ = st2.d;
    *ps++ = st2.e;

    *ps++ = st2.c;
    *ps++ = st2.e;
    *ps   = st2.f;

    return *this;
  }

  /*!
   * Calculate the deformation gradient tensor with 6 neighbourhood
   * using the current offset lattice position, the index positions of the neighbours, and the distances
   * FD method in star
   */

  template<class U, class IndexType>

  inline void DeformationGradientTensor(U *pox, U *poy, U *poz, IndexType ixm, IndexType iym, IndexType izm,
                                        IndexType ixp, IndexType iyp, IndexType izp, T ddx, T ddy, T ddz) {
    this->a(0) = (*(pox+ixp)-*(pox+ixm))*ddx+1.0;
    this->a(3) = (*(poy+ixp)-*(poy+ixm))*ddx;
    this->a(6) = (*(poz+ixp)-*(poz+ixm))*ddx;

    this->a(1) = (*(pox+iyp)-*(pox+iym))*ddy;
    this->a(4) = (*(poy+iyp)-*(poy+iym))*ddy+1.0;
    this->a(7) = (*(poz+iyp)-*(poz+iym))*ddy;

    this->a(2) = (*(pox+izp)-*(pox+izm))*ddz;
    this->a(5) = (*(poy+izp)-*(poy+izm))*ddz;
    this->a(8) = (*(poz+izp)-*(poz+izm))*ddz+1.0;
  }

  /*!
   * Calculate the deformation gradient tensor in center of cell
   * using the current offset lattice position and the index positions of the neighbours
   * FE method with linear form functions
   */

  template<class U, class IndexType>

  inline void DeformationGradientTensor(U *pox, U *poy, U *poz, IndexType ixp, IndexType iyp, IndexType izp) {
    this->a(0) =
      (*(pox+ixp)+*(pox+ixp+iyp)+*(pox+ixp+izp)+*(pox+ixp+iyp+izp)-*(pox)-*(pox+iyp)-*(pox+izp)-*(pox+iyp+izp))*.25+1.0;
    this->a(3) =
      (*(poy+ixp)+*(poy+ixp+iyp)+*(poy+ixp+izp)+*(poy+ixp+iyp+izp)-*(poy)-*(poy+iyp)-*(poy+izp)-*(poy+iyp+izp))*.25;
    this->a(6) =
      (*(poz+ixp)+*(poz+ixp+iyp)+*(poz+ixp+izp)+*(poz+ixp+iyp+izp)-*(poz)-*(poz+iyp)-*(poz+izp)-*(poz+iyp+izp))*.25;

    this->a(1) =
      (*(pox+iyp)+*(pox+ixp+iyp)+*(pox+iyp+izp)+*(pox+ixp+iyp+izp)-*(pox)-*(pox+ixp)-*(pox+izp)-*(pox+ixp+izp))*.25;
    this->a(4) =
      (*(poy+iyp)+*(poy+ixp+iyp)+*(poy+iyp+izp)+*(poy+ixp+iyp+izp)-*(poy)-*(poy+ixp)-*(poy+izp)-*(poy+ixp+izp))*.25+1.0;
    this->a(7) =
      (*(poz+iyp)+*(poz+ixp+iyp)+*(poz+iyp+izp)+*(poz+ixp+iyp+izp)-*(poz)-*(poz+ixp)-*(poz+izp)-*(poz+ixp+izp))*.25;

    this->a(2) =
      (*(pox+izp)+*(pox+ixp+izp)+*(pox+iyp+izp)+*(pox+ixp+iyp+izp)-*(pox)-*(pox+ixp)-*(pox+iyp)-*(pox+ixp+iyp))*.25;
    this->a(5) =
      (*(poy+izp)+*(poy+ixp+izp)+*(poy+iyp+izp)+*(poy+ixp+iyp+izp)-*(poy)-*(poy+ixp)-*(poy+iyp)-*(poy+ixp+iyp))*.25;
    this->a(8) =
      (*(poz+izp)+*(poz+ixp+izp)+*(poz+iyp+izp)+*(poz+ixp+iyp+izp)-*(poz)-*(poz+ixp)-*(poz+iyp)-*(poz+ixp+iyp))*.25+1.0;
  }

  inline const Tensor2 operator-(const Tensor2 &S) const {
    Tensor2 tmp;
    int i, j;

    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        tmp[i][j] = this->Get(i, j)-S[i][j];

    return tmp;
  }

  inline const Tensor2 operator*(Tensor2 p) const {
    Tensor2 temp;
    T *pTemp = temp.pElements;

    const T *pl = this->pElements;

    for (unsigned int  i = 0; i < 3; i++, pl += 3) {
      const T *pr = p.pElements;
      for (unsigned int j = 0; j < 3; j++, pr++) {
        T val = 0;
        const T *plk = pl, *plke = pl+3, *prk = pr;
        for (; plk < plke; plk++, prk += 3)
          val += *plk **prk;
        *pTemp++ = val;
      }
    }
    return temp;
  }

  template<class U> const kaPoint<U> operator*(kaPoint<U> p) const {
    const U x = this->Get(0)*p.x+this->Get(1)*p.y+this->Get(2)*p.z;
    const U y = this->Get(3)*p.x+this->Get(4)*p.y+this->Get(5)*p.z;
    const U z = this->Get(6)*p.x+this->Get(7)*p.y+this->Get(8)*p.z;

    return kaPoint<U>(x, y, z);
  }

  /*!
   * Polar Decomposition
   * A 2. order tensor is decomposed into a rotation R and symmetric scaling tensor U
   */

  inline void PolarDecomposition(Tensor2<T> &R, SymmetricTensor2<T> &U) {
    Tensor2 C(*this);  // C <- F^T F

    C.Transpose();
    C.MatrixNX(*this);

    T d[3];  // C -> S D^2 S
    T a[3][3], v[3][3];
    C.PutIntoArray((T *)&a[0][0]);
    bool rc = jacobi<T, 3>(a, d, v);
    assert(rc);
    Tensor2<T> S;
    S.SetFromArray((T *)&v[0][0]);

    assert(!S.testnan());

    // cerr << "S: " << S << "\n";

    U.a = sqrt(d[0]);  // U <- S D S
    U.b = 0.;
    U.c = 0.;
    U.d = sqrt(d[1]);
    U.e = 0.;
    U.f = sqrt(d[2]);
    U.CoordinateTransform(S);

    // cerr << "U: " << U << "\n";

    if (U.testnan()) {
      cerr << "d " << d[0] << '\t' << d[1] << '\t' << d[2] << '\t' << "\n";
      cerr << "S: " << S << "\n";
      cerr << "U: " << U << "\n";
    }

    SymmetricTensor2<T> UInv(U);  // R -< F U^-1
    rc = UInv.Inverse();
    assert(rc);
    Tensor2<T> UI(UInv);

    assert(!UInv.testnan());

    R = *this;
    R.MatrixNX(UI);

    // cerr << "R: " << R << "\n";
  }  // PolarDecomposition

  inline void CoordinateTransform(const kaMatrixN<T, 3> R) {
    MatrixNX(R);
    Tensor2<T> RT = R;
    RT.Transpose();
    XMatrixN(RT);
  }

  inline void InverseCoordinateTransform(const kaMatrixN<T, 3> R) {
    Tensor2<T> RT = R;
    RT.Transpose();
    MatrixNX(RT);
    XMatrixN(R);
  }

  inline bool testnan() {
    return ::isnan(this->a(0)) || ::isnan(this->a(1)) || ::isnan(this->a(2)) || ::isnan(this->a(3)) || ::isnan(this->a(
                                                                                                                 4)) ||
           ::isnan(this->a(5)) || ::isnan(this->a(6)) || ::isnan(this->a(7)) || ::isnan(this->a(8));
  }
};  // class Tensor2


template<class T>
inline ostream & operator<<(ostream &os, const Tensor2<T> &t) {
  os << t.Get(0, 0) << '\t' << t.Get(0, 1) << '\t' << t.Get(0, 2) << '\n';
  os << t.Get(1, 0) << '\t' << t.Get(1, 1) << '\t' << t.Get(1, 2) << '\n';
  os << t.Get(2, 0) << '\t' << t.Get(2, 1) << '\t' << t.Get(2, 2) << '\n';
  return os;
}

/*!template class SymmetricTensor2
 * tensor of 2. order represented by a 6 dimensional vector
 * compatible to kaMatrixN and Tensor2 */

template<class T>
class SymmetricTensor2 {
 public:
  T a, b, c, d, e, f;

  inline SymmetricTensor2() {}

  inline SymmetricTensor2(T aArg, T bArg, T cArg, T dArg, T eArg, T fArg)
  {a = aArg; b = bArg; c = cArg; d = dArg; e = eArg; f = fArg;}

  inline SymmetricTensor2(const kaMatrixN<T, 3> &m) {
    const T *pss = m.Pointer();

    a = pss[0];
    b = pss[1];
    c = pss[2];
    d = pss[4];
    e = pss[5];
    f = pss[8];
  }

  inline SymmetricTensor2(const SymmetricTensor2<T> &rs) {*this = rs;}

  inline SymmetricTensor2 & operator=(const kaMatrixN<T, 3> &m) {
    const T *pss = m.Pointer();

    a = pss[0];
    b = pss[1];
    c = pss[2];
    d = pss[4];
    e = pss[5];
    f = pss[8];

    return *this;
  }

  inline SymmetricTensor2 & operator=(const SymmetricTensor2<T> &rs) {
    a = rs.a;
    b = rs.b;
    c = rs.c;
    d = rs.d;
    e = rs.e;
    f = rs.f;

    return *this;
  }

  inline SymmetricTensor2 & operator+=(const SymmetricTensor2 &S) {
    a += S.a;
    b += S.b;
    c += S.c;
    d += S.d;
    e += S.e;
    f += S.f;

    return *this;
  }

  inline SymmetricTensor2 & operator+=(const kaMatrixN<T, 3> &rm) {
    const T *pss = rm.Pointer();

    a += pss[0];
    b += pss[1];
    c += pss[2];
    d += pss[4];
    e += pss[5];
    f += pss[8];

    return *this;
  }

  inline SymmetricTensor2 & operator/=(T dArg) {
    assert(dArg != 0.);
    dArg = 1./dArg;

    a *= dArg;
    b *= dArg;
    c *= dArg;
    d *= dArg;
    e *= dArg;
    f *= dArg;

    return *this;
  }

  inline SymmetricTensor2 & operator*=(T dArg) {
    a *= dArg;
    b *= dArg;
    c *= dArg;
    d *= dArg;
    e *= dArg;
    f *= dArg;

    return *this;
  }

  inline bool operator==(const SymmetricTensor2 &S) const
  {return a == S.a && b == S.b && c == S.c && d == S.d && e == S.e && f == S.f;}

  inline const SymmetricTensor2 operator*(T fArg) const {
    return SymmetricTensor2(a*fArg, b*fArg, c*fArg, d*fArg, e*fArg, f*fArg);
  }

  inline const SymmetricTensor2 operator*(const SymmetricTensor2 &S) const {
    return SymmetricTensor2(a*S.a+b*S.b+c*S.c, a*S.b+b*S.d+c*S.e, a*S.c+b*S.e+c*S.f, b*S.b+d*S.d+e*S.e,
                            b*S.c+d*S.e+e*S.f, c*S.c+e*S.e+f*S.f);
  }

  inline const SymmetricTensor2 operator+(const SymmetricTensor2 &S) const {
    return SymmetricTensor2(a+S.a, b+S.b, c+S.c, d+S.d, e+S.e, f+S.f);
  }

  inline const SymmetricTensor2 operator-(const SymmetricTensor2 &S) const {
    return SymmetricTensor2(a-S.a, b-S.b, c-S.c, d-S.d, e-S.e, f-S.f);
  }

  /*! CoordinateTransform
   * basis transformation of tensor 2. order: tensor' = R tensor R^T
   */

  inline void CoordinateTransform(const kaMatrixN<T, 3> R) {
    const T at = a, bt = b, ct = c, dt = d, et = e, ft = f;
    const T *R0 = R[0], *R1 = R[1], *R2 = R[2];
    const T  R00 = R0[0], R01 = R0[1], R02 = R0[2];
    const T  R10 = R1[0], R11 = R1[1], R12 = R1[2];
    const T  R20 = R2[0], R21 = R2[1], R22 = R2[2];

    a = at*R00*R00+2.*bt*R00*R01+dt*R01*R01+2.*ct*R00*R02+2.*et*R01*R02+ft*R02*R02;
    b = at*R00*R10+bt*R01*R10+ct*R02*R10+bt*R00*R11+dt*R01*R11+et*R02*R11+ct*R00*R12+et*R01*R12+ft*R02*R12;
    c = at*R00*R20+bt*R01*R20+ct*R02*R20+bt*R00*R21+dt*R01*R21+et*R02*R21+ct*R00*R22+et*R01*R22+ft*R02*R22;
    d = at*R10*R10+2.*bt*R10*R11+dt*R11*R11+2.*ct*R10*R12+2.*et*R11*R12+f*R12*R12;
    e = at*R10*R20+bt*R11*R20+ct*R12*R20+bt*R10*R21+dt*R11*R21+et*R12*R21+ct*R10*R22+et*R11*R22+ft*R12*R22;
    f = at*R20*R20+2.*bt*R20*R21+dt*R21*R21+2.*ct*R20*R22+2.*et*R21*R22+ft*R22*R22;
  }

  /*! InverseCoordinateTransform
   * basis transformation of tensor 2. order: tensor' = R^T tensor R
   */

  inline void InverseCoordinateTransform(const kaMatrixN<T, 3> R) {
    const T at = a, bt = b, ct = c, dt = d, et = e, ft = f;
    const T *R0 = R[0], *R1 = R[1], *R2 = R[2];
    const T  R00 = R0[0], R01 = R1[0], R02 = R2[0];
    const T  R10 = R0[1], R11 = R1[1], R12 = R2[1];
    const T  R20 = R0[2], R21 = R1[2], R22 = R2[2];

    a = at*R00*R00+2.*bt*R00*R01+dt*R01*R01+2.*ct*R00*R02+2.*et*R01*R02+ft*R02*R02;
    b = at*R00*R10+bt*R01*R10+ct*R02*R10+bt*R00*R11+dt*R01*R11+et*R02*R11+ct*R00*R12+et*R01*R12+ft*R02*R12;
    c = at*R00*R20+bt*R01*R20+ct*R02*R20+bt*R00*R21+dt*R01*R21+et*R02*R21+ct*R00*R22+et*R01*R22+ft*R02*R22;
    d = at*R10*R10+2.*bt*R10*R11+dt*R11*R11+2.*ct*R10*R12+2.*et*R11*R12+f*R12*R12;
    e = at*R10*R20+bt*R11*R20+ct*R12*R20+bt*R10*R21+dt*R11*R21+et*R12*R21+ct*R10*R22+et*R11*R22+ft*R12*R22;
    f = at*R20*R20+2.*bt*R20*R21+dt*R21*R21+2.*ct*R20*R22+2.*et*R21*R22+ft*R22*R22;
  }

  /*!SetRDRT
   * Construct 2. order tensor with diagonal matrix (KAnisotropyX, ..., KAnisotropyZ)
   * and rotational matrix (cphi, cheta)
   */

  inline void SetRDRT(unsigned char cphi, unsigned char ctheta, T KAnisotropyX, T KAnisotropyY, T KAnisotropyZ) {
    a = KAnisotropyX;
    d = KAnisotropyY;
    f = KAnisotropyZ;
    b = c = e = 0.;

    assert(cphi < 255);

    kaMatrixN<T, 3> R;

    // Seems that signs are swapped in the third column of R.
    // Does not make a difference when applied to diagonal matrices as done in this function, though.
    R.a(0) = SinTable.val[ctheta]*CosTable.val[cphi]; R.a(1) = -SinTable.val[cphi];
    R.a(2) = CosTable.val[ctheta]*CosTable.val[cphi];
    R.a(3) = SinTable.val[ctheta]*SinTable.val[cphi]; R.a(4) = CosTable.val[cphi];
    R.a(5) = CosTable.val[ctheta]*SinTable.val[cphi];
    R.a(6) = CosTable.val[ctheta];                    R.a(7) = 0.;                  R.a(8) = -SinTable.val[ctheta];

    CoordinateTransform(R);
  }

  inline void Null() {a = b = c = d = e = f = 0;}

  inline void Identity() {a = d = f = 1; b = c = e = 0;}

  inline void Transpose() {}

  inline bool Inverse(void) {
    T det = Determinant();

    if (!det)
      return false;

    det = 1./det;

    SymmetricTensor2<T> S;
    S.a   = -e*e + d*f;
    S.b   =  c*e - b*f;
    S.c   = -c*d + b*e;
    S.d   = -c*c + a*f;
    S.e   =  c*b - a*e;
    S.f   = -b*b + a*d;
    *this = S*det;

    return true;
  }

  inline T Determinant(void) const {return a*d*f + 2.0*b*e*c - d*c*c - e*e*a - f*b*b;}

  inline bool testnan() {
    return ::isnan(a) || ::isnan(b) || ::isnan(c) || ::isnan(d) || ::isnan(e) || ::isnan(f);
  }
};  // class SymmetricTensor2


template<class T>
inline const SymmetricTensor2<T> operator*(T fArg, const SymmetricTensor2<T> &tArg) {
  return SymmetricTensor2<T>(tArg.a*fArg, tArg.b*fArg, tArg.c*fArg, tArg.d*fArg, tArg.e*fArg, tArg.f*fArg);
}

template<class T>
inline ostream & operator<<(ostream &os, const SymmetricTensor2<T> &t) {
  os << t.a << '\t' << t.b << '\t' << t.c << '\t' << t.d << '\t' << t.e << '\t' << t.f << '\n';
  return os;
}

#endif  // ifndef TENSOR2_H
