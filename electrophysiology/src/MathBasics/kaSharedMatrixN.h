/*! \file kaSharedMatrixN.h
   \brief Class for handling of quadratic matrices in shared memory

   \author cw,fs,mm,os, IBT - Universität Karlsruhe (TH)
 */

#ifndef KASHAREDMATRIXN_H
#define KASHAREDMATRIXN_H

#include <kaMatrixMathN.h>

#define NOT_USED(x) (void)(x)

//! Class for handling of quadratic matrices in shared memory.
/*!
   kaSharedMatrixN includes variant methods for matrix mathematics and file-io.

   \author fs,os, IBT - Universität Karlsruhe (TH)
 */

template<class T> class kaSharedMatrixN : public kaMatrixMathN<T> {
 protected:
  // df

  unsigned int *pn;  //!< Pointer to number of rows and columns
  T *pElements;     //!< Pointer to elements

 public:
  //! Default constructor
  kaSharedMatrixN(unsigned int dimension = 0);

  //! Constructor with kaSharedMatrixN
  kaSharedMatrixN(const kaSharedMatrixN<T> &);

  //! Destructor
  virtual ~kaSharedMatrixN();

  //! Pointer to elements
  inline virtual T *Pointer() {return pElements;}

  //! Const pointer to elements
  inline virtual const T *Pointer() const {return pElements;}

  //! Dimension of matrix
  virtual inline unsigned int Dimension() const {return *pn;}

  //! Number of rows (=Dimension)
  inline unsigned int NumRows() const {return *pn;}

  //! Number of columns (=Dimension)
  inline unsigned int NumCols() const {return *pn;}

  kaSharedMatrixN & operator=(const kaSharedMatrixN<T> &);

  const kaSharedMatrixN operator*(const T);
  const kaSharedMatrixN operator/(const T);
  kaSharedMatrixN     & operator*=(const T);
  kaSharedMatrixN     & operator/=(const T);

  kaSharedMatrixN   operator+(const kaSharedMatrixN &) const;
  kaSharedMatrixN   operator-(const kaSharedMatrixN &) const;
  kaSharedMatrixN   operator*(const kaSharedMatrixN &) const;
  kaSharedMatrixN & operator+=(const kaSharedMatrixN &);
  kaSharedMatrixN & operator-=(const kaSharedMatrixN &);
  kaSharedMatrixN & operator*=(const kaSharedMatrixN &);
  void              XMatrixN(const kaSharedMatrixN &);
  void              MatrixNX(const kaSharedMatrixN &);

  //! Write dimension and matrix elements to file
  virtual void Save(FILE *) const;

  //! Read dimension and matrix elements from file
  virtual void Restore(FILE *);

  //! Create new matrix by allocation of memory
  void New(unsigned int dimension);

  //! Delete matrix by freeing of memory
  void Delete();

  //! Add matrix memory size to pointer
  void GetSize(unsigned int &p, unsigned int dimension) const;

  //! Attach matrix to memory block
  void Attach(char * &, unsigned int);

  //! Detach matrix from memory block
  void Detach();
  int  Inverse();
  void Rotate(T, T, T, T);
  void XTranslate(T, T, T);
  void XScale(T, T, T);
  void XRotate(T, T, T, T);
  void XFrustum(T, T, T, T, T, T);
  void XOrtho(T, T, T, T, T, T);
  void TranslateX(T, T, T);
  void ScaleX(T, T, T);
  void RotateX(T, T, T, T);
  void FrustumX(T, T, T, T, T, T);
  void OrthoX(T, T, T, T, T, T);
  void Viewport(unsigned int x, unsigned int y);
};  // class kaSharedMatrixN


template<class T>
inline kaSharedMatrixN<T>::kaSharedMatrixN(unsigned int dimension) {
#if KALATTICEDEBUG
  cerr << "kaSharedMatrixN<X>::kaSharedMatrixN(unsigned int dimension)" << endl;
#endif  // if KALATTICEDEBUG
  pn        = new unsigned int;
  *pn       = dimension;
  pElements = new T[dimension*dimension];
}

template<class T>
inline kaSharedMatrixN<T>::kaSharedMatrixN(const kaSharedMatrixN<T> &rm) {
#if KALATTICEDEBUG
  cerr << "kaSharedMatrixN<X>::kaSharedMatrixN(const kaSharedMatrixN<T> & rm)" << endl;
#endif  // if KALATTICEDEBUG
  pn        = NULL;
  pElements = NULL;
  unsigned int n = rm.Dimension();
  New(n);

  T *ps = pElements, *pe = pElements+n*n;
  const T *pss = rm.pElements;

  while (ps < pe)
    *ps++ = *pss++;
}

template<class T>
inline kaSharedMatrixN<T>::~kaSharedMatrixN() {
#if KALATTICEDEBUG
  cerr << "kaSharedMatrixN<X>::~kaSharedMatrixN" << endl;
#endif  // if KALATTICEDEBUG

  // delete pn;
  // delete[] pElements;
  Delete();
}

template<class T>
inline void kaSharedMatrixN<T>::New(unsigned int dimension) {
  Delete();
  pn        = new unsigned int;
  *pn       = dimension;
  pElements = new T[dimension*dimension];
}

template<class T>
inline void kaSharedMatrixN<T>::Delete() {
  if (pn)
    delete pn;
  pn = NULL;
  if (pElements)
    delete[] pElements;
  pElements = NULL;
}

template<class T>
void kaSharedMatrixN<T>::GetSize(unsigned int &s, unsigned int dimension) const {
  s += sizeof(unsigned int);
  s += dimension*dimension*sizeof(T);
}

template<class T>
void kaSharedMatrixN<T>::Attach(char * &rp, unsigned int dimension) {
  Delete();

  pn        = (unsigned int *)rp;
  *pn       = dimension;
  rp       += sizeof(unsigned int);
  pElements = (T *)rp;
  rp       += dimension*dimension*sizeof(T);
}

template<class T>
void kaSharedMatrixN<T>::Detach() {
  pn        = (unsigned int *)NULL;
  pElements = (T *)NULL;
}

template<class T>
void kaSharedMatrixN<T>::XMatrixN(const kaSharedMatrixN<T> &rm) {
  unsigned int n = Dimension();

  assert(rm.Dimension() == n);

  kaSharedMatrixN<T> mHelp(*this);

  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++) {
      T t = 0.0;
      for (int k = 0; k < n; k++)
        t += mHelp.Get(i, k)*rm.Get(k, j);
      kaMatrix<T>::a(i, j) = t;
    }
}

template<class T>
void kaSharedMatrixN<T>::MatrixNX(const kaSharedMatrixN<T> &rm) {
  unsigned int n = Dimension();

  assert(rm.Dimension() == n);

  kaSharedMatrixN<T> mHelp(*this);

  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++) {
      T t = 0.0;
      for (int k = 0; k < n; k++)
        t += rm.Get(i, k)*mHelp.Get(k, j);
      kaMatrix<T>::a(i, j) = t;
    }
}

template<class T>
void kaSharedMatrixN<T>::Save(FILE *hFile) const {
  fwrite(pn, sizeof(unsigned int), 1, hFile);  //!< dimension of matrix
  fwrite(pElements, sizeof(T), kaMatrix<T>::Size(), hFile);  //!< elements
}

template<class T>
void kaSharedMatrixN<T>::Restore(FILE *hFile) {
  unsigned int o;

  fread(&o, sizeof(unsigned int), 1, hFile);  //!< dimension of matrix

  New(o);
  fread(pElements, sizeof(T), kaMatrix<T>::Size(), hFile);  //!< elements
}

template<class T>
kaSharedMatrixN<T> & kaSharedMatrixN<T>::operator=(const kaSharedMatrixN<T> &rm) {
  unsigned int i, j;
  unsigned int n = Dimension();

  n = (rm.Dimension() < n ? rm.Dimension() : n);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      kaMatrix<T>::a(i, j) = rm.Get(i, j);

  return *this;
}

template<class T>
inline const kaSharedMatrixN<T> kaSharedMatrixN<T>::operator*(const T f) {
  unsigned int n = Dimension();

  kaSharedMatrixN<T> temp(n);

  T *ps = temp.pElements, *pse = (temp.pElements+n*n);
  const T *pss = pElements;

  while (ps < pse)
    *ps++ = *pss++ *f;

  return temp;
}

template<class T>
inline const kaSharedMatrixN<T> kaSharedMatrixN<T>::operator/(const T f) {
  return *this*(1/f);
}

template<class T>
inline kaSharedMatrixN<T> & kaSharedMatrixN<T>::operator*=(const T f) {
  unsigned int n = Dimension();
  T *ps = pElements, *pse = (pElements+n*n);

  while (ps < pse)
    *ps++ *= f;

  return *this;
}

template<class T>
inline kaSharedMatrixN<T> & kaSharedMatrixN<T>::operator/=(const T f) {
  return *this *= (1/f);
}

template<class T>
kaSharedMatrixN<T> kaSharedMatrixN<T>::operator+(const kaSharedMatrixN<T> &rm) const {
  unsigned int n = Dimension();

  kaSharedMatrixN<T> temp(n);

  unsigned int i, j;
  n = (rm.Dimension() < n) ? rm.Dimension() : n;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      temp.a(i, j) = kaMatrix<T>::Get(i, j)+rm.Get(i, j);

  return temp;
}

template<class T>
kaSharedMatrixN<T> kaSharedMatrixN<T>::operator-(const kaSharedMatrixN<T> &rm) const {
  unsigned int n = Dimension();

  kaSharedMatrixN<T> temp(n);

  unsigned int i, j;
  n = (rm.Dimension() < n) ? rm.Dimension() : n;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      temp.a(i, j) = kaMatrix<T>::Get(i, j)-rm.Get(i, j);

  return temp;
}

template<class T>
kaSharedMatrixN<T> kaSharedMatrixN<T>::operator*(const kaSharedMatrixN<T> &rm) const {
  unsigned int n = Dimension();

  n = (rm.Dimension() < n) ? rm.Dimension() : n;
  kaSharedMatrixN<T> temp(n);

  unsigned int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      T val = 0;
      for (k = 0; k < n; k++)
        val += kaMatrix<T>::Get(i, k)*rm.Get(k, j);
      temp.a(i, j) = val;
    }

  return temp;
}

template<class T>
kaSharedMatrixN<T> & kaSharedMatrixN<T>::operator+=(const kaSharedMatrixN<T> &rm) {
  unsigned int n = Dimension();

  n = (rm.Dimension() < n) ? rm.Dimension() : n;
  unsigned int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      kaMatrix<T>::a(i, j) += rm.Get(i, j);
  return *this;
}

template<class T>
kaSharedMatrixN<T> & kaSharedMatrixN<T>::operator-=(const kaSharedMatrixN<T> &rm) {
  unsigned int n = Dimension();
  unsigned int i, j;

  n = (rm.Dimension() < n) ? rm.Dimension() : n;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      kaMatrix<T>::a(i, j) -= rm.Get(i, j);
  return *this;
}

template<class T>
kaSharedMatrixN<T> &  kaSharedMatrixN<T>::operator*=(const kaSharedMatrixN<T> &rm) {
  unsigned int n = Dimension();
  unsigned int i, j, k;

  n = (rm.Dimension() < n) ? rm.Dimension() : n;

  kaSharedMatrixN<T> mHelp(*this);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      T t = 0.0;
      for (k = 0; k < n; k++)
        t += mHelp.Get(i, k)*rm.Get(k, j);
      kaMatrix<T>::a(i, j) = t;
    }

  return *this;
}

template<class T>
int kaSharedMatrixN<T>::Inverse() {
  unsigned int n = Dimension();
  unsigned int i, j, k, m;
  T t;

  kaSharedMatrixN<T> InverseMatrix(n);
  InverseMatrix.Identity();


  for (i = 0; i < n; i++) {
    m = i;
    for (j = i+1; j < n; j++)
      if (fabs(kaMatrix<T>::Get(j, i)) > fabs(kaMatrix<T>::Get(m, i)))
        m = j;

    for (k = i; k < n; k++) {
      t                    = kaMatrix<T>::Get(i, k);
      kaMatrix<T>::a(i, k) = kaMatrix<T>::Get(m, k);
      kaMatrix<T>::a(m, k) = t;
    }

    for (k = 0; k < n; k++) {
      t                     = InverseMatrix.Get(i, k);
      InverseMatrix.a(i, k) = InverseMatrix.Get(m, k);
      InverseMatrix.a(m, k) = t;
    }

    for (j = i+1; j < n; j++) {
      if (kaMatrix<T>::Get(i, i) == 0.0)
        return 0;

      t = kaMatrix<T>::Get(j, i)/kaMatrix<T>::Get(i, i);

      for (k = i; k < n; k++)
        kaMatrix<T>::a(j, k) -= kaMatrix<T>::Get(i, k)*t;

      for (k = 0; k < n; k++)
        InverseMatrix.a(j, k) -= InverseMatrix.Get(i, k)*t;
    }
  }
  i = n;

  do {
    i--;
    if (kaMatrix<T>::Get(i, i) == 0.0)
      return 0;

    for (k = 0; k < n; k++)
      InverseMatrix.a(i, k) /= kaMatrix<T>::Get(i, i);
    kaMatrix<T>::a(i, i) = 1.0;

    for (j = 0; j < i; j++) {
      t                    = kaMatrix<T>::Get(j, i);
      kaMatrix<T>::a(j, i) = 0.0;

      for (k = 0; k < n; k++)
        InverseMatrix.a(j, k) -= InverseMatrix.Get(i, k)*t;
    }
  } while (i);

  *this = InverseMatrix;
  return 1;
}  // >::Inverse

template<class T> void kaSharedMatrixN<T>::XTranslate(T x, T y, T z) {
  kaSharedMatrixN<T> temp(4);
  temp.Translate(x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::XScale(T x, T y, T z) {
  kaSharedMatrixN<T> temp(4);
  temp.Scale(x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::XRotate(T a, T x, T y, T z) {
  kaSharedMatrixN<T> temp(4);
  temp.Rotate(a, x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::XFrustum(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  kaSharedMatrixN<T> temp(4);
  temp.Frustum(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::XOrtho(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  kaSharedMatrixN<T> temp(4);
  temp.Ortho(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::TranslateX(T x, T y, T z) {
  kaSharedMatrixN<T> temp(*this);
  kaMatrixMathN<T>::Translate(x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::ScaleX(T x, T y, T z) {
  kaSharedMatrixN<T> temp(*this);
  this->Scale(x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::RotateX(T a, T x, T y, T z) {
  kaSharedMatrixN<T> temp(*this);
  Rotate(a, x, y, z);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::FrustumX(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  kaSharedMatrixN<T> temp(*this);
  Frustum(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T> void kaSharedMatrixN<T>::OrthoX(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  kaSharedMatrixN<T> temp(*this);
  Ortho(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T>
void kaSharedMatrixN<T>::Rotate(T a, T x, T y, T z) {
  unsigned int n = Dimension();

  assert(n == 4);
  NOT_USED(n);

  T iBetrag = 1.0/sqrt(x*x + y*y + z*z);
  x *= iBetrag;
  y *= iBetrag;
  z *= iBetrag;

  kaSharedMatrixN<T> uut(3);
  uut.a(0, 0) = x*x; uut.a(0, 1) = x*y; uut.a(0, 2) = x*z;
  uut.a(1, 0) = y*x; uut.a(1, 1) = y*y; uut.a(1, 2) = y*z;
  uut.a(2, 0) = z*x; uut.a(2, 1) = z*y; uut.a(2, 2) = z*z;

  kaSharedMatrixN<T> s(3);
  s.a(0, 0) = 0;  s.a(0, 1) = -z; s.a(0, 2) = y;
  s.a(1, 0) = z;  s.a(1, 1) = 0;  s.a(1, 2) = -x;
  s.a(2, 0) = -y; s.a(2, 1) = x;  s.a(2, 2) = 0;
  s        *= sin(a);

  kaSharedMatrixN<T> I(3);
  I.Identity();
  I -= uut;
  I *= cos(a);

  uut += s;
  uut += I;

  kaMatrix<T>::Null();
  *this                = uut;
  kaMatrix<T>::a(3, 3) = 1.0;
}  // >::Rotate

template<class T>
void kaSharedMatrixN<T>::Viewport(unsigned int x, unsigned int y) {
  double tx = (double)x/2.0;
  double ty = (double)y/2.0;
  double s  = tx < ty ? tx : ty;

  ScaleX(s, s, s);
  TranslateX(tx, ty, 0.0);
}

#endif  // ifndef KASHAREDMATRIXN_H
