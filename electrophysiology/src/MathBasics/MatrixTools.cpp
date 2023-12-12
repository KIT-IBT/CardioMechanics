/*
 * File: MatrixTools.cpp
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


#include <MatrixTools.h>
#include <algorithm>

class isort_compare {  // Cmp class for isort
  double *keys;

 public:
  inline int operator()(const int &a1, const int &a2) const {
    return fabs(*(keys+a1)) > fabs(*(keys+a2));
  }

  isort_compare(double *k) {
    keys = k;
  }
};


CMatrixStructure::CMatrixStructure(const int &maxn) {
  nz   = 0;
  MAXN = maxn;
  mstruct.resize(MAXN);
}

void CMatrixStructure::add(const int &i, const int &j) {
  if (i == j)
    return;

  set<int>::iterator iexist = mstruct[j].find(i);
  if (iexist == mstruct[j].end()) {
    pair<set<int>::iterator, bool> AddRowIndex = mstruct[j].insert(i);
    if (!AddRowIndex.second)
      throw(mt_excp("CMatrixStructure::add: cannot reserve place in the matrix structure"));
    else
      nz++;
  }
}

int CMatrixStructure::size() const {
  return sizeof(int)*(nz+mstruct.size());
}

void CMatrixStructure::Print(ostream &s) {
  vector<set<int>>::iterator j = mstruct.begin();
  int coln                     = 0;
  while (j != mstruct.end()) {
    s<<"Col: "<< coln <<" Rows: ";
    set<int>::iterator i = j->begin();
    while (i != j->end()) {
      s<<*i<<' ';
      i++;
    }
    s<<'\n';
    coln++;
    j++;
  }
}

CMatrix::CMatrix(const int &maxn, const int &nz) : CDIAG(maxn), COFDIAG(nz), CCLPTR(maxn+1), CRNUM(nz) {
  MAXN = maxn;
  NZ   = nz;
}

CMatrix::CMatrix(CMatrixStructure &ms) : CDIAG(ms.MAXN), COFDIAG(ms.nz), CCLPTR(ms.MAXN+1), CRNUM(ms.nz) {
  if (ms.size() == 0)
    throw mt_excp("CMatrix::empty CMatrixStructure! Fill it first.");
  MAXN = ms.MAXN;
  NZ   = ms.nz;
  int nzIndex = 0;
  for (int j = 0; j < MAXN; j++) {
    std::set<int> &entry        = ms.mstruct[j];
    std::set<int>::iterator row = entry.begin();
    CCLPTR[j] = nzIndex;
    while (row != entry.end()) {
      CRNUM[nzIndex] = *row;
      nzIndex++;
      row++;
    }
  }
  CCLPTR[MAXN] = NZ;
}

CMatrix::CMatrix(const CMatrix &a) : CDIAG(a.MAXN), COFDIAG(a.NZ), CCLPTR(a.MAXN+1), CRNUM(a.NZ) {
  MAXN    = a.MAXN;
  NZ      = a.NZ;
  CDIAG   = a.CDIAG;
  CCLPTR  = a.CCLPTR;
  COFDIAG = a.COFDIAG;
  CRNUM   = a.CRNUM;
}

CMatrix & CMatrix::operator=(const CMatrix &A) {
  if ((MAXN != A.MAXN) ||  (NZ != A.NZ))
    throw mt_excp("CMatrix::operator=: bad matrix size\n");
  CDIAG   = A.CDIAG;
  COFDIAG = A.COFDIAG;
  CCLPTR  = A.CCLPTR;
  CRNUM   = A.CRNUM;
  return *this;
}

int CMatrix::save(const char *name) {
  ofstream to(name, ios_base::binary);

  if (!to) {
    cerr<<"CMatrix::save: cannot open "<<name<<'\n';
    return 0;
  }
  char d[] = "CMatrixFile";
  to.write(d, sizeof("CMatrixFile"));
  to.write(reinterpret_cast<char *>(&MAXN), sizeof(int));
  to.write(reinterpret_cast<char *>(&NZ), sizeof(int));
  to.write(reinterpret_cast<char *>(&CDIAG[0]), sizeof(double)*MAXN);
  to.write(reinterpret_cast<char *>(&COFDIAG[0]), sizeof(double)*NZ);
  to.write(reinterpret_cast<char *>(&CCLPTR[0]), sizeof(int)*(MAXN+1));
  to.write(reinterpret_cast<char *>(&CRNUM[0]), sizeof(int)*NZ);
  to.close();
  return 1;
}

int CMatrix::read(const char *name) {
  ifstream from(name, ios_base::binary);

  if (!from) {
    cerr<<"CMatrix::read: cannot open "<<name<<'\n';
    return 0;
  }
  char d[20];
  from.read(d, sizeof("CMatrixFile"));
  if (strcmp(d, "CMatrixFile") != 0) {
    cerr<<name<<" -not a matrix file\n";
    return 0;
  }
  from.read((char *)&MAXN, sizeof(int));
  from.read((char *)&NZ, sizeof(int));
  cout<<"maxn: "<<MAXN<<'\n';
  cout<<"nz:   "<<NZ<<'\n';


  CDIAG.resize(MAXN);
  COFDIAG.resize(NZ);
  CCLPTR.resize(MAXN+1);
  CRNUM.resize(NZ);
  from.read(reinterpret_cast<char *>(&CDIAG[0]), sizeof(double)*MAXN);
  from.read(reinterpret_cast<char *>(&COFDIAG[0]), sizeof(double)*NZ);
  from.read(reinterpret_cast<char *>(&CCLPTR[0]), sizeof(int)*(MAXN+1));
  from.read(reinterpret_cast<char *>(&CRNUM[0]), sizeof(int)*NZ);

  from.close();
  return 1;
}  // CMatrix::read

int CMatrix::save_ascii(const char *name, const int symmetric) {
  ofstream os(name);

  if (!os)
    return 0;

  os<<setprecision(15);

  for (int j = 0; j < MAXN; j++) {
    if (CDIAG[j])
      os<<j+1<<' '<<j+1<<' '<<CDIAG[j]<<'\n';

    for (int i = CCLPTR[j]; i < CCLPTR[j+1]; i++) {
      os<<CRNUM[i]+1<<' '<<j+1<<' '<<COFDIAG[i]<<'\n';
      if (symmetric)
        os<<j+1<<' '<<CRNUM[i]+1<<' '<<COFDIAG[i]<<'\n';
    }
  }

  return 1;
}

/*** Some operators ***/
ostream & operator<<(ostream &s, CMatrix &m) {
  s<<"CDIAG: \n";
  for (int i = 0; i < m.CDIAG.size(); i++)
    s<<m.CDIAG[i]<<' ';
  s<<"\nCOFDIAG:\n";
  for (int i = 0; i < m.COFDIAG.size(); i++)
    s<<m.COFDIAG[i]<<' ';
  s<<"\nCRNUM:\n";
  for (int i = 0; i < m.CRNUM.size(); i++)
    s<<m.CRNUM[i]<<' ';
  s<<"\nCCLPTR:\n";
  for (int i = 0; i < m.CCLPTR.size(); i++)
    s<<m.CCLPTR[i]<<' ';
  return s<<'\n';
}

/***Some useful functions ****/

void SetColToZero(CMatrix &m, const int &col) {
  if (col > m.MAXN-1)
    throw mt_excp("CMatrix::SetColToZero: number of column greater than MAXN-1");
  int nb = m.CCLPTR[col];
  int ne = m.CCLPTR[col+1];
  for (int i = nb; i < ne; i++) m.COFDIAG[i] = 0.0;
  m.CDIAG[col] = 0.0;
}

void SetRowToZero(CMatrix &m, const int &row) {
  if (row < m.MAXN)
    m.CDIAG[row] = 0.0;
  for (int j = 0; j < m.MAXN; j++) {
    int nb = m.CCLPTR[j];
    int ne = m.CCLPTR[j+1];
    if (nb != ne) {
      while (ne-nb > 1) {  // bisection to find row
        int pos = (nb+ne)/2;
        if (m.CRNUM[pos] > row)
          ne = pos; else
          nb = pos;
      }
      if (m.CRNUM[nb] == row)
        m.COFDIAG[nb] = 0.0;
    }
  }
}

int JPICC(CMatrix &A) {
  /*
     Jones/Plassman incomplete Cholesky: column oriented
     ACM Transactions on Mathematical Software, Vol 21, No 1, March 1995, pp. 5-17
     FORTRAN version is available via netlib
   */


  const int N = A.MAXN;  // dimensionality of input matrix

  valarray<double> ta(N);  // a vector to keep the current column
  int talen = 0;  // actual size of ta

  valarray<int> itcol(N);  // keeps the row values of the current column

  valarray<int> ifirst(-1, N);  /* ifirst(j) points to the next value in  column j to use.
                                   ifirst also has a dual use. At step K, only  the first K-1 elements are used for the
                                   above purpose. For the last N-K elements, ifirst(j) indicates if a nonzero value
                                   exists in the position j of column K */
  valarray<int> List(-1, N);  /* List(j) points to the linked list of columns which will update column j */

  /* here begins ... */
  for (int k = 0; k < N; k++) {  // over all columns
    // load column k into ta
    talen = 0;
    int isk = A.CCLPTR[k];
    int iek = A.CCLPTR[k+1];
    for (int j = isk; j < iek; j++) {
      int row = A.CRNUM[j];
      ta[row]      = A.COFDIAG[j];
      itcol[talen] = row;
      ifirst[row]  = 0;
      talen++;
    }

    // take the sqrt of the diagonal of column k
    if (A.CDIAG[k] < 0)
      return -k; // factorization failed
    else
      A.CDIAG[k] = sqrt(A.CDIAG[k]);

    // update column k using the previous columns
    int j = List[k];
    while (j != -1) {
      int isj     = ifirst[j];
      int iej     = A.CCLPTR[j+1]-1;
      double lval = A.COFDIAG[isj];
      isj++;
      if (isj < iej) {
        ifirst[j] = isj;
        int iptr = j;
        j                  = List[j];
        List[iptr]         = List[A.CRNUM[isj]];
        List[A.CRNUM[isj]] = iptr;
      } else {j = List[j];}

      for (int i = isj; i <= iej; i++) {
        int row = A.CRNUM[i];
        if (ifirst[row] != -1) {ta[row] -= lval*A.COFDIAG[i];} else {
          ifirst[row]  = 0;
          itcol[talen] = row;
          ta[row]      = (-1.0)*lval*A.COFDIAG[i];
          talen++;
        }
      }
    }


    // update remaining diagonals using column k
    for (int j = 0; j < talen; j++) {
      int row = itcol[j];
      ta[row]      /= A.CDIAG[k];
      A.CDIAG[row] -= ta[row]*ta[row];
    }


    // find the lagest elements in column k
    int klen  = iek-isk;
    int count = (klen < talen) ? klen : talen;  // k=min(klen,talen)
    isort(talen, count, ta, itcol);

    // sort itcol
    int *first = &itcol[0];
    int *last  = first+count;
    sort(first, last);

    // put the lagest elements back into sparse data structure
    count = 0;
    for (int j = isk; j < iek; j++) {
      A.COFDIAG[j] = ta[itcol[count]];
      A.CRNUM[j]   = itcol[count];
      count++;
    }

    // ifirst and List keep track of where in column k we are
    if (isk < (iek-1)) {
      int iptr = A.CRNUM[isk];
      List[k]    = List[iptr];
      List[iptr] = k;
      ifirst[k]  = isk;
    }
    for (int j = 0; j < talen; j++) {
      ifirst[itcol[j]] = -1;
    }
  }
  return 0;
}  // JPICC

void ibsort(const int &n, const int &k, const valarray<double> &akeys, valarray<int> &indvec) {
  // internal variables
  int i, j;
  int itemp, curptr, right, left;
  double curmin, newval, curval, lval;

  // if the list is small or the number required is 0 then return
  if ((n <= 1) || (k <= 0))
    return;

  // heap sort the first k elements of the vector
  ihsort(k, indvec, akeys);

  // loop through the rest of the vector and find any elements that are lager than any
  // of the first k elements

  curmin = fabs(akeys[indvec[k-1]]);
  for (i = k; i < n; i++) {
    itemp  = indvec[i];
    newval = fabs(akeys[itemp]);
    if (newval > curmin) {
      // find position for the new value
      left = 0;
      lval = fabs(akeys[indvec[0]]);
      if (newval > lval) {
        curptr = 0;
        goto lbl2;
      }
      right  = k-1;
      curptr = k/2;
lbl1: if (right > left+1) {
        curval = fabs(akeys[indvec[curptr]]);
        if (curval < newval) {right = curptr;} else {
          left = curptr;
          lval = curval;
        }
        curptr = (right+left)/2;
        goto lbl1;
      }
      curptr = right;

      // shift sorted values and insert new value
lbl2: indvec[i] = indvec[k-1];
      for (j = k-1; j > curptr; j--) indvec[j] = indvec[j-1];
      indvec[curptr] = itemp;
      curmin         = fabs(akeys[indvec[k-1]]);
    }
  }
}  // ibsort

void ihsort(const int &len, valarray<int> &indvec, const valarray<double> &akeys) {
  int k, m, lheap, rheap, mid;
  int x;

  if (len <= 1)
    return;

  // build the heap
  mid = (len-1)/2;
  for (k = mid; k >= 0; k--) {
    x     = indvec[k];
    lheap = k;
    rheap = len-1;
    m     = lheap*2;
lbl1: if (m > rheap) {
      indvec[lheap] = x;
      goto lbl2;
    }
    if (m < rheap) {
      if (fabs(akeys[indvec[m]]) > fabs(akeys[indvec[m+1]]))
        m++;
    }
    if (fabs(akeys[x]) <= fabs(akeys[indvec[m]])) {
      m = rheap+1;
    } else {
      indvec[lheap] = indvec[m];
      lheap         = m;
      m             = 2*lheap;
    }
    goto lbl1;
lbl2:;
  }

  // sort the heap
  for (k = len-1; k >= 1; k--) {
    x         = indvec[k];
    indvec[k] = indvec[0];
    lheap     = 0;
    rheap     = k-1;
    m         = 1;
lbl4: if (m > rheap) {
      indvec[lheap] = x;
      goto lbl5;
    }
    if (m < rheap) {
      if (fabs(akeys[indvec[m]]) > fabs(akeys[indvec[m+1]]))
        m++;
    }
    if (fabs(akeys[x]) <= fabs(akeys[indvec[m]])) {
      m = rheap+1;
    } else {
      indvec[lheap] = indvec[m];
      lheap         = m;
      m             = 2*lheap;
    }
    goto lbl4;
lbl5:;
  }
}  // ihsort

void isort(const int &n, const int &k, valarray<double> &akeys, valarray<int> &indvec) {
  double *kk = &akeys[0];
  isort_compare Cmp(kk);
  int *start = &indvec[0];
  int *end   = &indvec[n-1];
  int *mid   = &indvec[k-1];

  partial_sort(start, mid+1, end+1, Cmp);
}

void CholFactSolve(const CMatrix &A, const valarray<double> &y, valarray<double> &x) {
  // no error check will be performed - think before use

  valarray<double> temp_y(y);

  /* first solve  L.z=y where z=Transp(L).x  here z overwrites x*/
  for (int i = 0; i < A.MAXN; i++) {
    x[i] = temp_y[i]/A.CDIAG[i];
    for (int k = A.CCLPTR[i]; k < A.CCLPTR[i+1]; k++)
      temp_y[A.CRNUM[k]] -= x[i]*A.COFDIAG[k];
  }

  /* now solve Transp(L).x=z */
  for (int i = A.MAXN-1; i >= 0; i--) {
    double sum = 0;
    for (int k = A.CCLPTR[i]; k < A.CCLPTR[i+1]; k++)
      sum += x[A.CRNUM[k]]*A.COFDIAG[k];
    x[i] = (x[i]-sum)/A.CDIAG[i];
  }
}

void Scale(CMatrix &A, valarray<double> &B) {
  if (A.MAXN != B.size())
    throw mt_excp("Scale: MAXN != B.size()");
  for (int i = 0; i < A.MAXN; i++) {
    if (A.CDIAG[i] != 0) {
      double factor = 1.0/A.CDIAG[i];
      A.CDIAG[i] = factor;
      B[i]      *= factor;
    } else {throw mt_excp("Scale: zero diagonal element encountered\n");}
  }
  for (int j = 0; j < A.MAXN; j++)
    for (int k = A.CCLPTR[j]; k < A.CCLPTR[j+1]; k++) {
      int row = A.CRNUM[k];
      A.COFDIAG[k] *= A.CDIAG[row];
    }
  A.CDIAG = 1.0;
}

void SymmScale(CMatrix &A, valarray<double> &B, valarray<double> &scale) {
  if (A.MAXN != B.size())
    throw mt_excp("Scale: MAXN != B.size()");
  for (int i = 0; i < A.MAXN; i++) {
    if (A.CDIAG[i] != 0) {
      double factor = 1.0/sqrt(A.CDIAG[i]);
      A.CDIAG[i] = factor;
      B[i]      *= factor;
      scale[i]   = factor;
    } else {throw mt_excp("SymmScale: zero diagonal element encountered\n");}
  }
  for (int j = 0; j < A.MAXN; j++)
    for (int k = A.CCLPTR[j]; k < A.CCLPTR[j+1]; k++) {
      int row = A.CRNUM[k];
      A.COFDIAG[k] *= A.CDIAG[row]*A.CDIAG[j];
    }
  A.CDIAG = 1.0;
}

int ISTDIC(CMatrix &A) {
  int isk, iek, isj, iej;
  int i, j, k;
  int row;
  double lval, t;
  int iptr;

  const int N = A.MAXN;  // dimensionality of input matrix

  valarray<double> ta(A.MAXN);  // a temporary work vector of length N to keep the current column
  valarray<int> ifirst(-1, A.MAXN);  /* ifirst(j) points to the next value in  column j to use.
                                        ifirst also has a dual use. At step K, only  the first K-1 elements are used for
                                           the
                                        above purpose. For the last N-K elements, ifirst(j) indicates if a nonzero value
                                        exists in the position j of column K */
  valarray<int> List(-1, A.MAXN);  /* List(j) points to the linked list of columns which will update column j */

  for (k = 0; k < N; k++) {  // over all columns
    // load column k into ta
    isk = A.CCLPTR[k];
    iek = A.CCLPTR[k+1];
    for (j = isk; j < iek; j++) {
      row         = A.CRNUM[j];
      ta[row]     = A.COFDIAG[j];
      ifirst[row] = 0;
    }

    // make sure the diagonal of k is okay and then take the sqrt
    if (A.CDIAG[k] < 0)
      return -k;
    else
      A.CDIAG[k] = sqrt(A.CDIAG[k]);

    // update column k using previous columns
    j = List[k];
    while (j != -1) {
      isj  = ifirst[j];
      iej  = A.CCLPTR[j+1]-1;
      lval = A.COFDIAG[isj];
      isj++;
      if (isj < iej) {
        ifirst[j]          = isj;
        iptr               = j;
        j                  = List[j];
        List[iptr]         = List[A.CRNUM[isj]];
        List[A.CRNUM[isj]] = iptr;
      } else {
        j = List[j];
      }
      for (i = isj; i <= iej; i++) {
        row = A.CRNUM[i];
        if (ifirst[row] != -1)
          ta[row] -= lval*A.COFDIAG[i];
      }
    }

    // ifirst and List keep track of where in column k we are
    if (isk < (iek-1)) {
      iptr       = A.CRNUM[isk];
      List[k]    = List[iptr];
      List[iptr] = k;
      ifirst[k]  = isk;
    }

    // update remaining diagonals using column k
    lval = A.CDIAG[k];
    for (j = isk; j < iek; j++) {
      row          = A.CRNUM[j];
      t            = ta[row];
      ifirst[row]  = -1;
      t           /= lval;
      A.CDIAG[row] = A.CDIAG[row]-t*t;
      A.COFDIAG[j] = t;
    }
  }

  return 0;
}  // ISTDIC
