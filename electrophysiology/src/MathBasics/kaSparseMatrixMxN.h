/*
 * File: kaSparseMatrixMxN.h
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


#pragma once

#include <valarray>
#include <vector>
#include <map>

#include <kaExceptions.h>
#include <kaMatrixBase.h>
#include <kaMatrixSparseCol.h>
#include <kaSharedMatrixMxN.h>
#include <kaBasicIO.h>

using namespace nskaGlobal;
template<class T>
class kaSparseMatrixMxN : public kaMatrixBase<T> {
 public:
  typedef std::valarray<std::map<size_t, T>> sparse_data_t;

  kaSparseMatrixMxN() : m_numberOfNonZeros(0), m_numberOfRows(0), m_numberOfCols(0) {}

  kaSparseMatrixMxN(const kaMatrixSparseCol<T> &systMatrix, const std::valarray<size_t> &listOfRows,
                    const std::valarray<size_t> &listOfColumns);

  kaSparseMatrixMxN(const kaSparseMatrixMxN<T> &matrix);

  kaSparseMatrixMxN(size_t cols, const sparse_data_t &data);

  //  ~kaSparseMatrixMxN();
  void assignData(size_t cols, const sparse_data_t &data);

  virtual unsigned int NumCols() const {return m_numberOfCols;}

  virtual unsigned int NumRows() const {return m_numberOfRows;}

  virtual const T Get(int row, int col) const;

  kaSharedMatrixMxN<T> *selfMultiplyByMatrix(kaMatrix<T> *matrix);  // returns matrix*this
  kaSharedMatrixMxN<T> *matrixMultiplyBySelf(const kaMatrix<T> *matrix) const;  // returns this*matrix; NOT TESTED!
  void                  matrixVectorProduct(const std::valarray<T> &x, std::valarray<T> &result) const;
  void                  exportASCII(const char *filename); // exports into ascii with columns containing row id, column
                                                           // id and value"
  virtual void Save(FILE *) const; // NOT TESTED!
  virtual void Restore(FILE *);    // NOT TESTED!

 protected:
  size_t m_numberOfNonZeros;
  size_t m_numberOfRows;
  size_t m_numberOfCols;

  std::valarray<T> m_values;          // length = m_numberOfNonZeros
  std::valarray<size_t> m_columns;    // length = m_numberOfNonZeros
  std::valarray<size_t> m_rowOffsets;  // length = m_numberOfRows+1

 private:
  virtual T *Pointer() {}

  virtual const T *Pointer() const {}

  virtual T *operator[](int i) {}

  virtual const T *operator[](int i) const {}

  virtual T & a(int, int) {}

  virtual T & a(int) {}

  virtual const T Get(int p) const {}

  virtual void Print(const char *format = "%7.3lf\t") const {}
};  // class kaSparseMatrixMxN

// debug: test the case if some elements in listOfRows or listOfColumns exceed the size of the systMatrix
template<class T>
kaSparseMatrixMxN<T>::kaSparseMatrixMxN(const kaMatrixSparseCol<T> &systMatrix, const std::valarray<size_t> &listOfRows,
                                        const std::valarray<size_t> &listOfColumns) {
  m_numberOfRows = listOfRows.size();
  m_numberOfCols = listOfColumns.size();

  // create reverse mask arrays: systMatrix -> this
  std::valarray<int> idxRows(-1, systMatrix.Dimension());
  std::valarray<int> idxCols(-1, systMatrix.Dimension());

  for (size_t i = 0; i < m_numberOfRows; i++) idxRows[listOfRows[i]] = (int)i;
  for (size_t i = 0; i < m_numberOfCols; i++) idxCols[listOfColumns[i]] = (int)i;

  // compute the number of non-zero elements: overall and in each row
  m_numberOfNonZeros = 0;

  std::valarray<size_t> elementsInRow((size_t)0, m_numberOfRows);

  sparse_data_t matrixStructure(m_numberOfRows);  // structure of the future matrix

  for (size_t i = 0; i < systMatrix.Dimension(); i++) {
    if ((idxRows[i] >= 0) && (idxCols[i] >= 0)) {
      m_numberOfNonZeros++;
      elementsInRow[idxRows[i]]++;
      matrixStructure[idxRows[i]][idxCols[i]] = systMatrix.CDIAG[i];
    }
    for (int j = systMatrix.CCLPTR[i]; j < systMatrix.CCLPTR[i+1]; j++) {
      if ((idxRows[i] >= 0) && (idxCols[systMatrix.CRNUM[j]] >= 0)) {
        m_numberOfNonZeros++;
        elementsInRow[idxRows[i]]++;
        matrixStructure[idxRows[i]][idxCols[systMatrix.CRNUM[j]]] = systMatrix.COFDIAG[j];
      }
      if ((idxRows[systMatrix.CRNUM[j]] >= 0) && (idxCols[i] >= 0)) {
        m_numberOfNonZeros++;
        elementsInRow[idxRows[systMatrix.CRNUM[j]]]++;
        matrixStructure[idxRows[systMatrix.CRNUM[j]]][idxCols[i]] = systMatrix.COFDIAG[j];
      }
    }
  }

  // create the data
  m_values.resize(m_numberOfNonZeros);
  m_columns.resize(m_numberOfNonZeros);
  m_rowOffsets.resize(m_numberOfRows+1);

  m_rowOffsets[0] = 0;
  for (size_t row = 0; row < m_numberOfRows; row++) {
    m_rowOffsets[row+1] = m_rowOffsets[row] + elementsInRow[row];
    size_t counter = 0;
    for (typename std::map<size_t, T>::iterator mit = matrixStructure[row].begin(); mit != matrixStructure[row].end();
         mit++) {
      size_t offset = m_rowOffsets[row]+(counter++);
      m_columns[offset] = mit->first;
      m_values[offset]  = mit->second;
    }
  }

  std::cerr << "kaSparseMatrixMxN<T>::kaSparseMatrixMxN: created with " << m_numberOfRows << "rows and " <<
    m_numberOfCols << "columns\n";
}

template<class T>
kaSparseMatrixMxN<T>::kaSparseMatrixMxN(const kaSparseMatrixMxN<T> &matrix) {
  m_numberOfNonZeros = matrix.m_numberOfNonZeros;
  m_numberOfRows     = matrix.m_numberOfRows;
  m_numberOfCols     = matrix.m_numberOfCols;

  m_values.resize(matrix.m_values.size());
  m_columns.resize(matrix.m_columns.size());
  m_rowOffsets.resize(matrix.m_rowOffsets.size());

  m_values     = matrix.m_values;
  m_columns    = matrix.m_columns;
  m_rowOffsets = matrix.m_rowOffsets;
}

//! constructor creating a new sparse matrix from the number of columns and a sparse_data_t object.
//! Number of rows is taken from the size of the sparse_data_t.
template<class T>
kaSparseMatrixMxN<T>::kaSparseMatrixMxN(size_t cols, const sparse_data_t &data) {
  assignData(cols, data);
}

template<class T> void kaSparseMatrixMxN<T>::assignData(size_t               cols,
                                                        const sparse_data_t &data) {
  m_numberOfRows     = data.size();
  m_numberOfCols     = cols;
  m_numberOfNonZeros = 0;

  for (size_t i = 0; i < m_numberOfRows; i++) m_numberOfNonZeros += data[i].size();

  // create the data
  m_values.resize(m_numberOfNonZeros);
  m_columns.resize(m_numberOfNonZeros);
  m_rowOffsets.resize(m_numberOfRows+1);

  m_rowOffsets[0] = 0;
  for (size_t row = 0; row < m_numberOfRows; row++) {
    m_rowOffsets[row+1] = m_rowOffsets[row] + data[row].size();
    size_t counter = 0;
    for (typename std::map<size_t, T>::const_iterator mit = data[row].begin(); mit != data[row].end(); mit++) {
      size_t offset = m_rowOffsets[row]+(counter++);
      m_columns[offset] = mit->first;
      m_values[offset]  = mit->second;
    }
  }
}

template<class T> const T kaSparseMatrixMxN<T>::Get(int row, int col) const {
  if ((row < 0) || (row >= m_numberOfRows))
    throw kaBaseException("kaSparseMatrixMxN<T>::Get: wrong row index %d", row);
  if ((col < 0) || (col >= m_numberOfCols))
    throw kaBaseException("kaSparseMatrixMxN<T>::Get: wrong col index %d", col);
  int ind_start = m_rowOffsets[row];
  int ind_end   = m_rowOffsets[row+1];

  for (int i = ind_start; i < ind_end; i++) {
    if (m_columns[i] == col)
      return m_values[i];
    else if (m_columns[i] > col)
      return (T)0;
  }

  return (T)0;
}

template<class T> kaSharedMatrixMxN<T> *kaSparseMatrixMxN<T>::selfMultiplyByMatrix(kaMatrix<T> *matrix) { // returns
                                                                                                          // matrix*this
  if (!matrix)
    throw kaBaseException("kaSparseMatrixMxN::selfMultiplyByMatrix: no matrix");
  if (matrix->NumCols() != m_numberOfRows)
    throw kaBaseException("kaSparseMatrixMxN::selfMultiplyByMatrix: wrong matrix size");

  kaSharedMatrixMxN<T> *result = new kaSharedMatrixMxN<T>(matrix->NumRows(), m_numberOfCols);
  result->Null();
  for (size_t row = 0; row < m_numberOfRows; row++)
    for (size_t col_idx = m_rowOffsets[row]; col_idx < m_rowOffsets[row+1]; col_idx++)
      for (int i = 0; i < matrix->NumRows(); i++)
        result->a(i, m_columns[col_idx]) += matrix->Get(i, row) * m_values[col_idx];

  return result;
}

template<class T> kaSharedMatrixMxN<T> *kaSparseMatrixMxN<T>::matrixMultiplyBySelf(const kaMatrix<T> *matrix) const// returns
                                                                                                                   // this*matrix
{
  if (!matrix)
    throw kaBaseException("kaSparseMatrixMxN::matrixMultiplyBySelf: no matrix");
  if (matrix->NumRows() != m_numberOfCols)
    throw kaBaseException("kaSparseMatrixMxN::matrixMultiplyBySelf: wrong matrix size");

  kaSharedMatrixMxN<T> *result = new kaSharedMatrixMxN<T>(m_numberOfRows, matrix->NumCols());
  result->Null();
  for (size_t row = 0; row < m_numberOfRows; row++)
    for (size_t col_idx = m_rowOffsets[row]; col_idx < m_rowOffsets[row+1]; col_idx++)
      for (int i = 0; i < matrix->NumCols(); i++)
        result->a(row, i) += m_values[col_idx] * matrix->Get(m_columns[col_idx], i);

  return result;
}

template<class T> void kaSparseMatrixMxN<T>::matrixVectorProduct(const std::valarray<T> &x,
                                                                 std::valarray<T>       &result) const {
  if (x.size() != m_numberOfCols)
    throw kaBaseException("kaSparseMatrixMxN<T>::matrixVectorProduct: wrong vector size");
  result.resize(m_numberOfRows);
  result = 0;

  for (size_t row = 0; row < m_numberOfRows; row++)
    for (size_t col_idx = m_rowOffsets[row]; col_idx < m_rowOffsets[row+1]; col_idx++)
      result[row] += m_values[col_idx] * x[m_columns[col_idx]];
}

template<class T> void kaSparseMatrixMxN<T>::exportASCII(const char *filename) {
  std::ofstream fd(filename);

  if (!fd)
    throw kaBaseException("kaSparseMatrixMxN::exportASCII: error opening file %s", filename);

  for (size_t row = 0; row < m_numberOfRows; row++)
    for (size_t col_idx = m_rowOffsets[row]; col_idx < m_rowOffsets[row+1]; col_idx++)
      fd << row << " " << m_columns[col_idx] << " " << m_values[col_idx] << "\n";
}

template<class T> void kaSparseMatrixMxN<T>::Save(FILE *fd) const {
  kaWrite("kaSparseMatrixMxN", sizeof("kaSparseMatrixMxN"), fd);
  kaWrite(m_numberOfNonZeros, fd);
  kaWrite(m_numberOfRows, fd);
  kaWrite(m_numberOfCols, fd);

  // Stroustrup says it's OK (section 22.4.3). TO TEST!
  kaWrite(&m_values[0], m_values.size(), fd);
  kaWrite(&m_columns[0], m_columns.size(), fd);
  kaWrite(&m_rowOffsets[0], m_rowOffsets.size(), fd);
}

template<class T> void kaSparseMatrixMxN<T>::Restore(FILE *fd) {
  char buf[32];

  kaRead(buf, sizeof("kaSparseMatrixMxN"), fd);
  if (strcmp(buf, "kaSparseMatrixMxN"))
    throw kaBaseException("kaSparseMatrixMxN::Restore: not a kaSparseMatrixMxN file");
  kaRead(&m_numberOfNonZeros, fd);
  kaRead(&m_numberOfRows, fd);
  kaRead(&m_numberOfCols, fd);

  m_values.resize(m_numberOfNonZeros);
  m_columns.resize(m_numberOfNonZeros);
  m_rowOffsets.resize(m_numberOfRows+1);

  // Stroustrup says it's OK (section 22.4.3). TO TEST!
  kaRead(&m_values[0], m_values.size(), fd);
  kaRead(&m_columns[0], m_columns.size(), fd);
  kaRead(&m_rowOffsets[0], m_rowOffsets.size(), fd);
}
