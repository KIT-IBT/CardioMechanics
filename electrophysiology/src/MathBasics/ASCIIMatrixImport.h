/*! \file ASCIIMatrixImport.h
   \brief Import matrices in ASCII format, understand comments (with '#').

   os May 18 2001
   os Mar 23 2003 modification

   The import subroutines understand comments ('#') in ascii files

   helper functions to read matrices form ASCII files
   (e. g. MATLAB .dat-files)

 */

#ifndef __ASCIIMatrixImport_h
#define __ASCIIMatrixImport_h

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

#include <kaExceptions.h>
#include <kaMatrix.h>

class ascii_import_excp : public kaBaseException {
 public:
  ascii_import_excp(const char *msg1, const char *msg2 = "") : kaBaseException("exception in %s: %s %s", __FILE__, msg1,
                                                                               msg2) {}
};

//! Analyse ASCII file to get matrix dimensions
void GetASCIIMatrixDimensions(const char *name, int &rows, int &cols);

//! Simply reads the text file.
template<class T> void ReadASCIIMatrix(const char *name, kaMatrix<T> &m, const int rows, const int cols);

// ====== Implementation ======
void GetASCIIMatrixDimensions(const char *name, int &rows, int &cols) {
  ifstream f(name);

  if (!f)
    throw ascii_import_excp("cannot open file ", name);


  int numrows(0);
  int numcols(0);

  string buf;
  bool   firstrow = 1;

  while (getline(f, buf)) {
    string::iterator pcomm = find(buf.begin(), buf.end(), '#');
    buf.erase(pcomm, buf.end());

    istringstream iss(buf);

    char c;
    while (iss.get(c)) // skip first spaces
      if (!isspace(c)) {
        iss.putback(c);
        break;
      }

    // determine number of columns
    if (iss) {
      if (firstrow) {
        double a;
        while (iss>>a) {numcols++;}
        if (!numcols)
          throw ascii_import_excp("cannot get any columns from the ascii file: ", name);
        firstrow = false;
      }
      numrows++;
    }
  }

  if (!numrows)
    throw ascii_import_excp("cannot read any rows from the ascii file :", name);

  rows = numrows;
  cols = numcols;
}  // GetASCIIMatrixDimensions

// Description:
// Simply reads the text file.
template<class T> void ReadASCIIMatrix(const char *name, kaMatrix<T> &m, const int rows, const int cols) {
  ifstream f(name);

  if (!f)
    throw ascii_import_excp("cannot open :", name);

  int i = 0, j = 0;  // row, col

  string buf;

  while (getline(f, buf)) {
    string::iterator pcomm = find(buf.begin(), buf.end(), '#');
    buf.erase(pcomm, buf.end());

    istringstream iss(buf);

    char c;
    while (iss.get(c)) // skip first spaces
      if (!isspace(c)) {
        iss.putback(c);
        break;
      }


    if (iss) {  // read the row
      double a;
      j = 0;
      while (iss>>a) {
        m.a(i, j) = a;
        j++;
      }
      i++;
    }
    if (i > rows)
      throw ascii_import_excp("more rows are present as allocated");
    if (j > cols)
      throw ascii_import_excp("more columns are present as allocated");
  }

  /*   int k=0; */
  /*   const int msize = m.Size(); */
  /*   for(int i=0; i<rows; i++){ */
  /*     for(int j=0; j<cols; j++){ */
  /*       T val; */
  /*       f>>val; */
  /*       // if(!f.good()) throw ascii_import_excp("error reading matrix"); */
  /*       if(k >= msize) throw ascii_import_excp("not enough place to read matrix"); */
  /*       m.a(i,j)=val; */
  /*       k++; */
  /*     } */
  /*   } */
}  // ReadASCIIMatrix

#endif  // ifndef __ASCIIMatrixImport_h
