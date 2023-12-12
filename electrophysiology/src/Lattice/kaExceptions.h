/*
 * File: kaExceptions.h
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


#ifndef KAEXCEPTIONS_H
#define KAEXCEPTIONS_H

#include <kaMachineOS.h>

//! Base and standard class for exception handling
/*!
   An object of class kaBaseException and derived classes is typically thrown in case of exception.
   The object includes a string describing details of the exception's cause.
   \n\n
   \author fs, IBT - UniversitÃ¤t Karlsruhe (TH)
 */

class kaBaseException {
  char errstr[256];  //!< Text describing details of the exception's cause.

 public:
  //! Default constructor
  kaBaseException() {
    errstr[0] = 0;
  }

  //! Constructor with arguments similiar to printf
  kaBaseException(const char *format, ...) {
    va_list args;

    va_start(args, format);
    vsnprintf(errstr, sizeof(errstr), format, args);
    va_end(args);
  }

  //! Get exception text
  inline const char *getText() {return errstr;}

  friend ostream & operator<<(ostream &s, const kaBaseException &e);
};

//! Output operator for class kaBaseException. The errstr is printed.
inline ostream & operator<<(ostream &s, const kaBaseException &e) {
  return s << e.errstr;
}

namespace nskaGlobal {
//! Class for handling of missing options in arguments of a program via exception
class MissingOption : public kaBaseException {
 public:
  MissingOption(const char *call) : kaBaseException("Command line option %s is missing", call) {}
};


//! Class for handling of invalid options in arguments of a program via exception
class IncorrectOption : public kaBaseException {
 public:
  IncorrectOption(const char *call) : kaBaseException("Incorrect use of command line option %s", call) {}
};


//! Class for handling of a dimension mismatch of objects via exception
/*!
   A dimension mismatch is found, if x-, y- or z- dimensions differ.
   Typically, the class is used to check objects of classes kaLattice, kaLatticeHeader etc..
 */
class DimensionMismatch : public kaBaseException {
 public:
  DimensionMismatch(const char *d1) : kaBaseException("Lattice dimension mismatch %s", d1) {}

  DimensionMismatch(const char *d1, const char *d2) : kaBaseException("Lattice dimension mismatch %s <-> %s", d1, d2) {}

  DimensionMismatch(const char *d1, int x, int y, int z) : kaBaseException("Lattice dimension mismatch %s <-> %d %d %d",
                                                                           d1, x, y, z) {}
};

//! Class for handling of a type mismatch of objects via exception
/*!
   A type mismatch is found, if the DataType differs.
   Typically, the class is used to check objects of classes kaLattice, kaLatticeHeader etc..
 */
class TypeMismatch : public kaBaseException {
 public:
  TypeMismatch(const char *d1) : kaBaseException("Lattice type mismatch %s", d1) {}

  TypeMismatch(const char *d1, const char *d2) : kaBaseException("Lattice type mismatch %s <-> %s", d1, d2) {}
};
}  // namespace nskaGlobal
#endif  // ifndef KAEXCEPTIONS_H
