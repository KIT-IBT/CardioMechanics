/*
 * File: kaDataTypes.h
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


#ifndef KADATATYPES_H
#define KADATATYPES_H

#include <kaMachineOS.h>

namespace nskaGlobal {
//! Enumeration for data types of lattices.
/*!
   dtLlong and dtULlong are not fully implemented yet.
   dtUnknown is used for non-standard types.

   \author cw,fs, IBT - UniversitÃ¤t Karlsruhe (TH)
 */
typedef enum {
  dtUchar   = 0,
  dtSchar   = 1,
  dtShort   = 2,
  dtUshort  = 3,
  dtFloat   = 4,
  dtInt     = 5,
  dtUint    = 6,
  dtLong    = 7,
  dtUlong   = 8,
  dtDouble  = 9,
  dtLlong   = 10,
  dtULlong  = 11,
  dtLdouble = 12,
  dtUnknown = 13
} DataType;

//! Number of data types
static const int DT_MAX_DATATYPES = 14;

//! Get data type with template type
template<class X> inline DataType GetTemplateDataType() {
  if (typeid(X) == typeid(uint8_t))
    return dtUchar;

  if (typeid(X) == typeid(int8_t))
    return dtSchar;

  if (typeid(X) == typeid(int16_t))
    return dtShort;

  if (typeid(X) == typeid(uint16_t))
    return dtUshort;

  if (typeid(X) == typeid(int32_t))
    return dtInt;

  if (typeid(X) == typeid(uint32_t))
    return dtUint;

  if (typeid(X) == typeid(int64_t))
    return dtLong;

  if (typeid(X) == typeid(uint64_t))
    return dtUlong;

  //    if (typeid(X)==typeid(int128_t)) return dtLLong;
  //    if (typeid(X)==typeid(uint128_t)) return dtULLong;
  if (typeid(X) == typeid(float))
    return dtFloat;

  if (typeid(X) == typeid(double))
    return dtDouble;

  if (typeid(X) == typeid(long double))
    return dtLdouble;

  return dtUnknown;
}  // GetTemplateDataType

//! Class for definition of data types used in lattices.
/*!
   data types, sizes and corresponding lattice file names are defined.

   \author cw,fs,idb, IBT - UniversitÃ¤t Karlsruhe (TH)
 */

class DataTypeInfo {
 public:
  //! Default constructor
  DataTypeInfo(void) {
    datasize[dtUchar]   = sizeof(uint8_t);
    datasize[dtSchar]   = sizeof(int8_t);
    datasize[dtShort]   = sizeof(int16_t);
    datasize[dtUshort]  = sizeof(uint16_t);
    datasize[dtInt]     = sizeof(int32_t);
    datasize[dtUint]    = sizeof(uint32_t);
    datasize[dtLong]    = sizeof(int64_t);
    datasize[dtUlong]   = sizeof(uint64_t);
    datasize[dtLlong]   = sizeof(int64_t) * 2;
    datasize[dtULlong]  = sizeof(uint64_t) * 2;
    datasize[dtFloat]   = sizeof(float);
    datasize[dtDouble]  = sizeof(double);
    datasize[dtLdouble] = sizeof(double) * 2;
    datasize[dtUnknown] = 0;

    strcpy(suffixlist[dtUchar],   ".lat");
    strcpy(suffixlist[dtSchar],   ".slat");
    strcpy(suffixlist[dtShort],   ".shlat");
    strcpy(suffixlist[dtUshort],  ".ushlat");
    strcpy(suffixlist[dtInt],     ".ilat");
    strcpy(suffixlist[dtUint],    ".uilat");
    strcpy(suffixlist[dtLong],    ".llat");
    strcpy(suffixlist[dtUlong],   ".ullat");
    strcpy(suffixlist[dtLlong],   ".lllat");
    strcpy(suffixlist[dtULlong],  ".ulllat");
    strcpy(suffixlist[dtFloat],   ".flat");
    strcpy(suffixlist[dtDouble],  ".dlat");
    strcpy(suffixlist[dtLdouble], ".ldlat");
    strcpy(suffixlist[dtUnknown], ".ulat");

    strcpy(dtypename[dtUchar], "unsigned char");
    strcpy(dtypename[dtSchar], "signed char");
    strcpy(dtypename[dtUshort], "unsigned short");
    strcpy(dtypename[dtShort], "short");
    strcpy(dtypename[dtUint], "unsigned int");
    strcpy(dtypename[dtInt], "int");
    strcpy(dtypename[dtUlong], "unsigned long");
    strcpy(dtypename[dtLong], "long");
    strcpy(dtypename[dtULlong], "unsigned long long");
    strcpy(dtypename[dtLlong], "long long");
    strcpy(dtypename[dtFloat], "float");
    strcpy(dtypename[dtDouble], "double");
    strcpy(dtypename[dtLdouble], "long double");
    strcpy(dtypename[dtUnknown], "unknown");
  }

  //! Get size of a data type
  inline unsigned int SizeOfDataType(DataType v) {
    return datasize[(int)v];
  }

  //! Guess data type by filename extension
  inline DataType GuessDataType(const char *filename) {
    if (!filename)
      return dtUnknown;

    DataType ret = dtUnknown;
    size_t   len = strlen(filename);
    for (int i = 0; i < DT_MAX_DATATYPES; i++) {
      char *suffix = suffixlist[i];
      size_t slen  = strlen(suffix);

      if (len >= slen)
        if (!strcmp(&filename[len-slen], suffix)) {
          ret = (DataType)i;
          break;
        }
    }
    return ret;
  }

  //! Get name corresponding to data type
  inline char *TypeName(DataType d) {
    return dtypename[d];
  }

  //! Get file extension corresponding to data type
  template<class X> const char *GetSuffix() {
    return suffixlist[(int)GetTemplateDataType<X>()];
  }

 protected:
  unsigned int datasize[DT_MAX_DATATYPES];
  char suffixlist[DT_MAX_DATATYPES][8];
  char dtypename[DT_MAX_DATATYPES][64];
};  // class DataTypeInfo

//! Enumeration for compression types of lattices.
typedef enum { ctNone = 0, ctGZIP = 1, ctRLE = 2, ctERLE = 3, ctUnknown = 5 } CompressionType;

//! Names corresponding enumeration for compression types of lattices.
static const char *compTypeName[] = {"no", "GZIP", "RLE", "ERLE", "unknown", 0};
}  // namespace nskaGlobal

#endif  // ifndef KADATATYPES_H
