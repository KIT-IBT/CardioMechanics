/**@file FNExtension.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef FNEXTENSION_H
#define FNEXTENSION_H

#include <kaMachineOS.h>

enum FNExtensionLatticeType {
  AS_LAT      =      0,
  AS_SHLAT    =    1,
  AS_ILAT     =     2,
  AS_FLAT     =     3,
  AS_DLAT     =      4,
  AS_TGA      =      5,
  AS_BMP      =      6,
  AS_EPS      =      7,
  AS_RGB      =      8,
  AS_COL      =      9,
  AS_BLD      =    10,
  AS_CCVE     =    11,
  AS_SCVE     =    12,
  AS_ICVE     =    13,
  AS_FCVE     =    14,
  AS_DCVE     =    15,
  AS_GZ       =    16,
  AS_FRE      =    17,
  AS_FRO      =    18,
  AS_RAW      =    19,
  AS_Z        =    20,
  AS_JPEG     =    21,
  AS_JFIF     =    22,
  AS_JPG      =    23,
  AS_TMO      =     24,
  AS_UNKNOWN  =  -1,
  AS_MAX      =    25,
  AS_UNSIGNED = 128,
  AS_ALLOC    =  64
};

class FNExtension {
 public:
  FNExtension(void) {
    strcpy(suffixlist[AS_LAT],   ".lat");
    strcpy(suffixlist[AS_SHLAT], ".shlat");
    strcpy(suffixlist[AS_ILAT],  ".ilat");
    strcpy(suffixlist[AS_FLAT],  ".flat");
    strcpy(suffixlist[AS_DLAT],  ".dlat");
    strcpy(suffixlist[AS_TGA],   ".tga");
    strcpy(suffixlist[AS_BMP],   ".bmp");
    strcpy(suffixlist[AS_EPS],   ".eps");
    strcpy(suffixlist[AS_RGB],   ".rgb");
    strcpy(suffixlist[AS_COL],   ".col");
    strcpy(suffixlist[AS_BLD],   ".bld");
    strcpy(suffixlist[AS_CCVE],  ".ccve");
    strcpy(suffixlist[AS_SCVE],  ".scve");
    strcpy(suffixlist[AS_ICVE],  ".icve");
    strcpy(suffixlist[AS_FCVE],  ".fcve");
    strcpy(suffixlist[AS_DCVE],  ".dcve");
    strcpy(suffixlist[AS_GZ],    ".gz");
    strcpy(suffixlist[AS_FRE],   ".fre");
    strcpy(suffixlist[AS_FRO],   ".fro");
    strcpy(suffixlist[AS_RAW],   ".raw");
    strcpy(suffixlist[AS_Z],     ".Z");
    strcpy(suffixlist[AS_JPEG],  ".jpeg");
    strcpy(suffixlist[AS_JFIF],  ".jfif");
    strcpy(suffixlist[AS_JPG],   ".jpg");
    strcpy(suffixlist[AS_TMO],   ".tmo");
  }

  void AddSuffix(char *out, const char *in, const char *suffix) {
    int len  = strlen(in);
    int slen = strlen(suffix);

    strcpy(out, in);

    if (len >= slen) {
      if (strcmp(&out[len-slen], suffix))
        strcat(out, suffix);
    } else {strcat(out, suffix);}
  }

  void SubSuffix(char *out, const char *in, const char *suffix) {
    int len  = strlen(in);
    int slen = strlen(suffix);

    if (len >= slen) {
      if (!strcmp(&in [len-slen], suffix)) {
        strncpy(out, in, len-slen);
        out[len-slen] = '\0';
      } else {strcpy(out, in);}
    } else {strcpy(out, in);}
  }

  bool CheckExistance(char *file) {
    FILE *test = fopen(file, "rb");
    bool  ret  = false;

    if (test) {
      ret = true;
      fclose(test);
    }

    return ret;
  }

  FNExtensionLatticeType CheckSuffixType(char *name) {
    int ret = AS_UNKNOWN;

    int len = strlen(name);

    for (int i = 0; i < AS_MAX; i++) {
      char *suffix = suffixlist[i];
      int   slen   = strlen(suffix);

      if (!strcmp(&name[len-slen], suffix)) {
        ret = i;
        break;
      }
    }

    return (FNExtensionLatticeType)ret;
  }

  char *GetSuffix(unsigned int i) {
    if (i < AS_MAX)
      return suffixlist[i];

    return NULL;
  }

 private:
  char suffixlist[AS_ALLOC][64];
};  // class FNExtension


#endif  // ifndef FNEXTENSION_H
