/*
 * File: kaMask.h
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

#ifndef KAMASK_H
#define KAMASK_H

#include <kaMachineOS.h>
#include <ctype.h>

//! Class for handling of masks for unsigned char data
/*!
   kaMask allows the definition of a 256-array of bools by a text string.
   The string includes numbers and ranges separated by commands, e.g. "1", "3,16", "2-4", "3,15-17,19"
 */

class kaMask {
 public:
  bool m[256];  //!< mask given by 256-array of bools
  kaMask(bool val = false) {
#if KALATTICEDEBUG
    fprintf(stderr, "kaMask::kaMask\n");
#endif  // if KALATTICEDEBUG
    Set(val);
  }

  kaMask(const char *def, bool val = true) {
#if KALATTICEDEBUG
    fprintf(stderr, "kaMask::kaMask\n");
#endif  // if KALATTICEDEBUG
    Set(false);
    Set(def, val);
  }

  void Set(const char *def, bool val = true) {
#if KALATTICEDEBUG
    fprintf(stderr, "kaMask::Set %s\n", def);
#endif  // if KALATTICEDEBUG
    int err = 0;
    if (def)
      while (*def && !err) {
        int v = 0;
        sscanf(def, "%d", &v);
        while (isdigit(*def) || *def == ' ')
          def++;
        int b = v;
        if (*def == '-') {
          def++;
          b = 255;
          sscanf(def, "%d", &b);
          while (isdigit(*def))
            def++;
        }
        if (*def == ',')
          def++;
        else if (*def)
          err = 1;
        if (v < 0)
          v = 0;
        if (v > 255)
          v = 255;
        if (b < 0)
          b = 0;
        if (b > 255)
          b = 255;

        if (!err)
          while (v <= b)
            m[v++] = val;
      }
  }  // Set

  void Set(bool val = false) {
    for (int j = 0; j < 256; j++)
      m[j] = val;
  }
};  // class kaMask

#endif  // ifndef KAMASK_H
