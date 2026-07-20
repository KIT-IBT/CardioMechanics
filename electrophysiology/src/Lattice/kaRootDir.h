/*
 * File: kaRootDir.h
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


#ifndef KAROOTDIR_H
#define KAROOTDIR_H

#include <kaMachineOS.h>

//! Resolves the location of the bundled model parameter files (electrophysiology/data).
/*!
   Defaults to the source directory baked in at build time (CM_SOURCE_DIR); the
   "kaRootDir" environment variable overrides it when set.
 */

class kaRootDir {
  std::string data;

 public:
  kaRootDir(const char *name) {
    const char *p = getenv("kaRootDir");
    std::string base = (p && *p) ? p : CM_SOURCE_DIR;
    data = base + "/electrophysiology/data/" + name;
  }

  const char *GetData() {return data.c_str();}
};  // class kaRootDir

#endif  // ifndef KAROOTDIR_H
