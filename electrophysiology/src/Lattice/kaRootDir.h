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

//! Class for handling of default root directories
/*!
   The default is hardcoded, but can be overwritten by setting the system variable "kaRootDir"

 */

class kaRootDir {
  std::string root;
  std::string data;
  std::string bin;

 public:
  kaRootDir(const char *name) {
    char *p = getenv("kaRootDir");

    root = p + std::string("/electrophysiology/") + name;
    data = p + std::string("/electrophysiology/data/") + name;

    # ifdef osLinux
    bin = p + std::string("/bin/linux/") + name;
    # endif  // ifdef osLinux
    # ifdef osMac
    bin = p + std::string("/bin/macosx/") + name;
    # endif  // ifdef osMac

  }

  const char *GetRoot() {return root.c_str();}
  const char *GetData() {return data.c_str();}
  const char *GetBin() {return bin.c_str();}
};  // class kaRootDir

#endif  // ifndef KAROOTDIR_H
