/*
 * File: Verbosity.h
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


#ifndef VERBOSITY_H
#define VERBOSITY_H

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

class Verbosity {
 public:
  class fakestream : public ostream {
   public:
    fakestream() : ios(0), ostream(0) {}
  };

  static ostream & Stream(unsigned int);

  static void Set(unsigned int level) {verbose = level;}

  static unsigned int Get() {return verbose;}

 private:
  static unsigned int verbose;
  static fakestream fake;
};

#endif  // ifndef VERBOSITY_H
