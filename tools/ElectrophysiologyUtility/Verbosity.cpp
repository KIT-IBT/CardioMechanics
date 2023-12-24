/*
 * File: Verbosity.cpp
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


#include "Verbosity.h"

/* initialize static variables */
unsigned int Verbosity::verbose = 0;
Verbosity::fakestream Verbosity::fake;

ostream & Verbosity::Stream(unsigned int level) {
  if (verbose >= level)
    return std::cerr;
  else
    return fake;
}
