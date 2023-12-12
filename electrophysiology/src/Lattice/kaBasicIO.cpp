/*
 * File: kaBasicIO.cpp
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


#include <kaBasicIO.h>

namespace nskaGlobal {
#ifdef KA_LITTLE_ENDIAN

void Swap(unsigned char *data, unsigned int size) {}

void Swap(signed char *data, unsigned int size) {}

void Swap(char *data, unsigned int size) {}

#endif  // ifdef KA_LITTLE_ENDIAN
}
