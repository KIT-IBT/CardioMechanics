/*! \file kaBasicIO.cpp
 \brief Basic functions for lattice-file input-output
 
 
 \author gs IBT - Karlsruhe Institute of Technology
 */

#include <kaBasicIO.h>

namespace nskaGlobal {
#ifdef KA_LITTLE_ENDIAN

void Swap(unsigned char *data, unsigned int size) {}

void Swap(signed char *data, unsigned int size) {}

void Swap(char *data, unsigned int size) {}

#endif  // ifdef KA_LITTLE_ENDIAN
}
