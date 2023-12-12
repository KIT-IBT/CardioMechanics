/*
 * File: kaBasicIO_intel.h
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


#ifdef KABASICIO_H
#undef KABASICIO_H

#ifdef KA_LITTLE_ENDIAN
# undef KA_LITTLE_ENDIAN
#else  // ifdef KA_LITTLE_ENDIAN
# define KA_LITTLE_ENDIAN
#endif  // ifdef KA_LITTLE_ENDIAN

#define kaWrite kaWrite_intel
#define kaRead kaRead_intel
#define Swap Swap_intel
#define Raw Raw_intel
#define GZIP GZIP_intel
#define kaRLE kaRLE_intel
#define ibufsize ibufsize_intel
#define obufsize obufsize_intel

#include <kaBasicIO.h>

#undef kaWrite
#undef kaRead
#undef Swap
#undef Raw
#undef GZIP
#undef kaRLE
#undef ibufsize
#undef obufsize

#ifdef KA_LITTLE_ENDIAN
# undef KA_LITTLE_ENDIAN
#else  // ifdef KA_LITTLE_ENDIAN
# define KA_LITTLE_ENDIAN
#endif  // ifdef KA_LITTLE_ENDIAN


#endif  // ifdef KABASICIO_H
