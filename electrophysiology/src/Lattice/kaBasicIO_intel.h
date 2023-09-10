/*
 *  kaBasicIO_intel.h
 *
 *
 *  Created by Dimitri Farina on 16.05.2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 *  file defines kaWrite_intel and kaRead_intel templates, which enable to write and read
 *  little-endian data on any platform
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
