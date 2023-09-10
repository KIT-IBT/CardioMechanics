/*! \file LatticeCreate.cpp
   \brief Functions for program LatticeCreate
   \date created ??/??/?? ??\n
        man     02/26/03 mm
   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
   \sa Synopsis \ref LatticeCreate
 */

// doxygen manpage LatticeCreate
/*! \page LatticeCreate LatticeCreate
    creates a kaLattice data structure

   \section SYNOPSIS_LatticeCreate SYNOPSIS
   LatticeCreate \<kalattice latticename.suffix\> [-dim \<w\> \<h\> \<d\>] [-fill \<value\>] [-rle | -gzip] [-version
 \<value\>]

   \section OPTIONS_LatticeCreate OPTIONS
   \param  "<kaLattice latticename.suffix>"  name of lattice (for suffix see \ref DESCRIPTION_LatticeCreate )
   \param  "-dim <w> <h> <d>" specify the lattice dimensions: width, height, depth
   \param  "-fill <value>" set all lattice entries to \<value\> (Default value 0)
   \param  "-rle | -gzip" set lattice compression to "rle" or "gzip"
   \param  "-version <value>" set version info

   \section DESCRIPTION_LatticeCreate DESCRIPTION
   LatticeCreate creates a lattice data structure of name "latticename".
   The suffix specifies the type of data to be stored.
   \par suffixes (see DataTypeInfo)

   \li \b .lat   unsigned char
   \li \b .slat   signed char
   \li \b .shlat  short
   \li \b .ushlat unsigned short
   \li \b .ilat   int
   \li \b .uilat  unsigned int
   \li \b .llat   long
   \li \b .ullat  unsigned long
   \li \b .lllat  long long
   \li \b .ulllat unsigned long long
   \li \b .flat   float
   \li \b .dlat   double
   \li \b .ldlat  long double
   \li \b .ulat   unknown

   The parameter -dim specifies the lattice dimensions. Width, heigth and depth can be set. (Default values are -dim 1 1
      1.)
   The lattice can be filled with -fill \<value\> option. The default value is 0.

   A Compression may be set using -rle or -gzip for compression algorithms.

   The version information can be set. (??)

   \section SOURCE_LatticeCreate SOURCE
   LatticeCreate.cpp

   \section SEEALSO_LatticeCreate SEE ALSO
   LatticeInfo LatticeHeaderInfo LatticeSet LatticeSize
 *
 */

#include <kaLattice.h>
#include <kaExceptions.h>
#include <kaParser.h>

//! Template function for creating lattice

/*!
   The arguments include the definition of lattice and an initial fill value
 */
template<class X> void LatticeCreate(kaLatticeCreate &cls, X fill) {
  kaLattice<X> T(cls);
  T.FillWith(fill);
  T.Save();
}

//! Function for parsing of command line arguments and call of function LatticeCreate
void LatticeCreateSub(int argc, char **argv) {
  nskaGlobal::FileCheck FC;

  FC.Check(argv[1]);
  if (FC.Exist())
    throw kaBaseException("Lattice %s already exists", argv[1]);

  kaLatticeCreate cls;
  strcpy(cls.Name, argv[1]);

  cls.xLattice = 1;
  cls.yLattice = 1;
  cls.zLattice = 1;

  kaParser CL(argc, argv);

  if (CL.IsOption("-dim")) {
    IndexType w = 0, h = 0, d = 0;
    if (sscanf(CL.GetOption("-dim", 3), "%d %d %d", &w, &h, &d) == 3) {
      cls.xLattice = w;
      cls.yLattice = h;
      cls.zLattice = d;
    }
  }

  cls.compression = ctNone;
  if (CL.IsOption("-rle"))
    cls.compression = ctRLE;
  if (CL.IsOption("-gzip"))
    cls.compression = ctGZIP;

  cls.writeVersion = ioNow;
  if (CL.IsOption("-version"))
    cls.writeVersion = (kaIOVersion)atoi(CL.GetOption("-version", 1));

  DataTypeInfo dti;

  switch (dti.GuessDataType(argv[1])) {
    case dtUchar: {
      uint8_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (uint8_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<uint8_t>(cls, fill);
    }
    break;
    case dtSchar: {
      int8_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (int8_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<int8_t>(cls, fill);
    }
    break;
    case dtShort: {
      int16_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (int16_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<int16_t>(cls, fill);
    }
    break;
    case dtUshort: {
      uint16_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (uint16_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<uint16_t>(cls, fill);
    }
    break;
    case dtInt: {
      int32_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (int32_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<int32_t>(cls, fill);
    }
    break;
    case dtUint: {
      uint32_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (uint32_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<uint32_t>(cls, fill);
    }
    break;
    case dtLong: {
      int64_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (int64_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<int64_t>(cls, fill);
    }
    break;
    case dtUlong: {
      uint64_t fill = 0;
      if (CL.IsOption("-fill"))
        fill = (uint64_t)atoi(CL.GetOption("-fill", 1));
      LatticeCreate<uint64_t>(cls, fill);
    }
    break;
    case dtFloat: {
      float fill = 0;
      if (CL.IsOption("-fill"))
        fill = (float)atof(CL.GetOption("-fill", 1));
      LatticeCreate<float>(cls, fill);
    }
    break;
    case dtDouble: {
      double fill = 0;
      if (CL.IsOption("-fill"))
        fill = (double)atof(CL.GetOption("-fill", 1));
      LatticeCreate<double>(cls, fill);
    }
    break;
    default:
      throw kaBaseException("Unknown or unsupported type of lattice %s", argv[1]);
  }  // switch
}  // LatticeCreateSub

//! Main function

/*! The function handles the printing of a help message, the call of function LatticeCreateSub and provides error
   handler via exception catching
 */
int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << argv[0] << " <kaLattice latticename.suffix>" << endl;
    cerr << "\t[-dim <w> <h> <d>]" << endl;
    cerr << "\t[-version <value>]" << endl;
    cerr << "\t[-rle | -gzip]" << endl;
    cerr << "\t[-fill <value>]" << endl;
    exit(-1);
  }

  try {
    LatticeCreateSub(argc, argv);
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}
