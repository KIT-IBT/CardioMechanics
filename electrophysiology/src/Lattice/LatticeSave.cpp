/*
 * File: LatticeSave.cpp
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


#include <kaLattice.h>

#include <kaParser.h>

//! Template function for saving of lattice
template<class X> void LatticeSave(kaLatticeCreate &cls, char *name) {
  kaLattice<X> T(cls);
  T.SetWriteVersion(cls.writeVersion);
  T.SetCompression(cls.compression);
  T.Save(name);
}

//! Main function

/*! The function handles the printing of a help message, the parsing of command line arguments, the call of template
   function LatticeSave and provides error handler via exception catching
 */
int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << argv[0] << " <lattice>" << endl;
    cerr << "\t[-lname <lattice>]" << endl;
    cerr << "\t[-raw][-rle][-gzip]" << endl;
    cerr << "\t[-version <1|2|3>]" << endl;
    return -1;
  }

  try {
    kaParser CL(argc, argv);

    kaLatticeCreate cls;
    strcpy(cls.Name, argv[1]);
    if (CL.IsOption("-raw"))
      cls.compression = ctNone;
    if (CL.IsOption("-rle"))
      cls.compression = ctRLE;
    if (CL.IsOption("-gzip"))
      cls.compression = ctGZIP;

    char *name = NULL;
    if (CL.IsOption("-lname"))
      name = CL.GetOption("-lname", 1);

    if (CL.IsOption("-version"))
      cls.writeVersion = (kaIOVersion)atoi(CL.GetOption("-version", 1));

    kaLatticeSharedHeader sh(cls);
    DataType dt = sh.lhDtype();

    switch (dt) {
      case dtUchar:  {
        LatticeSave<uint8_t>(cls, name);
      }
      break;
      case dtSchar:  {
        LatticeSave<int8_t>(cls, name);
      }
      break;
      case dtShort:  {
        LatticeSave<int16_t>(cls, name);
      }
      break;
      case dtUshort: {
        LatticeSave<uint16_t>(cls, name);
      }
      break;
      case dtInt:    {
        LatticeSave<int32_t>(cls, name);
      }
      break;
      case dtUint:   {
        LatticeSave<uint32_t>(cls, name);
      }
      break;
      case dtLong:   {
        LatticeSave<int64_t>(cls, name);
      }
      break;
      case dtUlong:  {
        LatticeSave<uint64_t>(cls, name);
      }
      break;
      case dtFloat:  {
        LatticeSave<float>(cls, name);
      }
      break;
      case dtDouble: {
        LatticeSave<double>(cls, name);
      }
      break;
      default:
        throw kaBaseException("Unknown or unsupported type of lattice %s", argv[1]);
    }  // switch
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
