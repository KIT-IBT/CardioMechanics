/*! \file kaLatticeHeader.h
   \brief Classes for handling of lattice headers.

   Lattice headers include information concerning size, type, geometry of lattices.
   Public interfaces are the classes kaLatticeHeader and kaLatticeSharedHeader

   \author cw,fs, IBT - Universität Karlsruhe (TH)

   19.12.2003 df: IOVersion -> kaIOVersion
 */
#ifndef KALATTICEHEADER_H
#define KALATTICEHEADER_H

#include <kaBasicIO.h>
#include <kaDataTypes.h>
#include <kaIPC.h>
#include <kaSharedMatrixN.h>

using namespace nskaGlobal;

typedef int32_t IndexType;

//! this is the magic id keyword for a file of type kaLattice

#define MAGIC_ID "LATTICE*IBT*KA"  // gs 16.02.04
// static const char * MAGIC_ID = "LATTICE*IBT*KA";
static const int32_t MAGIC_ID_LEN = 14;

static const int32_t SMALLEST_HEADER = 148;
static const int32_t LH_PARAMS       = 12;

static const int32_t MEM_HEADER = 256;  // is aligned to voxel elements of size 1, 2, 4, 8, 16, 32, 64, 128, 256

//! Enumeration for the existing different file versions.

typedef enum { ioUnknown = 0, ioInitial = 1, ioLast = 2, ioNow = 3 } kaIOVersion;


//! Class to build objects of class kaLattice

class kaLatticeCreate {
 public:
  //! Default constructor
  kaLatticeCreate(void) {
    Name[0]      = 0;
    xLattice     = -1;
    yLattice     = -1;
    zLattice     = -1;
    compression  = ctNone;
    writeVersion = ioNow;
  }

  char Name[256];              //!< Name of lattice.
  IndexType xLattice;          //!< Default size of lattice in x-direction. Used for creation of lattice.
  IndexType yLattice;          //!< Default size of lattice in y-direction. Used for creation of lattice.
  IndexType zLattice;          //!< Default size of lattice in z-direction. Used for creation of lattice.
  CompressionType compression;  //!< Default compression type of lattice. Used for creation of lattice.
  kaIOVersion writeVersion;    //!< Default type version of lattice. Used for creation of lattice.
};


/* *******************************************************************************************************/
/*                                                                                                       */
/* description of the lattice file header:                                                               */
/*                                                                                                       */
/* the file can be devided into two parts, the header and the data part:                                 */
/* the header contains information about                                                                 */
/* - the lattice magic id keyword                                                                        */
/* - the lattice file version                                                                            */
/* - the memory size, allocated as shared memory                                                         */
/* - the grid size,                                                                                      */
/* - the transformation matrix                                                                           */
/* - the type of the data                                                                                */
/* - the compression method                                                                              */
/* - the size of one voxel in bytes                                                                      */
/*   and                                                                                                 */
/* - the number of voxels.                                                                               */
/*                                                                                                       */
/* the exact assignment of the bytes is desribed below:                                                  */
/*                                                                                                       */
/* ------------------------------------------------------------------------------------------------------*/
/* data type       | variable name    | number of bytes needed      | description                        */
/* ------------------------------------------------------------------------------------------------------*/
/* char *          | MAGIC_ID         | 14 * sizeof(char)      = 14 | magic id keyword                   */
/* kaIOVersion     | lhReadVersion()  | sizeof(kaIOVersion)    =  4 | file version                       */
/* kaIOVersion     | lhWriteVersion() | sizeof(kaIOVersion)    =  4 | file version                       */
/* uint32_t        | SharedMemorySize | sizeof(uint32_t)       =  4 | size of shared memory block        */
/* uint32_t        | xLattice         | sizeof(uint32_t)       =  4 | width of data block                */
/* uint32_t        | yLattice         | sizeof(uint32_t)       =  4 | height of data block               */
/* uint32_t        | zLattice         | sizeof(uint32_t)       =  4 | depth of data block                */
/* uint32_t        | m.Dimension()    | sizeof(uint32_t)       =  4 | dimension of transformation matrix */
/* double *        | m.a(row,column)  | 16 * sizeof(double)    =128 | elements of matrix                 */
/* DataType        | lhDtype()        | sizeof(DataType)       =  4 | type of data                       */
/* CompressionType | lhCompression()  | sizeof(CompressionType)=  4 | compression method                 */
/* uint32_t        | lhElemSize()     | sizeof(uint32_t)   =      4 | size of one element                */
/* char *          | no name          |                          88 | alignment bytes (custom use)       */
/* uint32_t        | lhDsize()        | sizeof(uint32_t)   =      4 | number of voxels                   */
/* ------------------------------------------------------------------------------------------------------*/
/*                                                                                                       */
/*                                            total header size=274                                      */
/*                                                                                                       */
/* the data part contains the data in sequential order,                                                  */
/* beginning at the left, front, top corner of the data block.                                           */
/*                                                                                                       */
/* *******************************************************************************************************/


//! Class provides methods to examine a lattice header on disk and to store this information in private members

class kaLatticeHeader : virtual public nskaGlobal::DataTypeInfo, virtual public nskaIPC::IPCKey {
 public:
  //! Constructur for loading and creation of lattice header in shared memory
  kaLatticeHeader(kaLatticeCreate &cls) : nskaIPC::IPCKey(cls.Name, cls.xLattice != -1) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::kaLatticeHeader(kaLatticeCreate & cls)" << endl;
#endif  // if KALATTICEDEBUG
    PreInit();
    New();
  }

  //! Destructor
  virtual            ~kaLatticeHeader(void) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::~kaLatticeHeader(void)" << endl;
#endif  // if KALATTICEDEBUG
    Delete();
  }

  //! Initialization of components of lattice header
  void Init(char *pmem) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Init(char * pmem)" << endl;
#endif  // if KALATTICEDEBUG
    lh_data = pmem;
    char *pmatrix = lh_data + lh_memoffset[5];
    lh_matrix.Detach();
    lh_matrix.Attach(pmatrix, 4);
  }

  //! Initialization of components of lattice header
  void Define(kaLatticeCreate &cls, int32_t elementSize) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Define(kaLatticeCreate &, int32_t)" << endl;
#endif  // if KALATTICEDEBUG
    lhWriteVersion(ioNow);
    lhWidth(cls.xLattice);
    lhHeight(cls.yLattice);
    lhDepth(cls.zLattice);
    lhDtype(dtUnknown);
    lhCompression(ctNone);
    lhElemSize(elementSize);
    lhSMSize(calcSMSize());
    lhLoadFlag(false);
    lhDsize(calcDataSize());
    lhMatrix().Identity();
  }

  //! Uninitialization of components of lattice header
  void Uninit(void) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Uninit(void)" << endl;
#endif  // if KALATTICEDEBUG
    lh_matrix.Detach();
    lh_data = 0;
  }

  //! Memory and parameter initialization of components of lattice header
  void New() {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::New()" << endl;
#endif  // if KALATTICEDEBUG

    lh_data = (char *)lh_fallbackData;
    memset(lh_data, 0, sizeof(char) * MEM_HEADER);

    char *pmatrix = lh_data + lh_memoffset[5];
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::New() --> lh_data = "  << (long)lh_data << endl;
    cerr << "kaLatticeHeader::New() --> pmatrix = " << (long)pmatrix << endl;
#endif  // if KALATTICEDEBUG
    lh_matrix.Attach(pmatrix, 4);
    lh_matrix.Identity();
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::New() --> lh_data = "  << (long)lh_data << endl;
    cerr << "kaLatticeHeader::New() --> pmatrix = " << (long)pmatrix << endl;
#endif  // if KALATTICEDEBUG

    lhWriteVersion(ioUnknown);
    lhSMSize(0);
    lhWidth(0);
    lhHeight(0);
    lhDepth(0);
    lhDtype(dtUnknown);
    lhCompression(ctNone);
    lhElemSize(0);
    lhLoadFlag(false);
    lhDsize(0);
  }  // New

  //! Uninitialization of components of lattice header
  inline void Delete(void) {Uninit();}

  //! Read header from file
  inline bool headerRead(void) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::headerRead(void)" << endl;
#endif  // if KALATTICEDEBUG
    bool ret  = false;
    bool rmFD = false;

    if (!isReadFD()) {
      createReadFD();
      rmFD = true;
    }
    if (!isReadFD())
      return false;

    ret = Read();

    if (rmFD)
      closeReadFD();
    return ret;
  }

  //! Get lattice file version
  inline kaIOVersion lhReadVersion(void) const {return lh_readversion;}

  //! Get header size
  inline int32_t lhHeaderSize(void) const {return MEM_HEADER;}

  //! Get matrix describing geometry used e.g. in class kaLattice
  inline kaSharedMatrixN<double> & lhMatrix(void) {return lh_matrix;}

  //! Get lattice file version used for output
  inline kaIOVersion lhWriteVersion(void) const {return *((kaIOVersion *)(lh_data + lh_memoffset[0]));}

  //! Get size of shared memory block
  inline uint32_t lhSMSize(void) const {return *((uint32_t *)(lh_data + lh_memoffset[1]));}

  //! Get width used e.g. in class kaLattice
  inline uint32_t lhWidth(void) const {return *((uint32_t *)(lh_data + lh_memoffset[2]));}

  //! Get height used e.g. in class kaLattice
  inline uint32_t lhHeight(void) const {return *((uint32_t *)(lh_data + lh_memoffset[3]));}

  //! Get depth used e.g. in class kaLattice
  inline uint32_t lhDepth(void) const {return *((uint32_t *)(lh_data + lh_memoffset[4]));}

  //! Get dimension of matrix
  inline uint32_t lhMatrixDim(void) const {return *((uint32_t *)(lh_data + lh_memoffset[5]));}

  //! Get data type
  inline DataType lhDtype(void) const {return *((DataType *)(lh_data + lh_memoffset[6]));}

  //! Get compression type
  inline CompressionType lhCompression(void) const {return *(CompressionType *)(lh_data + lh_memoffset[7]);}

  //! Get size of an element
  inline uint32_t lhElemSize(void) const {return *((uint32_t *)(lh_data + lh_memoffset[8]));}

  //! Get size of an element
  inline uint8_t lhLoadFlag(void) const {return *((uint8_t *)(lh_data + lh_memoffset[9]));}

  inline uint8_t *lhAlignment(void) const {return (uint8_t *)(lh_data + lh_memoffset[10]);}

  inline uint32_t lhDsize(void) const {return *((uint32_t *)(lh_data + lh_memoffset[11]));}

  inline uint8_t *lhHeaderMem(void) const {return (uint8_t *)lh_data;}

  inline void lhWriteVersion(kaIOVersion v) {*((kaIOVersion *)(lh_data + lh_memoffset[0])) = v;}

  inline void lhSMSize(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[1])) = v;}

  inline void lhWidth(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[2])) = v;}

  inline void lhHeight(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[3])) = v;}

  inline void lhDepth(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[4])) = v;}

  inline void lhMatrix(double *v) {lh_matrix.SetFromArray(v);}

  inline void lhMatrix(kaSharedMatrixN<double> &v) {lh_matrix = v;}

  inline void lhDtype(DataType v) {*((DataType *)(lh_data + lh_memoffset[6])) = v;}

  inline void lhCompression(CompressionType v) {*((CompressionType *)(lh_data + lh_memoffset[7])) = v;}

  inline void lhElemSize(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[8])) = v;}

  inline void lhAlignment(uint8_t *palign) {memcpy(((uint8_t *)(lh_data + lh_memoffset[10])), palign, lh_alignment);}

  inline void lhDsize(uint32_t v) {*((uint32_t *)(lh_data + lh_memoffset[11])) = v;}

  inline void lhHeaderMem(uint8_t *pmem) {memcpy(lh_data, pmem, MEM_HEADER);}

  inline void lhDim(kaLatticeCreate &cls, int32_t esize = 1) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::lhDim(kaLatticeCreate & cls, int32_t esize = 1)" << endl;
#endif  // if KALATTICEDEBUG
    lhWidth(cls.xLattice);
    lhHeight(cls.yLattice);
    lhDepth(cls.zLattice);
    lhElemSize(esize);
  }

  //! Print file header in file
  void FileHeaderOut(FILE *out) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::FileHeaderOut(FILE * out)" << endl;
#endif  // if KALATTICEDEBUG
    fprintf(out, "name:               %s\n", getName());
    fprintf(out, "created:            %s\n", isCreated() ? "true" : "false");
    char help[256];
    switch (lhReadVersion()) {
      case ioUnknown:
        strcpy(help, "ioUnknown (0)"); break;
      case ioInitial:
        strcpy(help, "ioInitial (1)"); break;
      case ioLast:
        strcpy(help, "ioLast (2)");    break;
      case ioNow:
        strcpy(help, "ioNow (3)");     break;
      default:
        strcpy(help, "no valid information"); break;
    }
    fprintf(out, "read version:       %s\n", help);
    switch (lhWriteVersion()) {
      case ioUnknown:
        strcpy(help, "ioUnknown (0)"); break;
      case ioInitial:
        strcpy(help, "ioInitial (1)"); break;
      case ioLast:
        strcpy(help, "ioLast (2)");    break;
      case ioNow:
        strcpy(help, "ioNow (3)");     break;
      default:
        strcpy(help, "no valid information"); break;
    }
    fprintf(out, "write version:      %s\n", help);
    fprintf(out, "shared memory size: %u\n", lhSMSize());
    fprintf(out, "width:              %d\n", lhWidth());
    fprintf(out, "height:             %d\n", lhHeight());
    fprintf(out, "depth:              %d\n", lhDepth());
    fprintf(out, "transform matrix:\n");
    for (int i = 0; i < lh_matrix.Dimension(); i++) {
      fprintf(out, "| ");
      for (int j = 0; j < lh_matrix.Dimension(); j++)
        fprintf(out, "%9.6lf ", lh_matrix.a(i, j));
      fprintf(out, "|\n");
    }
    fprintf(out, "data type:          %s\n", TypeName(lhDtype()));
    fprintf(out, "compression:        %s\n", compTypeName[lhCompression()]);
    fprintf(out, "element size:       %d\n", lhElemSize());
    fprintf(out, "load flag:          %s\n", lhLoadFlag() ? "true" : "false");
    fprintf(out, "data size (voxels): %u\n", lhDsize());
  }  // FileHeaderOut

  inline int32_t calcSMSize(void) const {return lhHeaderSize() + calcDataSize();}

  inline int32_t calcDataSize(void) const {return calcNumVoxels() * lhElemSize();}

  inline int32_t calcNumVoxels(void) const {return lhWidth() * lhHeight() * lhDepth();}

  inline int32_t calcSizeOfFile(void) {
    if (headerRead())
      return calcSMSize();

    //     else if(Exist()) throw kaBaseException("Cannot calculate shared memory size from lattice file");
    return 0;
  }

  //! Compare header dimensions.
  inline void DimMatch(kaLatticeHeader *LH) {
    if (!LH)
      throw DimensionMismatch(getName());

    if ((lhWidth() != LH->lhWidth()) || (lhHeight() != LH->lhHeight()) || (lhDepth() != LH->lhDepth()) )
      throw DimensionMismatch(getName(), LH->getName());
  }

  //! Compare header dimensions.
  inline void DimMatch(IndexType x, IndexType y, IndexType z) {
    if ((lhWidth() != x) || (lhHeight() != y) || (lhDepth() != z) )
      throw DimensionMismatch(getName(), x, y, z);
  }

  //! Compare header dimensions.
  inline void DimMatch(kaLatticeCreate &cls) {
    if ((lhWidth() != cls.xLattice) || (lhHeight() != cls.yLattice) || (lhDepth() != cls.zLattice) )
      throw DimensionMismatch(getName(), cls.xLattice, cls.yLattice, cls.zLattice);
  }

  //! Compare header datatype.
  inline void TypeMatch(DataType d) {
    if (lhDtype() != d)
      throw TypeMismatch(getName());
  }

  //! Compare header datatype.
  inline void TypeMatch(kaLatticeHeader *LH) {
    if (lhDtype() != LH->lhDtype())
      throw TypeMismatch(getName(), LH->getName());
  }

 protected:
  void lhLoadFlag(bool v) {*((uint8_t *)(lh_data + lh_memoffset[9])) = v;}

  bool Read(void) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Read()" << endl;
#endif  // if KALATTICEDEBUG
    if (!FileCheck::Exist()) {
      // file does not exist
#if KALATTICEDEBUG
      cerr << "kaLatticeHeader::Read(): file does not exist." << endl;
#endif  // if KALATTICEDEBUG
      return false;
    }
    if (!FileCheck::Access()) {
#if KALATTICEDEBUG
      cerr << "kaLatticeHeader::Read(): file cannot be accessed." << endl;
#endif  // if KALATTICEDEBUG
      return false;
    }

    if (FileCheck::Size() < SMALLEST_HEADER) {
      // header incomplete
#if KALATTICEDEBUG
      cerr << "kaLatticeHeader::Read(): header is incomplete (too small)." << endl;
#endif  // if KALATTICEDEBUG
      return false;
    }

    // the file descriptor
    int fd = getReadFD();
    lseek(fd, 0, SEEK_SET);

    // check for magic id keyword
    char idcheck[MAGIC_ID_LEN];
    kaRead<char>(idcheck, MAGIC_ID_LEN, fd);

    if (!strncmp(idcheck, MAGIC_ID, MAGIC_ID_LEN)) {
      // new versions

      uint32_t nversion = (short)ioNow;
      kaRead<uint32_t>(&nversion, fd);
      if (nversion > (uint32_t)ioNow)
        fprintf(stderr, "Warning: the version information says this is a future version (%d). trying version %d ...\n",
                nversion, (uint32_t)ioNow);
      lhReadVersion((kaIOVersion)nversion);
      switch (lhReadVersion()) {
        case ioNow:
        default: {
          uint32_t ulDummy;
          double   mhelp[16];
          uint8_t *cDummy = new uint8_t[lh_alignment];
          uint8_t  bDummy;

          kaRead<uint32_t>(&ulDummy, fd); lhWriteVersion((kaIOVersion)ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhSMSize(ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhWidth(ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhHeight(ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhDepth(ulDummy);
          kaRead<uint32_t>(&ulDummy, fd);  // the matrix size is fixed to 4
          kaRead<double>(mhelp, 16, fd); lh_matrix.SetFromArray(mhelp);
          kaRead<uint32_t>(&ulDummy, fd); lhDtype((DataType)ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhCompression((CompressionType)ulDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhElemSize(ulDummy);
          kaRead<uint8_t>(&bDummy, fd);
          kaRead<uint8_t>(cDummy, lh_alignment, fd); lhAlignment(cDummy);
          kaRead<uint32_t>(&ulDummy, fd); lhDsize(ulDummy);
          delete[] cDummy;
        }
        break;
      }
    } else {
      // old versions

      uint32_t ulDummy;
      lhReadVersion(ioLast);
      lhCompression(ctNone);
      lhWriteVersion(ioNow);
      lseek(fd, 0, SEEK_SET);
      kaRead<uint32_t>(&ulDummy, fd); lhWidth(ulDummy);
      kaRead<uint32_t>(&ulDummy, fd); lhHeight(ulDummy);
      kaRead<uint32_t>(&ulDummy, fd); lhDepth(ulDummy);
      kaRead<uint32_t>(&ulDummy, fd);  // the matrix size is fixed to 4
      double mhelp[16];
      kaRead<double>(mhelp, 16, fd);
      lh_matrix.SetFromArray(mhelp);
      long xyz     = lseek(fd, 0, SEEK_CUR);
      long fileend = lseek(fd, 0, SEEK_END);
      long rest    = fileend - xyz - sizeof(uint32_t);
      lseek(fd, xyz, SEEK_SET);
      uint32_t volsize = lhWidth() * lhHeight() * lhDepth();
      uint32_t hdsize  = 0;
      kaRead<uint32_t>(&hdsize, fd);
      if (hdsize != volsize) {kaRead<uint32_t>(&hdsize, fd); rest -= sizeof(uint32_t);}
      if (hdsize != volsize) {
        fprintf(stderr, "failure reading header (%d (hdsize) != %d (volsize))\n", hdsize, volsize);
        lhDtype(dtUchar);
        lhDsize(volsize);
        lhElemSize(SizeOfDataType(lhDtype()));
      } else {
        lhElemSize((uint32_t)rest / volsize);

        // fprintf(stderr, "lh_elemsize = rest / volsize = %ld / %ld = %ld\n", rest, volsize, lh_elemsize);
        lhDtype(GuessDataType(getName()));
        lhDsize(hdsize);
      }
      lhSMSize(calcSMSize());
    }
#if KALATTICEDEBUG == 2
    FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 2
    return true;
  }  // Read

  bool Write(int fd) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Write(int fd)" << endl;
#endif  // if KALATTICEDEBUG
#if KALATTICEDEBUG == 3
    FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 3

    switch (lhWriteVersion()) {
      case ioInitial:
        return false;

        break;
      case ioLast: {
        kaWrite<uint32_t>(lhWidth(), fd);
        kaWrite<uint32_t>(lhHeight(), fd);
        kaWrite<uint32_t>(lhDepth(), fd);
        kaWrite<uint32_t>(4, fd);  // the value lh_matrix.Dimension() value  must be always 4
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            kaWrite<double>(lh_matrix.a(i, j), fd);
        uint32_t headbytes  = 3 * sizeof(uint32_t) + sizeof(uint32_t) + 4 * 4 * sizeof(double);
        unsigned char dummy = 0;
        while ((headbytes - 0 - sizeof(uint32_t)) % lhElemSize()) {headbytes++; kaWrite<unsigned char>(dummy, fd);}
        kaWrite<uint32_t>(lhDsize(), fd);
      }
      break;
      case ioNow: {
        kaWrite<char>((char *)MAGIC_ID, MAGIC_ID_LEN, fd);
        kaWrite<uint32_t>((uint32_t)ioNow, fd);
        kaWrite<uint32_t>((uint32_t)lhWriteVersion(), fd);
        kaWrite<uint32_t>(lhSMSize(), fd);
        kaWrite<uint32_t>(lhWidth(), fd);
        kaWrite<uint32_t>(lhHeight(), fd);
        kaWrite<uint32_t>(lhDepth(), fd);
        kaWrite<uint32_t>(4, fd);  // the value lh_matrix.Dimension() value must be always 4
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            kaWrite<double>(lh_matrix.a(i, j), fd);
        kaWrite<uint32_t>((uint32_t)lhDtype(), fd);
        kaWrite<uint32_t>((uint32_t)lhCompression(), fd);
        kaWrite<uint32_t>((uint32_t)lhElemSize(), fd);
        kaWrite<uint8_t>(0, fd);
        kaWrite<uint8_t>(lhAlignment(), lh_alignment, fd);
        kaWrite<uint32_t>(lhDsize(), fd);
      }
      break;
      case ioUnknown:
      default:
        return false;

        break;
    }  // switch
    return true;
  }  // Write

  uint32_t lh_memoffset[LH_PARAMS];
  uint32_t lh_bytes[LH_PARAMS];

 private:
  void PreInit() {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::PreInit(void)" << endl;
#endif  // if KALATTICEDEBUG
    lh_data = (char *)lh_fallbackData;

    Offsets();

    lhReadVersion(ioUnknown);
    lhLoadFlag(false);
    lhWriteVersion(ioNow);
    lhCompression(ctRLE);
    lhDsize(0);
  }

  void Offsets(void) {
#if KALATTICEDEBUG
    cerr << "kaLatticeHeader::Offsets(void)" << endl;
#endif  // if KALATTICEDEBUG
    lh_bytes[0] = sizeof(uint32_t);  // write version
    lh_bytes[1] = sizeof(uint32_t);  // shared memory size
    lh_bytes[2] = sizeof(uint32_t);  // width
    lh_bytes[3] = sizeof(uint32_t);  // height
    lh_bytes[4] = sizeof(uint32_t);  // depth
    lh_bytes[5] = sizeof(uint32_t)
      +16*sizeof(double);            // matrix
    lh_bytes[6]  = sizeof(uint32_t); // data type
    lh_bytes[7]  = sizeof(uint32_t); // compression
    lh_bytes[8]  = sizeof(uint32_t); // element size
    lh_bytes[9]  = sizeof(uint8_t);  // load flag
    lh_bytes[10] = MEM_HEADER;       // alignment (custom use)
    int i = 0;
    for (i = 0; i < 10; i++)
      lh_bytes[10] -= lh_bytes[i]; // alignment (custom use)
    lh_bytes[11]  = sizeof(uint32_t); // data size
    lh_bytes[10] -= lh_bytes[11];    // alignment (custom use)

    lh_memoffset[0] = 0;             // writeversion
    for (i = 1; i < LH_PARAMS; i++)
      lh_memoffset[i] = lh_memoffset[i-1] + lh_bytes[i-1];

    lh_alignment = lh_bytes[10];
#if KALATTICEDEBUG
    fprintf(stderr, "lattice header memory consumption and offsets:\n");
    for (i = 0; i < LH_PARAMS; i++)
      fprintf(stderr, "lh_bytes[%2d] = %3d, lh_memoffset[%2d] = %3d\n", i, lh_bytes[i], i, lh_memoffset[i]);
#endif  // if KALATTICEDEBUG
  }  // Offsets

  void lhMatrixDimension(uint32_t v = 4) {*(uint32_t *)(lh_data + lh_memoffset[5]) = v;}

  void lhReadVersion(kaIOVersion v) {lh_readversion = v;}

  char *lh_data;                                     //!< Data of lattice header
  double lh_fallbackData[MEM_HEADER/sizeof(double)];  //! Tricky: alignment through double array !

  kaSharedMatrixN<double> lh_matrix;  //!< Matrix commonly linked in data of lattice header

  kaIOVersion lh_readversion;

  // uint32_t    lh_elemsize;
  uint32_t lh_alignment;
};  // class kaLatticeHeader


//! Class for basic handling of a kaLattice header in shared memory
/*!
   Internal use only!
   \n\n
   \author \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

class kaLatticeSMHeader : public kaLatticeHeader, public nskaIPC::IPCShm {
 public:
  //! Constructor
  kaLatticeSMHeader(kaLatticeCreate &cls, int32_t esize = 1, bool shm_check = false)  /////////////////
    : kaLatticeHeader(cls), nskaIPC::IPCShm(cls.Name, cls.xLattice != -1, shm_check) {  /////////////////
#if KALATTICEDEBUG
    cerr << "kaLatticeSMHeader(kaLatticeCreate &, int32_t)" << endl;
#endif  // if KALATTICEDEBUG
    using namespace nskaIPC;

    if (!FileCheck::Access())
      throw kaBaseException("Lattice file %s cannot be accessed", cls.Name);

    SharedMemorySize = 0;


    // protect the header initialization process
    try {
      int smsTrials = 0;
      SharedMemorySize = calcSizeOfShm(); smsTrials++;
#if KALATTICEDEBUG
      cerr << "kaLatticeSMHeader SharedMemorySize = " << SharedMemorySize << endl;
#endif  // if KALATTICEDEBUG

      // if size is 0, try from file
      if (SharedMemorySize == 0) {
        SharedMemorySize = calcSizeOfFile(); smsTrials++;
#if KALATTICEDEBUG
        cerr << "kaLatticeSMHeader SharedMemorySize = " << SharedMemorySize << endl;
#endif  // if KALATTICEDEBUG

        // if size is 0, calculate it from the data in the kaLatticeCreate object
        if (SharedMemorySize == 0) {
          if (isCreated()) {
            lhDim(cls, esize);
            SharedMemorySize = calcSMSize(); smsTrials++;
#if KALATTICEDEBUG
            cerr << "kaLatticeSMHeader: SharedMemorySize = " << SharedMemorySize << endl;
#endif  // if KALATTICEDEBUG
          } else {
            throw kaBaseException(
                    "Lattice file %s exists, but has no valid header and is not allocated in shared memory", cls.Name);
          }
        }
      }

      // create the shared memory
      if (!shmCreate(SharedMemorySize))
        throw InvalidIPC("Lattice file %s cannot be used to create shared memory", cls.Name);

      // initialize the header
      if (smsTrials == 2)
        memcpy(getAddress(), lhHeaderMem(), lhHeaderSize());
      kaLatticeHeader::Init(getAddress());
      if (numAttach() <= 1)
        // if the header is not already read from disk
        if (smsTrials != 2)
          // read header data from disk
          if (!kaLatticeHeader::headerRead())
            // if reading form disk failed, initialize header from the kaLatticeCreate object
            kaLatticeHeader::Define(cls, esize);

      // leave lattice file locked
    } catch (...) {
      IPCShm::cleanUp();
      IPCKey::cleanUp();
      kaLatticeHeader::Delete();
      throw;
    }
  }

 protected:
  uint32_t SharedMemorySize;
};  // class kaLatticeSMHeader

//! Class for handling of a kaLattice header in shared memory
/*!
   If the kaLattice header is already in shared memory, kaLatticeSharedHeader attaches directly.
   Otherwise, the kaLattice header is loaded in shared memory, which is set sufficiently large to include the whole
      lattice data.
   Therefore, a secondary process can successfully create an object of type kaLattice, while an object of type
      kaLatticeSharedHeader is already existing.
   \n\n
   \author \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */
class kaLatticeSharedHeader :  public kaLatticeSMHeader {
 public:
  //! Constructor. The argument cls delivers name and default sizes for creation of lattice.
  kaLatticeSharedHeader(kaLatticeCreate &cls, int32_t esize = 1) : nskaIPC::IPCKey(cls.Name, false), kaLatticeSMHeader(
      cls, esize) {
#if KALATTICEDEBUG
    cerr << "kaLatticeSharedHeader(kaLatticeCreate & cls, int32_t esize = 1)" << endl;
#endif  // if KALATTICEDEBUG

    // lattice file is still locked. Unlock and close it finally.
    unlockReadFD();
    closeReadFD();
  }
};

#endif  // ifndef KALATTICEHEADER_H
