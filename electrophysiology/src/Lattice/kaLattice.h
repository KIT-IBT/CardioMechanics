/*! \file kaLattice.h
   \brief Template class kaLattice

   Public interface for io and access of 3D array (lattice) in shared memory.
   A lattice includes a header data and a flattend 3D array, i.e. 1D array.

   \author cw,fs, IBT - Universität Karlsruhe (TH)
   \version 3.0
 */
#ifndef KALATTICE_H
#define KALATTICE_H

#include <kaMachineOS.h>
#include <kaPoint.h>
#include <kaSharedArray.h>
#include <kaLatticeHeader.h>

//! Class for handling of 3D array (lattice) in shared memory
/*!
   The array's size is given by three sizes in x-, y-, and z-direction.
   kaLattice allows to specify the type of array's elements via template X.
   The 3D array is flattened in a 1D array.
   The position, orientation and scaling of the array are defined by a homogenous matrix.

   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

template<class X>
class kaLattice : public kaLatticeSMHeader, public nskaGlobal::Raw<X>, public nskaGlobal::kaRLE<X>,
                  public nskaGlobal::GZIP<X> {
 public:
  kaSharedMatrixN<double> m;  //!< Homogeneous matrix (4x4) describing geometry of lattice.
  kaSharedArray<X> lat;  //!< 1D-array of lattice's elements, flattened 3D-array.

  IndexType xLattice;        //!< Size of lattice in x-direction.
  IndexType yLattice;        //!< Size of lattice in y-direction.
  IndexType zLattice;        //!< Size of lattice in z-direction.
  IndexType xyLattice;       //!< Size of a xy-slice. xyLattice=xLattice * yLattice.
  IndexType xyzLattice;  //!< Size of whole lattice. xyzLattice=xLattice * yLattice * zLattice.

  //! Constructor. The argument cls delivers name, default size and compression etc.
  kaLattice(kaLatticeCreate &cls, bool shm_check = false);  ////////////////////////////////////////////////////////////

  //! Destructor.
  virtual ~kaLattice(void);

  //! Get the value of the array's element at position (x,y,z).
  inline X & val(IndexType x, IndexType y, IndexType z)
  {return lat[x+y*xLattice+z*xyLattice];}
  //! restore lattice from disk to shared memory.
  bool Restore(void);
  //! save lattice from shared memory to disk.
  void Save(const char *NameArg = NULL);

  //! Set compression type for writing to disk.
  inline void SetCompression(CompressionType ct = ctNone) {
    lhCompression(ct);
  }

  //! Get compression type for writing to disk.
  inline CompressionType GetCompression(void)
  {return lhCompression();}

  //! Set version for writing of lattice to disk.
  inline void SetWriteVersion(kaIOVersion io = ioNow)
  {lhWriteVersion(io); *(kaIOVersion *)getAddress() = io;}

  //! Get version for writing of lattice to disk.
  inline kaIOVersion GetWriteVersion(void)
  {return lhWriteVersion();}

  //! Get version of lattice read from disk.
  inline kaIOVersion GetReadVersion(void)
  {return lhReadVersion();}

  //! Get data type of array elements.
  inline DataType GetDataType(void)
  {return lhDtype();}
  //! Copy lattice's elements from lattice.
  void Copy(kaLattice<X> *l);
  //! Copy region of lattice's elements from lattice.
  void Copy(kaLattice<X> *l, unsigned char, IndexType x1, IndexType y1, IndexType z1, IndexType x2, IndexType y2,
            IndexType z2, IndexType inc);
  //! Fill lattice's elements with value.
  void FillWith(X r);
  //! Fill region of lattice's elements with value.
  void FillWith(X r, IndexType x1, IndexType y1, IndexType z1, IndexType x2, IndexType y2, IndexType z2, IndexType inc);
  //! Get mean value of lattice elements. Only usefill with scalar types!
  double Mean();
  //! Generate meanfree lattice elements. Only usefill with signed data types!
  void MeanFree();
  //! Get gradient of lattice's elements at position x, y, z
  const kaPointDouble Sobel(IndexType x, IndexType y, IndexType z);

  //! Check of positions out of range and handling of errors via exception
  void CheckRange(IndexType x, IndexType y, IndexType z) {
    if ((x < 0) || (y < 0) || (z < 0) || (x >= xLattice) || (y >= yLattice) || (z >= zLattice) )
      throw kaBaseException("kaLattice access (%d,%d,%d) out of range 0-%d 0-%d 0-%d", x, y, z, xLattice-1, yLattice-1,
                            zLattice-1);
  }

 private:
  void Init(char *pCharArg);
  void Define(void);
  void Uninit(void);

  inline void PreInit(void) {
    xLattice  = yLattice = zLattice = 0;
    xyLattice = xyzLattice = 0;
  }
};  // class kaLattice


typedef kaLattice<uint8_t> Lattice;

/* ******************************************************************************************************* */
/*                                                                                                         */
/* description of the memory setup:                                                                        */
/*                                                                                                         */
/* the memory setup is close to the header structure.                                                      */
/* the following table shows the memory structure;                                                         */
/*                                                                                                         */
/* ------------------------------------------------------------------------------------------------------- */
/* data type       | variable name    | number of bytes needed        | description                        */
/* ------------------------------------------------------------------------------------------------------- */
/* (kaIOVersion)     | lhWriteVersion() | sizeof(kaIOVersion)       =   4 | file version                       */
/* uint32_t        | lhSMSize()       | sizeof(uint32_t)        =   4 | size of shared memory block        */
/* uint32_t        | xLattice         | sizeof(uint32_t)        =   4 | width of data block                */
/* uint32_t        | yLattice         | sizeof(uint32_t)        =   4 | height of data block               */
/* uint32_t        | zLattice         | sizeof(uint32_t)        =   4 | depth of data block                */
/* uint32_t        | m.Dimension()    | sizeof(uint32_t)        =   4 | dimension of transformation matrix */
/* double *        | m.a(row,column)  | 16 * sizeof(double)     = 128 | elements of matrix                 */
/* DataType        | lhDtype()        | sizeof(DataType)        =   4 | type of data                       */
/* CompressionType | lhCompression()  | sizeof(CompressionType) =   4 | compression method                 */
/* uint32_t        | lhElemSize()     | sizeof(uint32_t)        =   4 | size of one element                */
/* bool            | lhLoadFlag()     | sizeof(bool)            =   4 | flag for restore status of lattice */
/* --------------- | ---------------- | ----------------------------- | ---------------------------------- */
/*                 |                  |                     sum = 164 |                                    */
/* char *          | no name          |                            88 | alignment bytes (custom use)       */
/* uint32_t        | lhDsize()        | sizeof(uint32_t)        = 4   | number of voxels                   */
/* X *             | &lat[0]          | lh_dsize * sizeof(X)    = ?   | voxel data (aligned to 256 byte)   */
/* ------------------------------------------------------------------------------------------------------- */
/*                                                                                                         */
/* the memory header contains 256 byte: the first 164 bytes and the last 4 bytes are needed at the moment. */
/* the bytes from byte 164 up to byte 251 are free for user purposes.                                      */
/* the voxel data section starts at byte 256.                                                              */
/*                                                                                                         */
/* ******************************************************************************************************* */


template<class X>
kaLattice<X>::kaLattice(kaLatticeCreate &cls, bool shm_check) : nskaIPC::IPCKey(cls.Name, cls.xLattice != -1),
  kaLatticeSMHeader(cls, sizeof(X), shm_check) {  //////////////////////
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::kaLattice(kaLatticeCreate & cls)" << endl;
#endif  // if KALATTICEDEBUG

  // the file is still locked, leave it locked to protect the data initialization process

  // initialize the lattice data
  // restore the data if not already loaded
  try {
    PreInit();
    Init(getAddress());
    if (!lhLoadFlag())
      if (!Restore()) Define();
    lhLoadFlag(true);
    if (isCreated()) {
      lhDtype(GetTemplateDataType<X>());
      SetCompression(cls.compression);
    }
    if (GetTemplateDataType<X>() != lhDtype())
      throw TypeMismatch(TypeName(GetTemplateDataType<X>()), cls.Name);
  } catch (...) {
    nskaIPC::IPCKey::cleanUp();
#if KALATTICEDEBUG
    cerr << "kaLattice<X>::kaLattice(kaLatticeCreate & cls)::catch(...)" << endl;
#endif  // if KALATTICEDEBUG
    Uninit();
    /*mm,os otherwise if typemismatch==true ->malloc errors!; set all pointers to NULL before possible calls of
       destructors. Destructor is not called if constuctor failed!*/
    throw;
  }

  // the file is still locked. Unlock and close it finally
  nskaIPC::IPCKey::cleanUp();

#if KALATTICEDEBUG == 2
  kaLatticeHeader::FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 2
}

template<class X>
kaLattice<X>::~kaLattice(void) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::~kaLattice(void)" << endl;
#endif  // if KALATTICEDEBUG
  Uninit();
}

template<class X>
void kaLattice<X>::Init(char *psharedmem) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Init(char * psharedmem)" << endl;
#endif  // if KALATTICEDEBUG

  // initialize the private variables from the shared memory
  char *pmem    = psharedmem + lh_memoffset[0];
  char *pmatrix = psharedmem + lh_memoffset[5];
  char *pdata   = psharedmem + lh_memoffset[11];

  kaLatticeHeader::Init(pmem);
  xLattice   = lhWidth();
  yLattice   = lhHeight();
  zLattice   = lhDepth();
  xyLattice  = xLattice * yLattice;
  xyzLattice = xyLattice * zLattice;
  m.Attach(pmatrix, 4);

  lat.Attach(pdata, xyzLattice);

#if KALATTICEDEBUG == 2
  kaLatticeHeader::FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 2
}

template<class X>
void kaLattice<X>::Define(void) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Define(void)" << endl;
#endif  // if KALATTICEDEBUG

  // define the shared memory from the private variables
  memcpy(getAddress(), lhHeaderMem(), lhHeaderSize());

  xLattice   = lhWidth();
  yLattice   = lhHeight();
  zLattice   = lhDepth();
  xyLattice  = xLattice * yLattice;
  xyzLattice = xyLattice * zLattice;

  m.Identity();
}

template<class X>
void kaLattice<X>::Uninit(void) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Uninit(void)" << endl;
#endif  // if KALATTICEDEBUG

  // uninitialize matrix and array ...
  // kaLatticeHeader::Uninit();  --conform mit fs, mm
  m.Detach();
  lat.Detach();
}

template<class X>
bool kaLattice<X>::Restore(void) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Restore(void)" << endl;
#endif  // if KALATTICEDEBUG
  bool ret = false;

  // load data from disk into already allocated shared memory

  bool rmFD = false;
  if (!isReadFD()) {createReadFD(); rmFD = true;}
  if (!isReadFD())
    return false;

  int fd = getReadFD();

  char *pmem = getAddress();

  // restore the lattice header
  bool ret2 = kaLatticeHeader::Read();
  if (!ret2) {
    if (rmFD)
      closeReadFD(); return false;
  }

  // depending on the compression type restore the data respectively
  switch (lhReadVersion()) {
    case ioLast:  // Lattice Version 2
    case ioNow:  // Lattice Version 3
    {
      pmem += lhHeaderSize();

      switch (lhCompression()) {
        case ctNone: {
          Raw<X> *pIO = (Raw<X> *)this;
          pIO->Read((X *)pmem, calcNumVoxels(), fd);
        }
        break;
        case ctRLE: {
          kaRLE<X> *pIO = (kaRLE<X> *)this;
          pIO->Read((X *)pmem, calcNumVoxels(), fd);
        }
        break;
        case ctERLE:
          if (rmFD)
            closeReadFD();
          return false;

          break;
        case ctGZIP: {
          GZIP<X> *pIO = (GZIP<X> *)this;
          pIO->Read((X *)pmem, calcNumVoxels(), fd);
        }
        break;
        default:
          if (rmFD)
            closeReadFD();
          return false;
      }  // switch

      ret = true;
    }
    break;
    case ioInitial:  // Lattice Version 1
    case ioUnknown:
#if KALATTICEDEBUG
      cerr << "restore operation failed (lattice version " << lhReadVersion() << ") ! ]" << endl;
#endif  // if KALATTICEDEBUG
      if (rmFD)
        closeReadFD();
      throw kaBaseException("Lattice %s is of ancient file format and cannot be read", getName());
      break;
  }  // switch
  if (rmFD)
    closeReadFD();

#if KALATTICEDEBUG == 2
  kaLatticeHeader::FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 2
  return ret;
}  // >::Restore

template<class X>
void kaLattice<X>::Save(const char *nameArg) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Save(char *nameArg)" << endl;
#endif  // if KALATTICEDEBUG
  if (!nameArg)
    nameArg = getName();

  int fd = open(nameArg, O_WRONLY | O_CREAT | O_TRUNC, 0640);
  if (fd < 0)
    throw kaBaseException("Lattice file \"%s\" cannot be written", nameArg);

  kaLockFile(fd);

  if (!kaLatticeHeader::Write(fd)) {
    kaUnlockFile(fd);
    close(fd);
    throw kaBaseException("Cannot write header in lattice file \"%s\"!", nameArg);
  }

  switch (lhWriteVersion()) {
    case ioInitial:
      kaUnlockFile(fd);
      close(fd);
      throw kaBaseException("Cannot write in lattice file \"%s\" with ancient lattice file format!", nameArg);
      break;
    case ioLast:
      lhCompression(ctNone);
    case ioNow:
      switch (lhCompression()) {
        case ctNone: {
          Raw<X> *pIO = (Raw<X> *)this;
          pIO->Write(&lat[0], calcNumVoxels(), fd);
        }
        break;
        case ctRLE: {
          kaRLE<X> *pIO = (kaRLE<X> *)this;
          pIO->Write(&lat[0], calcNumVoxels(), fd);
        }
        break;
        case ctGZIP: {
          GZIP<X> *pIO = (GZIP<X> *)this;
          pIO->Write(&lat[0], calcNumVoxels(), fd);
        }
        break;
        case ctERLE:
        case ctUnknown:
          break;
      }
      break;
    case ioUnknown:
      kaUnlockFile(fd);
      close(fd);
      throw kaBaseException("Cannot write in file \"%s\" with unknown lattice file version!", nameArg);
      break;
  }  // switch
  kaUnlockFile(fd);
  close(fd);
#if KALATTICEDEBUG == 2
  kaLatticeHeader::FileHeaderOut(stderr);
#endif  // if KALATTICEDEBUG == 2
}  // >::Save

template<class X>
void kaLattice<X>::Copy(kaLattice<X> *l) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Copy(kaLattice<X> * l)" << endl;
#endif  // if KALATTICEDEBUG

  DimMatch(l);

  X *p = &lat[0];
  X *q = &lat[xyzLattice];
  X *r = &l->lat[0];

  while (p < q) *p++ = *r++;
}

template<class X>
void kaLattice<X>::Copy(kaLattice<X> *l, unsigned char, IndexType x1, IndexType y1, IndexType z1, IndexType x2,
                        IndexType y2, IndexType z2, IndexType inc) {
#if KALATTICEDEBUG
  cerr <<
    "kaLattice<X>::Copy(kaLattice<X> * l, unsigned char, IndexType x1, IndexType y1, IndexType z1, IndexType x2, IndexType y2, IndexType z2, IndexType inc)"
       << endl;
#endif  // if KALATTICEDEBUG

  DimMatch(l);

  CheckRange(x1, y1, z1);
  CheckRange(x2, y2, z2);

  for (IndexType z = z1; z < z2; z += inc)
    for (IndexType y = y1; y < y2; y += inc) {
      X *p = &val(x1, y, z);
      X *q = &val(x2, y, z);
      X *r = &l->val(x1, y, z);
      while (p < q) {
        *p = *r;
        p += inc;
        r += inc;
      }
    }
}

template<class X>
void kaLattice<X>::FillWith(X r) {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::FillWith(X r)" << endl;
#endif  // if KALATTICEDEBUG

  X *p = &lat[0];
  X *q = &lat[xyzLattice];

  while (p < q) *p++ = r;
}

template<class X>
void kaLattice<X>::FillWith(X r, IndexType x1, IndexType y1, IndexType z1, IndexType x2, IndexType y2, IndexType z2,
                            IndexType inc) {
#if KALATTICEDEBUG
  cerr <<
    "kaLattice<X>::FillWith(X r, IndexType x1, IndexType y1, IndexType z1, IndexType x2, IndexType y2, IndexType z2, IndexType inc)"
       << endl;
#endif  // if KALATTICEDEBUG

  CheckRange(x1, y1, z1);
  CheckRange(x2, y2, z2);

  for (IndexType z = z1; z < z2; z += inc)
    for (IndexType y = y1; y < y2; y += inc) {
      X *p = &val(x1, y, z);
      X *q = &val(x2, y, z);
      while (p < q) {
        *p = r;
        p += inc;
      }
    }
}

template<class X>
double kaLattice<X>::Mean() {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::Mean()" << endl;
#endif  // if KALATTICEDEBUG

  X *p = &lat[0];
  X *q = &lat[xyzLattice];

  double Sum = 0;
  while (p < q) Sum += *p++;
  Sum /= xyzLattice;

  return Sum;
}

template<class X>
void kaLattice<X>::MeanFree() {
#if KALATTICEDEBUG
  cerr << "kaLattice<X>::MeanFree()" << endl;
#endif  // if KALATTICEDEBUG

  X *p = &lat[0];
  X *q = &lat[xyzLattice];

  double _Mean = Mean();

  p = &lat[0];
  while (p < q) *p++ -= _Mean;
}

template<class X>
const kaPointDouble kaLattice<X>::Sobel(IndexType x, IndexType y, IndexType z) {
#if KALATTICEDEBUG == 2
  cerr << "kaLattice<X>::Sobel(IndexType x, IndexType y, IndexType z)" << endl;
#endif  // if KALATTICEDEBUG == 2

  kaPointDouble g(0, 0, 0);
  static const float maskz[27] =
  {-0.5, -1, -0.5, -1, -2, -1, -0.5, -1, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 0.5, 1, 2, 1, 0.5, 1, 0.5};
  static const float masky[27] =
  {-0.5, -1, -0.5, 0, 0, 0, 0.5, 1, 0.5, -1, -2, -1, 0, 0, 0, 1, 2, 1, -0.5, -1, -0.5, 0, 0, 0, 0.5, 1, 0.5};
  static const float maskx[27] =
  {-0.5, 0, 0.5, -1, 0, 1, -0.5, 0, 0.5, -1, 0, 1, -2, 0, 2, -1, 0, 1, -0.5, 0, 0.5, -1, 0, 1, -0.5, 0, 0.5};

  int i = 0;
  for (IndexType gz = z-1; gz <= z+1; gz++)
    for (IndexType gy = y-1; gy <= y+1; gy++)
      for (IndexType gx = x-1; gx <= x+1; gx++, i++)
        if ((gz >= 0) && (gz < zLattice))
          if ((gy >= 0) && (gy < yLattice))
            if ((gx >= 0) && (gx < xLattice)) {
              X val = lat[gx+gy*xLattice+gz*xyLattice];
              g.x += maskx[i] * val;
              g.y += masky[i] * val;
              g.z += maskz[i] * val;
            }

  return g;
} // >::Sobel

#endif  // ifndef KALATTICE_H
