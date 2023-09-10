/*! \file kaBasicIO.h
   \brief Basic functions for lattice-file input-output

   Internal use only!
   Known bugs: osWin32 and osLinux versions of Raw, RLE, and GZIP read
   and write methods as well as the functions kaWrite and
   kaRead modify data in memory temporarily. This may interfere with
   further (shared) memory operations.

   \author cw,idb IBT - Universität Karlsruhe (TH), fs CVRTI - University of Utah
 */

#ifndef KABASICIO_H
#define KABASICIO_H

#include <kaMachineOS.h>
#include <kaExceptions.h>
#include <zlib.h>

//! Basic functions for kaLattice
namespace nskaGlobal {
#ifdef KA_LITTLE_ENDIAN

//! Swap function allows to adapt to platform dependent byteorder of datatypes
extern void Swap(unsigned char *data, unsigned int size);
extern void Swap(signed char *data, unsigned int size);
extern void Swap(char *data, unsigned int size);

template<class X> void Swap(X *data, unsigned int size) {
  if (sizeof(X) == 1)
    return;

  X *pdata = data;
  X *edata = pdata + size;

  while (pdata < edata) {
    unsigned char *adata = (unsigned char *)pdata;
    unsigned char *bdata = adata + sizeof(X) - 1;
    for (unsigned int i = 0; i < sizeof(X)/2; i++) {
      char help = *adata;
      *adata++ = *bdata;
      *bdata-- = help;
    }
    pdata++;
  }
}

#endif  // ifdef KA_LITTLE_ENDIAN

//! Platform independent basic write function
template<class X> inline void kaWrite(X *data, int size, size_t nitems, FILE *stream) {
#ifdef KA_LITTLE_ENDIAN
  Swap(data, nitems);
  size_t rc = fwrite((void *)data, size, nitems, stream);
  Swap(data, nitems);
#else  // ifdef KA_LITTLE_ENDIAN
  size_t rc = fwrite((void *)data, size, nitems, stream);
#endif  // ifdef KA_LITTLE_ENDIAN
  if (rc != nitems)
    throw kaBaseException("kaWrite. Cannot write");
}

//! Platform independent basic write function
template<class X> inline void kaWrite(X *data, size_t nitems, FILE *stream) {
  kaWrite(data, sizeof(X), nitems, stream);
}

//! Platform independent basic write function
template<class X> inline void kaWrite(X data, FILE *stream) {
  kaWrite(&data, 1, stream);
}

//! Platform independent basic write function
template<class X> inline void kaWrite(X *data, size_t nitems, int fd) {
#ifdef KA_LITTLE_ENDIAN
  Swap(data, nitems);
  size_t rc = write(fd, (void *)data, nitems*sizeof(X));
  Swap(data, nitems);
#else  // ifdef KA_LITTLE_ENDIAN
  size_t rc = write(fd, (void *)data, nitems*sizeof(X));
#endif  // ifdef KA_LITTLE_ENDIAN
  if (rc != nitems*sizeof(X))
    throw kaBaseException("kaWrite. Cannot write");
}

//! Platform independent basic write function
template<class X> inline void kaWrite(X data, int fd) {
  kaWrite(&data, 1, fd);
}

//! Platform independent basic read function
template<class X> inline void kaRead(X *data, int size, size_t nitems, FILE *stream) {
  size_t rc = fread((void *)data, size, nitems, stream);

  if (rc != nitems)
    throw kaBaseException("kaRead. Cannot read");

#ifdef KA_LITTLE_ENDIAN
  Swap(data, nitems);
#endif  // ifdef KA_LITTLE_ENDIAN
}

template<class X> inline void kaRead(X *data, size_t nitems, FILE *stream) {
  kaRead(data, sizeof(X), nitems, stream);
}

//! Platform independent basic read function
template<class X> inline void kaRead(X *data, FILE *stream) {
  kaRead(data, 1, stream);
}

//! Platform independent basic read function
template<class X> inline void kaRead(X *data, size_t nitems, int fd) {
  size_t rc = read(fd, (void *)data, sizeof(X)*nitems);

#ifdef KA_LITTLE_ENDIAN
  Swap(data, nitems);
#endif  // ifdef KA_LITTLE_ENDIAN
  if (rc != sizeof(X)*nitems)
    throw kaBaseException("kaRead. Cannot read");
}

//! Platform independent basic read function
template<class X> inline void kaRead(X *data, int fd) {
  kaRead(data, 1, fd);
}

//! Class for writing and reading of raw data to and from disk
/*!
   \author cw,fs, IBT - Universität Karlsruhe (TH)
 */

template<class X>
class Raw {
 public:
  inline void Write(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "w");

    fseek(stream, pos, SEEK_SET);

    Write(data, nitems, stream);

    fclose(stream);
  }

  inline void Write(X *data, unsigned int nitems, FILE *stream) {
    kaWrite(data, nitems, stream);
  }

  inline void Read(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "r");

    fseek(stream, pos, SEEK_SET);

    Read(data, nitems, stream);

    fclose(stream);
  }

  inline void Read(X *data, unsigned int nitems, FILE *stream) {
    kaRead(data, nitems, stream);
  }
};  // class Raw

//! Read buffer size for class GZIP
const int ibufsize = 32768;  // 32*1024;

//! Write buffer size for class GZIP
const int obufsize = 39337;  // (int)(32*1024*1.2+16);

//! Class for writing and reading of gzip-compressed data to and from disk.
/*!
   The class uses functions of zlib (described at http://www.gzip.org/zlib/)

   \author fs, IBT - Universität Karlsruhe (TH)
 */

template<class X>
class GZIP {
 public:
  inline void Write(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "w");

    fseek(stream, pos, SEEK_SET);

    Write(data, nitems, stream);

    fclose(stream);
  }

  inline void Write(X *data, unsigned int nitems, FILE *stream) {
    int len = nitems*sizeof(X), pos;
    Bytef buf[obufsize], *src = (Byte *)data;

#ifdef KA_LITTLE_ENDIAN
    Swap(data, nitems);
#endif  // ifdef KA_LITTLE_ENDIAN

    for (pos = 0; pos < len;) {
      int sizew = ibufsize;
      if (pos+sizew > len)
        sizew = len-pos;

      uLongf olen = obufsize;
      int rc      = compress(buf, &olen, (Bytef *)&src[pos], sizew);
      if (rc != Z_OK)
        throw kaBaseException("GZIP::Write");
      kaWrite<int>(olen, stream);
      rc = fwrite(buf, 1, olen, stream);
      if (rc != olen)
        throw kaBaseException("GZIP::Write");
      pos += sizew;
    }

#ifdef KA_LITTLE_ENDIAN
    Swap(data, nitems);
#endif  // ifdef KA_LITTLE_ENDIAN
  }  // Write

  inline void Read(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "r");

    fseek(stream, pos, SEEK_SET);

    Read(data, nitems, stream);

    fclose(stream);
  }

  inline void Read(X *data, unsigned int nitems, FILE *stream) {
    int len = nitems*sizeof(X), pos;
    Bytef buf[obufsize], *dest = (Byte *)data;

    for (pos = 0; pos < len;) {
      int rlen;
      kaRead<int>(&rlen, stream);
      if (rlen >= obufsize)
        throw kaBaseException("GZIP::Read");
      int rc = fread(buf, 1, rlen, stream);
      if (rc <= 0)
        throw kaBaseException("GZIP::Read");
      uLongf sizew = len-pos;

      rc = uncompress((Bytef *)&dest[pos], &sizew, buf, rc);
      if (rc != Z_OK)
        throw kaBaseException("GZIP::Read");
      pos += sizew;
    }

#ifdef KA_LITTLE_ENDIAN
    Swap(data, nitems);
#endif  // ifdef KA_LITTLE_ENDIAN
  }
};  // class GZIP


//! Class for writing and reading of rle-compressed data to and from disk
/*!
   rle (run length encoding) compression takes usage of repetitions.
   N-times repeating data elements (a,a,a, ..., a) are stored as (a,N).

   \author cw, IBT - Universität Karlsruhe (TH), fs - CVRTI - University of Utah
 */

template<class X> class kaRLE {
 public:
  inline void Write(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "w");

    fseek(stream, pos, SEEK_SET);

    Write(data, nitems, stream);

    fclose(stream);
  }

  inline void Write(X *data, unsigned int nitems, FILE *stream) {
    dataptr    = data;
    maxdataptr = data + nitems;
    Write(stream);
  }

  inline void Read(X *data, unsigned int nitems, int fd) {
    int pos      = lseek(fd, 0, SEEK_CUR);
    FILE *stream = fdopen(dup(fd), "r");

    fseek(stream, pos, SEEK_SET);

    Read(data, nitems, stream);

    fclose(stream);
  }

  inline void Read(X *data, unsigned int nitems, FILE *stream) {
    dataptr    = data;
    maxdataptr = data + nitems;
    Read(stream);
  }

 private:
  X *dataptr, *maxdataptr;

  int Search(int pos, int MAX_RETURN = 127) {
    X *pdataptr  = &dataptr[pos];
    X *ppdataptr = pdataptr;

    pdataptr++;

    while ((pdataptr < maxdataptr) && (pdataptr-ppdataptr <= MAX_RETURN)) {
      if (*pdataptr != *ppdataptr)
        break;
      pdataptr++;
    }

    return pdataptr-ppdataptr-1;
  }

  void Write(FILE *stream) {
    const int RLE_MIN_CRUNCH = 3;

    int nbytes = maxdataptr - dataptr;

    for (int i = 0; i < nbytes; i++) {
      int n = Search(i);

      if (n >= RLE_MIN_CRUNCH) {
        putc(128+n, stream);
        kaWrite(&dataptr[i], 1, stream);
        i += n;
      } else {
        int k    = i;
        int maxk = min(i+128, nbytes);

        while (k < maxk) {
          int l = Search(k, RLE_MIN_CRUNCH+1);
          k++;
          if (l >= RLE_MIN_CRUNCH)
            break;
        }

        int m = min(k-i-1, 127);

        putc(m, stream);
        kaWrite(&dataptr[i], m+1, stream);
        i += m;
      }
    }
  }  // Write

  void Read(FILE *stream) {
    int i = 0;

    while (!feof(stream)) {
      int ch = getc(stream);

      if (ch > 127) {
        int n = ch-127;
        if (dataptr+n > maxdataptr)
          n = maxdataptr-dataptr;

        X tmp;
        kaRead(&tmp, 1, stream);
        for (; n--; i++)
          dataptr[i] = tmp;
      } else {
        int n = ch+1;
        if (dataptr+n > maxdataptr)
          n = maxdataptr-dataptr;

        kaRead(&dataptr[i], n, stream);
        i += n;
      }
    }
  }
};  // class kaRLE
}  // namespace nskaGlobal
#endif  // ifndef KABASICIO_H
