/*
 * File: ShMemory.h
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


#ifndef SHMEMORY_H
#define SHMEMORY_H

#include <kaMachineOS.h>
#include <kaFileCheck.h>
#include <kaExceptions.h>

using namespace nskaGlobal;

#ifdef osWin32
# include <afxwin.h> // MFC-Kern- und -Standardkomponenten
typedef unsigned long shmatt_t;
#endif  // ifdef osWin32

typedef enum { shmSnooze = 0, shmDefine = 1, shmInit = 2 } shmNewAction;


class ShMIDErr : public kaBaseException {
 public:
  ShMIDErr() : kaBaseException("ShMemory: Could not get shared memory id!") {}
};


class ShMFileIOErr : public kaBaseException {
 public:
  ShMFileIOErr(const char *file) : kaBaseException("ShMemory: Could not access file \"%s\"!", file) {}
};


class ShMInvalidSize : public kaBaseException {
 public:
  ShMInvalidSize(int size) : kaBaseException("ShMemory: Shared memory size is invalid (%d)!", size) {}
};


class ShMInvalidFile : public kaBaseException {
 public:
  ShMInvalidFile(const char *fname) : kaBaseException("ShMemory: File name %s is invalid!", fname) {}
};


class ShMAlreadyAllocated : public kaBaseException {
 public:
  ShMAlreadyAllocated(char *fname, unsigned int size, char *mem) : kaBaseException(
      "ShMemory: Memory already allocated (\"%s\", %ld (%d bytes))!", fname, mem, size) {}
};


class ShMemory : virtual public nskaGlobal::FileCheck {
 public:
  ShMemory(const char *file) : FileCheck() {PreInit(); New(file);}

  ShMemory(void) : FileCheck() {PreInit();}

  ~ShMemory(void) {Delete();}

  inline shmNewAction New(const char *file) {
    try {
      setupAction = shmSnooze;
      if (shmSize == 0)
        throw ShMInvalidSize(0);

      if (memory)
        throw ShMAlreadyAllocated(filename, shmSize, memory);
      strcpy(filename, file);
      if (strlen(filename) == 0)
        throw ShMInvalidFile(filename);

      Check(filename);

      if (!Access())
        throw ShMFileIOErr(filename);
      if (!Exist()) {
        FILE *out = fopen(filename, "wb"); if (!out)
          throw ShMFileIOErr(filename); fclose(out);
      }

#ifdef osWin32
      shmID = NULL;
      shmID = CreateFileMapping((HANDLE)0xFFFFFFFF,
                                NULL,
                                PAGE_READWRITE,
                                0,
                                shmSize,
                                filename);

      if (shmID == NULL)
        throw kaBaseException("CreateFileMapping Error!");

      shmInfo.shm_nattch = 1;

      memory = (char *)MapViewOfFile(shmID,
                                     FILE_MAP_WRITE | FILE_MAP_READ,
                                     0, 0, shmSize);
      if (memory == NULL)
        throw kaBaseException("MapViewOfFile Error!");

      success     = true;
      setupAction = shmInit;
      return setupAction;

#else  // ifdef osWin32
      key_t key = kaFtok(filename, 0);
# if !defined(NO_SHM)
      shmID = shmget(key, shmSize, 0666 | IPC_CREAT);
      if (shmID != -1) {
        memory = (char *)shmat(shmID, NULL, 0);
        shmctl(shmID, IPC_STAT, (struct shmid_ds *)&shmInfo);
        success = true;
        if (shmInfo.shm_nattch == 1)
          setupAction = shmDefine; else
          setupAction = shmInit;
      } else {
        shmID = shmget(key, 0, 0);
        if (shmID) {
          memory = (char *)shmat(shmID, NULL, 0);
          shmctl(shmID, IPC_STAT, (struct shmid_ds *)&shmInfo);
          success     = true;
          setupAction = shmInit;
        } else {throw ShMIDErr();}
      }
      return setupAction;

# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32
    } catch (ShMIDErr &ie) {
      Destroy();
      throw;
    } catch (...) {
      throw;
    }
    return shmSnooze;
  }  // New

  inline void Delete(void) {
#ifndef osWin32
# if !defined(NO_SHM)
    shmdt(memory);
# endif  // if !defined(NO_SHM)
#endif  // ifndef osWin32
    Destroy();
    PreInit();
  }

  inline char *GetMem(void) {return memory;}

  inline struct kaShmDS & GetShMInfo(void) {return shmInfo;}

  inline unsigned int GetSize(void) {return shmSize;}

  inline char *GetName(void) {return filename;}

  inline long Align(long v, long bytes) {
    if (bytes == 0)
      return 0; unsigned int i = 0; while (v%bytes) {v++; i++;}
    return i;
  }

  inline long Align(char *pmem, long bytes) {long v = (long)pmem; return Align(v, bytes);}

  inline void ResetSize(void) {
    if (!success)
      shmSize = 0;
  }

  inline unsigned int ExpandSize(unsigned int items, unsigned int size = 1) {
    if (!success) {
      if ((size > 0) && (items > 0)) {
        shmSize += Align(shmSize, size);
        shmSize += items * size;
      }
    }
    return shmSize;
  }

  inline shmNewAction HowToSetup(void) {return setupAction;}

 private:
  inline void Destroy(void) {
#ifdef osWin32
    if (memory) {
      UnmapViewOfFile(memory);
      memory = NULL;
    }

    if (shmID) {
      CloseHandle(shmID);
      shmID = NULL;
    }
#else  // ifdef osWin32
# if !defined(NO_SHM)
    shmctl(shmID, IPC_STAT, (struct shmid_ds *)&shmInfo);
    if (shmInfo.shm_nattch == 0)
      shmctl(shmID, IPC_RMID, 0);
# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32
  }

  inline void PreInit(void) {
    success     = false;
    shmID       = 0;
    memory      = NULL;
    filename[0] = 0;
    ResetSize();
    setupAction = shmSnooze;
  }

  char *memory;
  char filename[256];
  struct kaShmDS shmInfo;
  unsigned int shmSize;
#ifdef osWin32
  HANDLE shmID;
#else  // ifdef osWin32
  int shmID;
#endif  // ifdef osWin32
  shmNewAction setupAction;

  bool success;
};  // class ShMemory


#endif  // ifndef SHMEMORY_H
