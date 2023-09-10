/**@file SharedMemory.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef SHAREDMEMORY_H
#define SHAREDMEMORY_H

#include <kaMachineOS.h>
#include <kaExceptions.h>

#ifdef osWin32
# include <afxwin.h> // MFC-Kern- und -Standardkomponenten
#endif  // ifdef osWin32

class SharedMemory {
 public:
  char Name[512];

#ifndef osWin32
  int ID;
#else  // ifndef osWin32
  HANDLE ID;
  HANDLE attID;
  char *attAddress;
  unsigned long size;
  char mapname[512];
  char attmapname[512];
#endif  // ifndef osWin32
  char *Address;
  struct shmid_ds Info;

  SharedMemory(const char *n) {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::SharedMemory %s", n);
#endif  // if KALATTICEDEBUG

    Address = NULL;
    strcpy(Name, n);
  }

  ~SharedMemory() {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::~SharedMemory");
#endif  // if KALATTICEDEBUG
  }

  int Create(unsigned int Size, int i = 0) {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Create");
#endif  // if KALATTICEDEBUG

#ifdef osWin32

    size = Size;

    // create the mapping names
    CreateMappingID(Name, mapname);
    strcpy(attmapname, mapname);
    strcat(attmapname, "attached");

    // create file mapping for the file
    ID = CreateFileMapping((HANDLE)0xFFFFFFFF, NULL, PAGE_READWRITE, 0, size, mapname);

    if (GetLastError() == ERROR_ALREADY_EXISTS) {
      ID = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, mapname);

      // map file exists -> read number of clients and add one
      attID        = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, attmapname);
      attAddress   = (char *)MapViewOfFile(attID, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, sizeof(char));
      *attAddress += 1;
    } else {
      // create file mapping for number of clients
      attID       = CreateFileMapping((HANDLE)0xFFFFFFFF, NULL, PAGE_READWRITE, 0, sizeof(char), attmapname);
      attAddress  = (char *)MapViewOfFile(attID, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, sizeof(char));
      *attAddress = 1;
    }

    if (ID == NULL)
      throw kaBaseException("SharedMemory::Create CreateFileMapping Error");

    Info.shm_segsz  = size;
    Info.shm_nattch = *attAddress;

# if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Create ID=%i Size=%d\n", ID, (int)Size);
# endif  // if KALATTICEDEBUG

    return (int)ID;

#else  // ifdef osWin32
    key_t key = kaFtok(Name, i);
    if (key == -1) {
# if KALATTICEDEBUG
      fprintf(stderr, "SharedMemory::Create key\n");
# endif  // if KALATTICEDEBUG
      return -1;
    }
# if !defined(NO_SHM)
    ID = shmget(key, (size_t)Size, 0666 | IPC_CREAT);
# else  // if !defined(NO_SHM)
    Address = (char *)malloc(Size);
# endif  // if !defined(NO_SHM)
    if (ID == -1) if (errno == EINVAL) {
# if KALATTICEDEBUG
        fprintf(stderr, "SharedMemory::Create error EINVAL\n");
# endif  // if KALATTICEDEBUG
      }
# if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Create ID=%i key=%i Size=%d\n", ID, key, Size);
# endif  // if KALATTICEDEBUG

    return ID;

#endif  // ifdef osWin32
  }  // Create

  int GetID(int i = 0) {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::GetID %d", i);
#endif  // if KALATTICEDEBUG

#ifdef osWin32
    return (int)ID;

#else  // ifdef osWin32
    key_t key = kaFtok(Name, i);
# if !defined(NO_SHM)
    ID = shmget(key, 0, 0);
# endif  // if !defined(NO_SHM)

# if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::GetID  ID=%i key=%i\n", ID, key);
# endif  // if KALATTICEDEBUG

    return ID;

#endif  // ifdef osWin32
  }

  void Attach() {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Attach");
#endif  // if KALATTICEDEBUG

#ifdef osWin32
    Address = (char *)MapViewOfFile(ID, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, size);
#else  // ifdef osWin32
# if !defined(NO_SHM)
    Address = (char *)shmat(ID, NULL, 0);
# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32
    if (!Address)
      throw kaBaseException("SharedMemory::Attach: Address==NULL");

#if KALATTICEDEBUG
    fprintf(stderr, "SM attach ID=%i  Address %ld\n", ID, Address);
#endif  // if KALATTICEDEBUG
  }

  void Control() {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Control");
#endif  // if KALATTICEDEBUG

#ifdef osWin32

    // check for the numbers of attached clients
    Info.shm_nattch = *attAddress;
# if KALATTICEDEBUG
    fprintf(stderr, "SM control ID=%i\n", ID);
# endif  // if KALATTICEDEBUG

#else  // ifdef osWin32
# if !defined(NO_SHM)
    shmctl(ID, IPC_STAT, &Info);
#  if KALATTICEDEBUG
    fprintf(stderr, "SM control ID=%i Attache Level=%i Size in bytes=%i\n", ID, Info.shm_nattch, Info.shm_segsz);
#  endif  // if KALATTICEDEBUG
# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32
  }

  void Detach() {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Detach");
#endif  // if KALATTICEDEBUG

#ifdef osWin32

    // remove this client by decreasing the number of clients
    // make Address and ID to NULL;

    *attAddress    -= 1;
    Info.shm_nattch = *attAddress;
    if (Info.shm_nattch > 0) {
      Address    = NULL;
      ID         = NULL;
      attID      = NULL;
      attAddress = NULL;
    }
#else  // ifdef osWin32
# if !defined(NO_SHM)
    shmdt(Address);
# else  // if !defined(NO_SHM)
    if (Address) {
      free(Address);
      Address = NULL;
    }

# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32

#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Detach ID=%i\n", ID);
#endif  // if KALATTICEDEBUG
  }  // Detach

  void Destroy() {
#if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Destroy");
#endif  // if KALATTICEDEBUG

#ifdef osWin32

    // check if number of clients are zero, if yes -> destroy
    if (Info.shm_nattch == 0) {
      // destroy
      if (Address) {
        UnmapViewOfFile(Address);
        Address = NULL;
      }

      if (ID) {
        CloseHandle(ID);
        ID = NULL;
      }

      if (attAddress) {
        UnmapViewOfFile(attAddress);
        attAddress = NULL;
      }

      if (attID) {
        CloseHandle(attID);
        attID = NULL;
      }
    }
# if KALATTICEDEBUG
    fprintf(stderr, "SharedMemory::Destroy destroyed ID=%i\n", ID);
# endif  // if KALATTICEDEBUG

#else  // ifdef osWin32
# if !defined(NO_SHM)
    shmctl(ID, IPC_STAT, &Info);
    if (Info.shm_nattch == 0) {
      shmctl(ID, IPC_RMID, 0);
#  if KALATTICEDEBUG
      fprintf(stderr, " SharedMemory::Destroy destroyed ID=%i\n", ID);
#  endif  // if KALATTICEDEBUG
    } else {
#  if KALATTICEDEBUG
      fprintf(stderr, "SharedMemory::Destroy not destroyed ID=%i, Attach Level=%i\n", ID, Info.shm_nattch);
#  endif  // if KALATTICEDEBUG
    }
# endif  // if !defined(NO_SHM)
#endif  // ifdef osWin32
  }  // Destroy
};  // class SharedMemory

#endif  // ifndef SHAREDMEMORY_H
