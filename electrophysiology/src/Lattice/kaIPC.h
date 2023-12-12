/*
 * File: kaIPC.h
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


#ifndef KAIPC_H
#define KAIPC_H

#include <kaMachineOS.h>
#include <kaExceptions.h>
#include <kaFileCheck.h>

//! Basic handling of IPC mechanisms (key, shared memory, and semphore)
namespace nskaIPC {

  //! Class for exception handling in IPC classes
  /*!
    An object of class InvalidIPC is typically thrown in case of errors with IPC functions.
    The object includes a string describing details of the error's cause.
    \n\n
    \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
  */
  //Class definition for InvalidIPC was removed due to a problem with inheritance from kaBaseException and the ... macro that was used to pass arguments \author dk091, gs, tf102
  #define InvalidIPC kaBaseException

  //! Class for exception handling in class IPCkey
  /*!
    An object of class InvalidIPC is thrown in case of non existance of file and creation flag equal false.
    The object includes a string describing details of the error's cause.
    \n\n
    \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
  */

	//Class definition for FileNotExisting was removed due to a problem with inheritance from kaBaseException and the ... macro that was used to pass arguments \author dk091, gs, tf102
	#define FileNotExisting kaBaseException

  //! Class for creating of an IPC key.
  /*!
    The creation of an IPC key uses a file given by name.
    The file is created if not already existing.
    The file is exclusively opened by locking.
    \n\n
    \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
  */

  class IPCKey : virtual public nskaGlobal::FileCheck
    {
    public:
      //! Default constructor.
      IPCKey(void) {
        init();
      }

      //! Constructor with file name and file creation flag.
      IPCKey(const char * fname, bool create) {
#if KALATTICEDEBUG
	cerr << "IPCKey::IPCKey(" << fname << ", " << create << ")" << endl;
#endif
	init();
	keyCreate(fname, create);
      }

      //! Destructor
      virtual ~IPCKey(void) {
#if KALATTICEDEBUG
	cerr << "IPCKey::~IPCKey()" << endl;
#endif
	cleanUp();
	if (name) free(name);
      }

      //! Create ipc key with file name and file creation flag.
      inline void keyCreate(const char *fname, bool create) {
#if KALATTICEDEBUG
	cerr << "IPCKey::keyCreate(" << fname << ", " << create << ")" << endl;
#endif
	if (name) free(name);
	//strdup allocate memory with malloc() and copy the char* to that memory location and return the address to the memory, memory can be deallocated with free()
	name=strdup(fname);
	Check(name);

	if(!Exist()) {
	  if (!create)
	    throw FileNotExisting("File %s does not exist!", fname);

	  //O_EXCL
	  //When used with O_CREAT, if the file already exists it is an error and the open() will fail.
	  //In this context, a symbolic link exists, regardless of where it points to. O_EXCL is broken on NFS file systems;
	  //programs which rely on it for performing locking tasks will contain a race condition.
	  //The solution for performing atomic file locking using a lockfile is to create a unique file on the same file system (e.g., incorporating hostname and pid), use link(2) to make a link to the lockfile.
	  //If link() returns 0, the lock is successful. Otherwise, use stat(2) on the unique file to check if its link count has increased to 2, in which case the lock is also successful.
	  readFD = open(name, O_CREAT | O_EXCL, 0640);
	  if(readFD >= 0)
	    created = true;
	}
	//if already exists try to open it
	createReadFD();
	if(readFD < 0) throw InvalidIPC("File %s is not readable!", fname);
	//exclusive locking of readFD
	kaLockFile(readFD);

	key = kaFtok(name, 0);
	if(key == -1) throw InvalidIPC("Bad key for file %s!", fname);
      }

      //! Get ipc key created by function keyCreate.
      inline key_t  getKey(void) const { return key; }

      //! Get file name associated to ipc key.
      inline const char * getName(void) const { return name; }

      //! Get creation flag concerning file of ipc key.
      inline bool isCreated(void) const { return created; }

    protected:
      inline bool isReadFD(void) const  { return ((readFD < 0) ? false : true); }
      inline int createReadFD(void) {
	if(readFD < 0)
	  readFD = open(name, O_RDONLY);
	return readFD;
      }
      inline void unlockReadFD(void) {
	if(readFD >= 0)
	  kaUnlockFile(readFD);
      }
      inline void closeReadFD(void) {
	if(readFD >= 0) {
	  close(readFD);
	  readFD=-1;
	}
      }
      inline int getReadFD(void) const { return readFD; }

      inline void cleanUp(void) {
	unlockReadFD();
	closeReadFD();
      }

    private:
      inline void init(void) {
	name=NULL;
	key = -1;
	readFD = -1;
	created = false;
      }

      key_t  key;
      char   *name;
      int    readFD;
      bool   created;
    };

#ifndef osWin32
  const int maxSemsPerID = 25;

  //! Class for creating of an IPC semaphore.
  /*!
    The creation of an IPC semaphore uses an IPC key.
    \n\n
    \author cw,fs, IBT - Universität Karlsruhe (TH)
  */
  class IPCSem : virtual public IPCKey
  {
  public:
    //! Default constructor.
    IPCSem(void) : IPCKey() { init(); }

    //! Constructor with file name and file creation flag.
    IPCSem(const char * fname, bool create) : IPCKey(fname, create) { init(); }

    //! Destructor.
    virtual ~IPCSem(void) { cleanUp(); }

    //! Create n semaphores.
    inline void semCreate(int n = 1)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::semCreate()" << endl;
#endif
      if(n < 1)
	throw InvalidIPC("Invalid numbers of semaphores (too low, < 1)");
      if(n >= maxSemsPerID-1)
	throw InvalidIPC("Invalid numbers of semaphores (too high, >= %d)", maxSemsPerID-1);

      numSem = n+1;
      if(getKey() >= 0)
        {
          semid = semget(getKey(), numSem, IPC_CREAT | 0666 | IPC_EXCL);
          if(semid < 0)
            {
              semid = semget(getKey(), 0, 0666);
              if(semid < 0)
		throw InvalidIPC("Semaphore not created. Bad ID.");
            }
          else
            {
              for(int i = 1; i < numSem; i++) vOperator(i);
#if KALATTICEDEBUG
              for(int ii = 0; ii < numSem; ii++)
                {
                  forthArgVal = semctl(semid, ii, GETVAL);
                  cerr << "nskaIPC::IPCSem::semCreate(int n = 1): semaphore id: " << semid;
                  cerr << ", #: " << ii << ", value: " << forthArgVal << endl;
                }
#endif
            }
        }
      else
	throw InvalidIPC("Semaphore not created. Bad key.");

      vOperator(0);
#if KALATTICEDEBUG
      forthArgVal = semctl(semid, 0, GETVAL);
      cerr << "nskaIPC::IPCSem::semCreate(int n = 1): semaphore id: " << semid << ", #: " << 0 << ", value: " << forthArgVal << endl;
      forthArgVal = semctl(semid, 0, GETNCNT);
      cerr << "nskaIPC::IPCSem::semCreate(int n = 1): semaphore id: " << semid << ", #: " << 0 << ", semncnt: " << forthArgVal << endl;
#endif
    }

    //! P()
    inline bool semDeplete(int i, int flag = 0)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::semDeplete(" << i << ", " << flag << ")" << endl;
#endif
      i++;
      if(i >= 1 && i < numSem)
	return pOperator(i, flag);
      else
	throw InvalidIPC("Invalid semaphore number.");
    }

    //! V()
    inline bool semRevive(int i, int flag = 0)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::semRevive(" << i <<")" << endl;
#endif
      i++;
      if(i >= 1 && i < numSem)
	return vOperator(i, flag);
      else
	throw InvalidIPC("Invalid semaphore number.");
    }

  protected:
    inline void       cleanUp(void)
    {
      forthArgVal = semctl(semid, 0, GETVAL);
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::cleanUp(void): semaphore id: " << semid << ", #: " << 0 << ", value: " << forthArgVal << endl;
#endif
      if(forthArgVal <= 1)
        {
          if(semid >= 0)
            {
#if KALATTICEDEBUG
              cerr << "nskaIPC::IPCSem::cleanUp(void): removing semaphore identifier" << endl;
#endif
              semctl(semid, 0, IPC_RMID);
            }
        }
      else
        {
          pOperator(0);
#if KALATTICEDEBUG
          forthArgVal = semctl(semid, 0, GETVAL);
          cerr << "nskaIPC::IPCSem::cleanUp(void): semaphore id: " << semid << ", #: " << 0 << ", value: " << forthArgVal << endl;
#endif
        }
    }

  private:
    inline void       init(void)
    {
      semid = -1;
      numSem = 0;
    }

    inline bool       pOperator(int i, int flag = 0)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::pOperator(int i): i = " << i << endl;
      forthArgVal = semctl(semid, i, GETVAL);
      cerr << "nskaIPC::IPCSem::pOperator(int i): value: " << forthArgVal << endl;
#endif

      struct sembuf   semoperation;
      semoperation.sem_num = i;
      semoperation.sem_op = -1;
      semoperation.sem_flg = flag;
      if(semop(semid, &semoperation, 1)) return false;
      return true;
    }

    inline bool       vOperator(int i, int flag = 0)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCSem::vOperator(int i): i = " << i << endl;
      forthArgVal = semctl(semid, i, GETVAL);
      cerr << "nskaIPC::IPCSem::vOperator(int i): value: " << forthArgVal << endl;
#endif
      struct sembuf   semoperation;
      semoperation.sem_num = i;
      semoperation.sem_op =  1;
      semoperation.sem_flg = flag;
      if(semop(semid, &semoperation, 1)) return false;
      return true;
    }

    int semid;
    int numSem;
    int forthArgVal;
    //uint16_t forthArray[maxSemsPerID];
  };
#endif

  //! Class for creating of an IPC shared memory area.
  /*!
    The creation of an IPC shared memory area uses an IPC key.
    \n\n
    \author \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
  */
  class IPCShm : virtual public IPCKey
  {
  public:
    //! Default constructor.
    IPCShm(bool shm_check=false) : IPCKey() { init(shm_check); } /////////////////

    //! Constructor with file name and file creation flag.
    IPCShm(const char * fname, bool create, bool shm_check=false) : IPCKey(fname, create) { init(shm_check); } /////////////////

    //! Destructor
    virtual ~IPCShm(void) { cleanUp(); }

    //! Create shared memory with file name, file creation flag and size
    inline bool shmCreate(const char * fname, bool create, unsigned int s = 0)
    {
#if KALATTICEDEBUG
      cerr<<"shmCreate, fname="<<fname<<", create="<<create<<endl;
#endif
      keyCreate(fname, create);
      return shmCreate(s);
    }

    //! Create shared memory with size
    inline bool shmCreate(unsigned int s = 0)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::shmCreate(" << s << ")" << endl;
#endif
      if (noshm==true) { /////////////////
	shmsize=s;
	Address=(char *)malloc(s);
	return true;
      }
#ifdef osWin32
     //create the mapping names
     CreateMappingID(getName(), mapname);
     strcpy(attmapname, mapname);
     strcat(attmapname,"attached");

     //create file mapping for the file
     shmid = CreateFileMapping((HANDLE) 0xFFFFFFFF, NULL, PAGE_READWRITE, 0, s, mapname);

     if (GetLastError() == ERROR_ALREADY_EXISTS) {
       shmid = OpenFileMapping(FILE_MAP_ALL_ACCESS, 0, mapname);
       //map file exists -> read number of clients and add one
       attID = OpenFileMapping(FILE_MAP_ALL_ACCESS, 0, attmapname);
       attAddress = (char *)MapViewOfFile(attID, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, sizeof(char));
       *attAddress ++;
     }
     else {
       //create file mapping for number of clients
       attID = CreateFileMapping((HANDLE) 0xFFFFFFFF, NULL, PAGE_READWRITE, 0, sizeof(char), attmapname);
       attAddress = (char *)MapViewOfFile(attID, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, sizeof(char));
       *attAddress = 1;
     }

     if (shmid == NULL)
         throw InvalidIPC(("Shared memory not created. Bad ID.");

     shmInfo.shm_segsz = s;
     shmInfo.shm_nattch = *attAddress;

#else
#if !defined(NO_SHM)

     if(getKey() < 0)
       throw InvalidIPC("Shared memory %s not created. %s.", getName(), strerror(errno));

      shmid = shmget(getKey(), s, IPC_CREAT | 0666 | IPC_EXCL);
      if(shmid < 0)
        shmid = shmget(getKey(), 0, 0666);
      if(shmid < 0)
        throw InvalidIPC("Shared memory %s not created. %s.", getName(), strerror(errno));
#endif
#endif

      shmAttach();
      shmGetInfo();
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::shmCreate(" << s << ") shmid = " << shmid << endl;
#endif
      return true;
    }

    //! Get size of allocated shared memory block
    inline int calcSizeOfShm(void)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::calcSizeOfShm(void)" << endl;
#endif
      int ret = 0;

      if (noshm==true)  /////////////////
	return shmsize;

#ifdef osWin32
      CreateMappingID(getName(), mapname);
      HANDLE hMap = OpenFileMapping(FILE_MAP_ALL_ACCESS, 0, mapname);
      if (hMap != NULL)
         ret = Size();
#else
#if !defined(NO_SHM)

      if(getKey() < 0)
	throw InvalidIPC("Shared memory size not determined. Bad key.");
      int tmpid = shmget(getKey(), 0, 0666);
      if(tmpid >= 0)
        {
          kaShmDS tmpInfo; // df: changed from struct shmid_ds 27.06.2007

          shmctl(tmpid, IPC_STAT, (struct shmid_ds*)&tmpInfo);
          ret = (int)tmpInfo.shm_segsz;

          shmctl(tmpid, IPC_STAT, (struct shmid_ds*)&tmpInfo);
          if(tmpInfo.shm_nattch == 0) shmctl(tmpid, IPC_RMID, 0);
        }
#endif
#endif
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::calcSizeOfShm(void): ret = " << ret << endl;
#endif
      return ret;
    }

    //! return start address of shared memory block
    inline char *getAddress(void) { return Address; }

    //!Get number of attached clients
    inline int numAttach(void) const {
      if (noshm==true) /////////////////
	return Address ? 1 : 0;
#if !defined(NO_SHM)
      return shmInfo.shm_nattch;
#else
      return 0;
#endif
    }

  private:

    //! Attach to shared memory block
    inline void shmAttach(void)
    {
#ifdef osWin32
      Address = (char *) MapViewOfFile(shmid, FILE_MAP_ALL_ACCESS, 0, 0, size());
#else
#if !defined(NO_SHM)
      Address = (char *) shmat(shmid, 0, 0);

      if (Address==(char *)-1)
	throw InvalidIPC("Attachment to shared memory %s failed. %s.", getName(), strerror(errno));

#endif
#endif
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::shmAttach(void) Address = " << (long) Address << endl;
#endif
    }

    //!Check for shared memory info
    inline void shmGetInfo(void)
      {
#ifdef osWin32
	shmInfo.shm_nattch = *attAddress;
#else
#if !defined(NO_SHM)
	shmctl(shmid, IPC_STAT, (struct shmid_ds*)&shmInfo);
#endif
#endif
      }

    //! Detach from shared memory block
    inline void shmDetach(void)
      {
#ifdef osWin32
        *attAddress -= 1;
        shmInfo.shm_nattch = *attAddress;
        if (shmInfo.shm_nattch > 0) {
           Address = NULL;
           shmid = NULL;
           attID = NULL;
           attAddress = NULL;
        }
#else
#if !defined(NO_SHM)
	shmdt(getAddress());
#endif
#endif
    }

    //! Detach from shared memory block and deallocate if no further clients are attached
    inline void shmDestroy(void)
    {
#if KALATTICEDEBUG
      cerr << "nskaIPC::IPCShm::shmDestroy(void)" << endl;
#endif
#ifdef osWin32
     //check if number of clients are zero, if yes -> destroy
     if(shmInfo.shm_nattch==0) {
         //destroy
         if (Address) {
            UnmapViewOfFile (Address);
            Address = NULL;
         }

         if (shmid) {
            CloseHandle (shmid);
            shmid = NULL;
         }

         if (attAddress) {
            UnmapViewOfFile (attAddress);
            attAddress = NULL;
         }

         if (attID) {
            CloseHandle (attID);
            attID = NULL;
         }
     }
#else
#if !defined(NO_SHM)
      shmGetInfo();
      if(numAttach() == 0) shmctl(shmid, IPC_RMID, 0);
#endif
#endif
    }


    //! return size of shared memory block
#if !defined(NO_SHM)
    inline int size(void) const { return (int)shmInfo.shm_segsz; }
#endif


  protected:
    inline void       cleanUp(void)
    {
      if (noshm==true) { /////////////////
	if (Address) {
	  free(Address);
	  Address=NULL;
	}
      }
#if !defined(NO_SHM)
      else {
	shmDetach();
	shmDestroy();
      }
#endif
    }

  private:
    inline void       init(bool shm_check) /////////////////
    {
      Address = 0;
#ifdef osWin32
      attAddress = 0;
      shmid = NULL;
      attID = NULL;
#else
      shmid = -1;
#endif
#if !defined(NO_SHM)
      shmInfo.shm_nattch = 0;
      shmInfo.shm_segsz = 0;
#endif

	  noshm=false;
	  if( (shm_check==true) || (getenv("kaNoShm")) )
	  {noshm=true;}

#if defined(NO_SHM)
      noshm=true;
#endif

	  shmsize=0;
    }

#ifdef osWin32
    HANDLE            shmid;
    char              mapname[512];
    HANDLE            attID;
    char              *attAddress;
    char              attmapname[512];
#else
    int               shmid;
#endif
    unsigned int               shmsize;
#if !defined(NO_SHM)
    struct kaShmDS    shmInfo; // df: changed from struct shmid_ds 27.06.2007
#endif
    char *            Address;
	bool              noshm;
  };
}


#endif
