/*
 * File: kaMachineOS.h
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


#ifndef MACHINEOS_H
#define MACHINEOS_H

#ifndef KALATTICEDEBUG
# define KALATTICEDEBUG 0
#endif  // ifndef KALATTICEDEBUG

#ifdef WIN32
# define osWin32
# define _STANDARD_C_PLUS_PLUS
#endif  // ifdef WIN32

#if defined(linux) || defined(__linux) || defined(__linux__)
# define osLinux
# ifndef _STANDARD_C_PLUS_PLUS
#  define _STANDARD_C_PLUS_PLUS
# endif  // ifndef _STANDARD_C_PLUS_PLUS
# define __STDC_LIMIT_MACROS
# ifndef _ISOC99_SOURCE
#  define _ISOC99_SOURCE
# endif  // ifndef _ISOC99_SOURCE
# ifndef _SVID_SOURCE
#  define _SVID_SOURCE
# endif  // ifndef _SVID_SOURCE
# ifndef _BSD_SOURCE
#  define _BSD_SOURCE
# endif  // ifndef _BSD_SOURCE
#endif  // if defined(linux) || defined(__linux) || defined(__linux__)

#ifdef __DARWIN__
# ifndef __APPLE__
#  define __APPLE__
# endif  // ifndef __APPLE__
#endif  // ifdef __DARWIN__

#ifdef __APPLE__
# define osMac
# define _STANDARD_C_PLUS_PLUS
# ifdef __x86_64__
#  define i386 1
# endif  // ifdef __x86_64__
#endif  // ifdef __APPLE__

/* macro 'i386' is deprecated in newer c++ compilers */
#if !defined(i386) && (defined(__i386) || defined(__i386__) || defined(__x86_64__) )
# define i386 1
#endif  // if !defined(i386) && (defined(__i386) || defined(__i386__) || defined(__x86_64__) )

#ifdef AIXV3
# define osAIX
# define _STANDARD_C_PLUS_PLUS
# undef hz
# ifndef MAXFLOAT
#  define MAXFLOAT    3.4e+38
# endif  // ifndef MAXFLOAT

#endif  // ifdef AIXV3

#ifdef __sparc
# define osSolaris
#endif  // ifdef __sparc

#ifdef __sgi
# define osIRIX
# ifndef _SGI_MP_SOURCE
#  define _SGI_MP_SOURCE
# endif  // ifndef _SGI_MP_SOURCE
#endif  // ifdef __sgi


#if !defined(osIRIX) && !defined(osSolaris) && !defined(osAIX) && !defined(osMac) && !defined(osLinux) && \
  !defined(osWin32) && !defined(osAIX)
# error  MachineType is not defined
#endif  // if !defined(osIRIX) && !defined(osSolaris) && !defined(osAIX) && !defined(osMac) && !defined(osLinux) &&
        // !defined(osWin32) && !defined(osAIX)

#if defined(osWin32) || (defined(osLinux) && defined(i386)) || (defined(osMac) && defined(i386))
# define KA_LITTLE_ENDIAN
#else  // df+ - otherwise ppc-architecture is broken
# undef LITTLE_ENDIAN
#endif  // if defined(osWin32) || (defined(osLinux) && defined(i386)) || (defined(osMac) && defined(i386))

#if (defined(osLinux) && !defined(i386))
# define NO_SHM
typedef __int64          int64_t;
typedef unsigned __int64 __uint64_t;
#endif  // if (defined(osLinux) && !defined(i386))

#ifdef osWin32
# define VC_EXTRALEAN
# include "StdAfx.h"
# include <winbase.h>
# include <winnt.h>
#endif  // ifdef osWin32

#define _BSD_COMPAT

#ifdef osMac
# define __STDC_LIMIT_MACROS

// this macros is used in stdint.h:
//      included in mach_types.h
//      included in pthread.h
// to define such things as INT64_MAX etc.
#endif  // ifdef osMac

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#ifndef osWin32

# if !defined(_STANDARD_C_PLUS_PLUS) && defined(osLinux)  // df+ osLinux condition 03.06.2003

template<class X> inline X min(const X &a, const X &b) {return a < b ? a : b;}

template<class X> inline X max(const X &a, const X &b) {return a > b ? a : b;}

# endif  // if !defined(_STANDARD_C_PLUS_PLUS) && defined(osLinux)

# ifdef __cplusplus
#  include <cmath>
#  undef HUGE
#  define HUGE HUGE_VAL
#  if __cplusplus >= 201103L
#   undef MAXFLOAT
#   define MAXFLOAT HUGE_VALF
#  endif  // if __cplusplus >= 201103L
# elif defined(osAIX)
#  define isnanf isnan
#  include </usr/include/math.h>
#  undef acos
#  undef asin
#  undef atan
#  undef atan2
#  undef cos
#  undef exp
#  undef fabs
#  undef log
#  undef log10
#  undef sin
#  undef sqrt
#  undef tan
# else  // ifdef __cplusplus
#  include <math.h>
# endif  // ifdef __cplusplus

# include <signal.h>
# include <sys/file.h>
# include <sys/ipc.h>
# if !defined(osMac)  // || (__GNUC__ == 4)
#  include <sys/msg.h>
# endif  // if !defined(osMac)
# include <libgen.h>
# include <sys/sem.h>
# include <unistd.h>
# ifndef NODIRENT
#  include <dirent.h>
# endif  // ifndef NODIRENT
#endif  // ifndef osWin32


#ifdef _STANDARD_C_PLUS_PLUS
# include <iostream>
# include <typeinfo>
# include <iomanip>
# include <fstream>
# include <string>
# include <map>
# include <vector>
using namespace std;
#else  // ifdef _STANDARD_C_PLUS_PLUS
# include <iostream.h>
# include <typeinfo.h>
# include <iomanip.h>
# include <fstream.h>
# include <map.h>
# include <vector.h>
#endif  // ifdef _STANDARD_C_PLUS_PLUS

#if !defined(osMac) && !defined(osAIX) && !defined(osLinux)
# include <sys/prctl.h>
#endif  // if !defined(osMac) && !defined(osAIX) && !defined(osLinux)

#define kaShmDS shmid_ds

#ifndef NOPTHREADS
# include <pthread.h>
# include <strings.h>
typedef pthread_t kaptid;
typedef void     *kacallfctn;

inline int kaExit() {
  return pthread_cancel(pthread_self());
}

inline kaptid kaSproc(kacallfctn (*fctn)(void *), void *args = NULL) {
  kaptid ptid;

  if (!pthread_create(&ptid, NULL, fctn, args))
    return ptid;
  else
    return (kaptid)-1;
}

inline int kaKill(kaptid pid, int signal) {
  int result = -1;

  switch (signal) {
    case SIGKILL:
      result = pthread_cancel(pid);
    default:
      cerr << "pthread kaKill doesn\'t support signal\n";
      break;
  }
  return result;
}

#else  // ifndef NOPTHREADS
# ifndef osWin32
typedef pid_t kaptid;
typedef void  kacallfctn;
# endif  // ifndef osWin32
#endif  // ifndef NOPTHREADS

#ifndef osWin32
int kaBlockProc(pid_t);
# if defined(osMac) && !defined(i386)
int kaUnblockProc(pid_t, int);
# endif  // if defined(osMac) && !defined(i386)
int kaUnblockProc(pid_t);
#endif  // ifndef osWin32
void kaLockFile(int, int LOCK_STYLE = LOCK_EX);
void kaUnlockFile(int);

#ifdef NO_SHM
struct shmid_ds {
  uint32_t shm_nattch;  // used only for shminfo
  int shm_segsz;       // size of segment in bytes
};

#endif  // ifdef NO_SHM

#ifdef osWin32

# include <io.h>
# include <float.h>
# include <direct.h>
# include <limits.h>
# include <sys/locking.h>
# include <basetsd.h>

# ifndef HUGE
#  define HUGE HUGE_VAL
# endif  // ifndef HUGE

# ifndef M_PI
#  define M_PI        3.14159265358979323846264338327950288419716939937511
# endif  // ifndef M_PI
# ifndef MAXFLOAT
#  define MAXFLOAT    FLT_MAX
# endif  // ifndef MAXFLOAT

# define vsnprintf _vsnprintf

typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
typedef __int8           int8_t;
typedef __int16          int16_t;
typedef __int32          int32_t;
typedef __int64          int64_t;
int strcasecmp(const char *string1, const char *string2);
int strncasecmp(const char *s1, const char *s2, size_t n);

typedef int32_t  key_t;
typedef uint32_t shmatt_t;

struct shmid_ds {
  shmatt_t shm_nattch;  /* used only for shminfo */
  int shm_segsz;       /* size of segment in bytes */
};
void CreateMappingID(const char *fname, char *cID);

# define LOCK_EX _LK_LOCK
# define LOCK_UN _LK_UNLCK
int   flock(int fd, int operation);
key_t kaFtok(const char *fname, int id);


#else  // osWin32

# if !defined(NO_SHM)
#  include <sys/shm.h>
# endif  // if !defined(NO_SHM)
# include <sys/wait.h>

# if defined(osLinux) || defined(osAIX)
union semun {
  int val;                   // <= value for SETVAL
  struct semid_ds *buf;      // <= buffer for IPC_STAT & IPC_SET
  unsigned short int *array;  // <= array for GETALL & SETALL
  struct seminfo *__buf;     // <= buffer for IPC_INFO
};
# endif  // if defined(osLinux) || defined(osAIX)


# ifdef osLinux
#  include <stdint.h>
#  include <values.h>
#  include <list>

inline int unordered(double x, double y) {
  //  extern int __builtin_isunordered(double, double);
  //  return __builtin_isunordered (x, y);
  return isunordered(x, y);
}

// df 30.06.2008: CLK_TCK is not CLOCKS_PER_SEC (kernel 2.6.24)
// moved definition here
#  ifndef CLK_TCK
#   define CLK_TCK sysconf(_SC_CLK_TCK)
#  endif  // CLK_TCK

#  ifndef SEM_R
#   define SEM_R 0400
#  endif  // SEM_R
#  ifndef SEM_A
#   define SEM_A 0200
#  endif  // SEM_A
# endif  // ifdef osLinux

# ifdef osIRIX
#  include <inttypes.h>
#  ifdef NOPTHREADS

inline int kaExit() {
  return kill(getpid(), SIGTERM);
}

inline kaptid kaSproc(kacallfctn (*fctn)(void *), void *args = NULL) {
  return sproc(fctn, PR_SADDR | PR_SFDS, args);
}

inline int kaKill(kaptid pid, int signal) {
  return kill(pid, signal);
}

#  endif  // ifdef NOPTHREADS
# endif  // ifdef osIRIX
key_t kaFtok(const char *fname, int id);

#endif  // osWin32

#if defined(osMac) || defined(osLinux) || defined(osAIX)

# ifndef NODIRENT
typedef dirent dirent_t;
# endif  // ifndef NODIRENT

#endif  // if defined(osMac) || defined(osLinux) || defined(osAIX)

#if defined(osMac) || defined(osLinux)
extern int mkdirp(const char *pathname, mode_t mode);

#endif  // osMac || osLinux

#ifdef osMac

# include <mach/mach.h>
# include <float.h> // declared DBL_MAX

# define MAXDOUBLE       DBL_MAX
/*
   inline double drand48(void)
   {
   return rand()/(double)RAND_MAX;
   }

   inline void srand48(long seedval)
   {
   srand(seedval);
   }
 */

# define      fpclass(x) ( (sizeof(x) == sizeof(double) ) ? \
                           __fpclassifyd(x) :               \
                           (sizeof(x) == sizeof(float) ) ?  \
                           __fpclassifyf(x) :               \
                           __fpclassify(x) )

# ifdef IBMXLC

inline int isnan(double x) {
  return __isnand(x);
}

# else  // ifdef IBMXLC
#  if (__GNUC__ < 4)

inline int unordered(double x, double y) {
  return __builtin_isunordered(x, y);
}

inline int isnan(double x) {
  return __isnand(x);
}

#  endif  // if (__GNUC__ < 4)
# endif  // ifdef IBMXLC

# define fcos cos
# define fsin sin
# define acosf acos
# define fsqrt sqrtf

# ifndef NOPTHREADS
extern int kaBlockProc(pthread_t pid);
extern int kaUnblockProc(pthread_t pid);
# endif  // NOPTHREADS

# ifdef __ppc64__

// df+ 27.06.2007: hack to enable shared memory on 64-bit systems. (TO TEST ON x86_64)
struct __ipc_perm_old {
  __uint16_t cuid;  /* Creator's user ID */
  __uint16_t cgid;  /* Creator's group ID */
  __uint16_t uid;  /* Owner's user ID */
  __uint16_t gid;  /* Owner's group ID */
  mode_t mode;     /* Read/Write permission */
  __uint16_t seq;  /* Reserved for internal use */
  key_t key;       /* Reserved for internal use */
};
struct __shmid_ds_old {
  struct __ipc_perm_old shm_perm;  /* [XSI] Operation permission value */
  size_t shm_segsz;    /* [XSI] Size of segment in bytes */
  pid_t  shm_lpid;     /* [XSI] PID of last shared memory op */
  pid_t  shm_cpid;     /* [XSI] PID of creator */
  short  shm_nattch;   /* [XSI] Number of current attaches */
  time_t shm_atime;    /* [XSI] Time of last shmat() */
  time_t shm_dtime;    /* [XSI] Time of last shmdt() */
  time_t shm_ctime;    /* [XSI] Time of last shmctl() change */
  void  *shm_internal; /* reserved for kernel use */
};
#  undef kaShmDS
#  define kaShmDS __shmid_ds_old
# endif  // __ppc64__
// apple's sem.h is missing one variant in X.2:
// #if !(__GNUC__ == 3 && __GNUC_MINOR__ > 2)
//
// inline int semctl(int semid, int semnum, int cmd)
// {
//      union semun arg;
//      arg.val = 0;
//      return semctl(semid, semnum, cmd, arg);
// }
// #endif

/*
 *      Here goes DF's implemenation of message queue mechanism,
 *      not implemented in Mach/Darwin. FIFO channels are used
 *      instead.
 *
 *      msg... functions are emulated. One instance of class
 *      msg_queue_mgr_t is used for it. It contains an indexed
 *      map of fifo channels, index denominates the msg queue key.
 */

// #if (__GNUC__ < 4)

# include <string>

# define MSG_FIFO_DIR "apple_msg_fifo"

const int32_t MESSAGE_BODY_LENGTH  = 256;
const int32_t MESSAGE_TYPE_DEFAULT = 0xDF24;

# define MsgSizeMax 1000

namespace msg_queues {
/*--------------------------------------message_t--------------------------------------*/
class MSGExcInvMessage {};

# if defined(__x86_64__) || defined(__ppc64__)
struct MsgBuf {
  int32_t msgtyp;
  char buf[MsgSizeMax];
};
# else  // if defined(__x86_64__) || defined(__ppc64__)
struct MsgBuf {
  int  msgtyp;
  char buf[MsgSizeMax];
};
# endif  // if defined(__x86_64__) || defined(__ppc64__)

typedef MsgBuf message_t;

/*------------------------------------end message_t------------------------------------*/

class msg_queue_mgr_t {
 public:
  /*    Constructor. Does nothing significant.  */
  msg_queue_mgr_t();

  /*    Destructor. Destroys the fifo channels map.     */
  ~msg_queue_mgr_t();

  /*    new_fifo_channel: if there's no channel with this id - creates one      */
  int new_fifo_channel(int channel_id, const char *filename, int perm = 0666);

  /*    delete_fifo_channel - safe procedure to erase existing fifo channel     */
  int delete_fifo_channel(int channel_id);

  /*    send_message and receive_message - interface routines to deal with fifo buffers */
  int send_message(int channel_id, const char *message, int flags);
  int receive_message(int channel_id, char *message, int flags);

  int get_error_code() {return m_error;}

  /*    Called from kaFtok.     */
  void set_new_filename(const char *filename);
  /************************************************************************************
  *                             Attn!!! This should be called before kaFtok!!!
  *                                         *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *
  *                                                                                                                   *
  *     Utility function to create a fifo file (and delete non-fifo file on its place)  *
  *     Used in BlackBoardServer and replacing original creation of a file for msgq id  *
  ************************************************************************************/
  static int create_fifo_file(const char *filename, int perm = 0666);

 protected:
  /*    check if there's a fifo channel with this id already allocated  */
  bool fifo_channel_exists(int channel_id);

  class bad_folder_exception {};
  class file_error_exception {};

  /*    Class for handling 1 fifo channel. Opens the fifo file in constructor,  */
  /*    closes it in destructor. <DF>
   */
  class fifo_channel_t {
   public:
    /*  Constructor. Creates fifo file and opens it for read & write.   */
    fifo_channel_t(string filename, int perm);

    /*  Destructor. Closes the fifo file if it was open.        */
    ~fifo_channel_t();

    /*Write to channel; emulates msgsnd() UNIX function. Data length is fixed and max!  */
    int put_to_channel(const char *data, int msg_flags);

    /*Read from channel; emulates msgrcv() UNIX function. Data length is fixed and max! */
    int take_from_channel(char *data, int msg_flags);

   protected:
    int m_fd;
    string m_filename;
  };

  /*************************************        end of class fifo_channel_t     *************************************/

  map<int, fifo_channel_t *> m_fifo_channels;
  int m_error;

  string m_last_filename;  // buffer changed from kaFtok. Last value is used for fifo buffer filename.
};  // class msg_queue_mgr_t

extern msg_queue_mgr_t msg_queue_mgr;

/*      returns 0 if alles OK   */
inline int msgsnd(int mq_id, const void *message, size_t msg_size, int flags) {
  if (message == NULL)
    return -1;

  char msg[sizeof(message_t)];
  memset(msg, 0, sizeof(message_t));

  // #if defined (__x86_64__) || defined (__ppc64__)
  //      memcpy(msg, message, msg_size + sizeof(long));
  // #else
  memcpy(msg, message, msg_size + sizeof(int32_t));

  // #endif

  return msg_queue_mgr.send_message(mq_id, msg, flags);
}

/*      If all OK, returns amount of received data, which is a fixed value      */
inline int msgrcv(int32_t mq_id, void *message, uint32_t msg_size, int32_t msg_type, int32_t flags) {
  char buf[sizeof(message_t)];

  memset(buf, 0, sizeof(message_t));

  int err = 0;
  if ((err = msg_queue_mgr.receive_message(mq_id, buf, flags)) != -1) {
    message_t *msg = (message_t *)buf;

    /*      while (msg -> msgtyp > msg_type && msg_type != 0)
          {
            msg_queue_mgr.send_message(mq_id, buf, flags);
            msg_queue_mgr.receive_message(mq_id, buf, flags);
          }*/
    while (true) {
      if (msg_type == 0)
        break;
      if ((msg_type > 0) && (msg->msgtyp == msg_type) )
        break;
      if ((msg_type < 0) && (msg->msgtyp <= -msg_type) )
        break;
      msg_queue_mgr.send_message(mq_id, buf, flags);
      msg_queue_mgr.receive_message(mq_id, buf, flags);
    }

    // #if defined (__x86_64__) || defined (__ppc64__)
    //      memcpy(message, buf, msg_size + sizeof(long));
    // #else
    memcpy(message, buf, msg_size + sizeof(int32_t));

    // #endif
    return msg_size;
  }
  return err;
}  // msgrcv

/*      always returns mq_id! (or -1 if fails)  */
inline int msgget(key_t mq_id, int perm) {
  return msg_queue_mgr.new_fifo_channel(mq_id, NULL, (perm & 0777));
}

inline int msgctl(int msqid, int cmd, ...) {
  if (cmd == IPC_RMID)
    return msg_queue_mgr.delete_fifo_channel(msqid);

  return -1;
}
}  // namespace msg_queues

using namespace msg_queues;

// #endif // __GNUC__ < 4

#endif  // ifdef osMac

#ifdef osAIX
# undef hz
#endif  // ifdef osAIX

#endif  // _MACHINEOS_H
