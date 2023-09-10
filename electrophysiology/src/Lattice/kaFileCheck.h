/*! \file kaFileCheck.h
   \brief Basic functions for file handling

   \author cw,fs,cw, IBT - Universität Karlsruhe (TH)
 */

#ifndef KAFILECHECK_H
#define KAFILECHECK_H

#include <kaMachineOS.h>

namespace nskaGlobal {
#ifdef osLinux
# ifdef __SVR4_I386_ABI_L1__
typedef struct timestruc timespec_t;
typedef long             blkcnt_t;
# else  // ifdef __SVR4_I386_ABI_L1__
typedef time_t timespec_t;
# endif  // ifdef __SVR4_I386_ABI_L1__
#endif  // ifdef osLinux

#ifdef osAIX
typedef time_t timespec_t;
#endif  // ifdef osAIX

#if defined(osMac)
typedef time_t timespec_t;
typedef long   blkcnt_t; /* blocks in a file */
#endif  // osMac

#ifdef osWin32
typedef unsigned short mode_t;
typedef unsigned short nlink_t;
typedef short          uid_t;
typedef short          gid_t;
typedef short          blkcnt_t;
typedef long           timespec_t;
# define kaStat _stat
#else  // ifdef osWin32
# define kaStat stat
#endif  // ifdef osWin32

//! Class for checking of file existance and attributes
/*!
   \author \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

class FileCheck {
 public:
  //! Default constructor
  FileCheck(void) {
    fc_error  = 0;
    fc_access = false;
    fc_exist  = false;

    fc.st_dev   = (dev_t)0;
    fc.st_ino   = (ino_t)0;
    fc.st_mode  = (mode_t)0;
    fc.st_nlink = (nlink_t)0;
    fc.st_uid   = (uid_t)0;
    fc.st_gid   = (gid_t)0;
    fc.st_rdev  = (dev_t)0;
    fc.st_size  = (off_t)0;
#ifndef osWin32
    fc.st_blksize = (long)0;
    fc.st_blocks  = (blkcnt_t)0;
#endif  // ifndef osWin32
#ifdef osLinux
# ifdef __SVR4_I386_ABI_L1__
    strcpy(fc.st_fstype, "");
# else  // ifdef __SVR4_I386_ABI_L1__

    // do nothing: st_fstype does not exist !
# endif  // ifdef __SVR4_I386_ABI_L1__
#else  // ifdef osLinux
# if !defined(osAIX) && !defined(osMac)
#  ifndef osWin32
    strcpy(fc.st_fstype, "");
#  endif  // ifndef osWin32
# endif  // if !defined(osAIX) && !defined(osMac)
#endif  // ifdef osLinux
  }

  //! Get attributes of file
  void Check(const char *filename) {
    if (!filename) {
      fc_exist  = false;
      fc_access = false;
      fc_error  = EFAULT;
      return;
    }
    if (::kaStat(filename, &fc) == -1) {
      fc_error = errno;
      switch (fc_error) {
        case EACCES:
        case EINTR:
#ifndef osWin32
        case ETIMEDOUT:
        case ELOOP:
# ifndef osMac
        case EMULTIHOP:
        case ENOLINK:
# endif  // ifndef osMac
        case EOVERFLOW:
#endif  // ifndef osWin32
          fc_access = false;
          fc_exist  = true;
          break;
        case EFAULT:
        case ENAMETOOLONG:
        case ENOTDIR:
          fc_access = false;
          fc_exist  = false;
          break;
        case ENOENT:
          fc_access = true;
          fc_exist  = false;
          break;
      }  // switch
    } else {
      fc_error  = 0;
      fc_access = true;
      fc_exist  = true;
    }
  }  // Check

  //! Get errors of file access.
  int Error(void) {return fc_error;}

  //! Returns true, if file exists, otherwise false.
  bool Exist(void) {return fc_exist;}

  //! Returns true, if file can be accessed, otherwise false.
  bool Access(void) {return fc_access;}

#ifndef osWin32

  bool IsBlockSpecialFile(void) {return S_ISBLK(fc.st_mode);}

  bool IsCharSpecialFile(void) {return S_ISCHR(fc.st_mode);}

  bool IsDirectory(void) {return S_ISDIR(fc.st_mode);}

  bool IsPipe(void) {return S_ISFIFO(fc.st_mode);}

  bool IsRegularFile(void) {return S_ISREG(fc.st_mode);}

  bool IsSymbolicLink(void) {return S_ISLNK(fc.st_mode);}

  bool IsSocket(void) {return S_ISSOCK(fc.st_mode);}

#else  // ifndef osWin32

  bool IsBlockSpecialFile(void) {return false;}

  bool IsCharSpecialFile(void) {return fc.st_mode == _S_IFCHR ? true : false;}

  bool IsDirectory(void) {return fc.st_mode == _S_IFDIR ? true : false;}

  bool IsPipe(void) {return false;}

  bool IsRegularFile(void) {return fc.st_mode == _S_IFREG ? true : false;}

  bool IsSymbolicLink(void) {return false;}

  bool IsSocket(void) {return false;}

#endif  // ifndef osWin32

  dev_t Device(void) {return fc.st_dev;}

  ino_t Inode(void) {return fc.st_ino;}

  mode_t Mode(void) {return fc.st_mode;}

  nlink_t NumberOfLinks(void) {return fc.st_nlink;}

  uid_t UserID(void) {return fc.st_uid;}

  gid_t GroupID(void) {return fc.st_gid;}

  dev_t SpecialDevice(void) {return fc.st_rdev;}

  off_t Size(void) {return fc.st_size;}

#ifdef osLinux
# ifdef __SVR4_I386_ABI_L1__
  timespec_t & Lastaccess(void) {return fc.st_atim;}

  timespec_t & LastModified(void) {return fc.st_mtim;}

  timespec_t & LastStatusChange(void) {return fc.st_ctim;}

# else  // ifdef __SVR4_I386_ABI_L1__
  timespec_t & Lastaccess(void) {return fc.st_atime;}

  timespec_t & LastModified(void) {return fc.st_mtime;}

  timespec_t & LastStatusChange(void) {return fc.st_ctime;}

# endif  // ifdef __SVR4_I386_ABI_L1__
#else  // ifdef osLinux
# if defined(osAIX) || defined(osMac) || (osWin32)
  timespec_t & Lastaccess(void) {return fc.st_atime;}

  timespec_t & LastModified(void) {return fc.st_mtime;}

  timespec_t & LastStatusChange(void) {return fc.st_ctime;}

# else  // if defined(osAIX) || defined(osMac) || (osWin32)
  timespec_t & Lastaccess(void) {return fc.st_atim;}

  timespec_t & LastModified(void) {return fc.st_mtim;}

  timespec_t & LastStatusChange(void) {return fc.st_ctim;}

# endif  // if defined(osAIX) || defined(osMac) || (osWin32)
#endif  // ifdef osLinux

#ifndef osWin32

  int BlockSize(void) {return fc.st_blksize;}

  blkcnt_t BlockCount(void) {return (blkcnt_t)fc.st_blocks;}

#else  // ifndef osWin32

  int BlockSize(void) {return 512;}

  blkcnt_t BlockCount(void) {return fc.st_size / 512;}

#endif  // ifndef osWin32
  struct kaStat & GetStatus(void) {return fc;}

 private:
  struct kaStat fc;
  int fc_error;
  bool fc_access;
  bool fc_exist;
};  // class FileCheck
}  // namespace nskaGlobal

#endif  // ifndef KAFILECHECK_H
