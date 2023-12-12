/*
 * File: kaMachineOS.cpp
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


#include "kaMachineOS.h"
#include "kaExceptions.h"

#ifdef osMac

// #if (__GNUC__ < 4)
msg_queues::msg_queue_mgr_t msg_queues::msg_queue_mgr;  // global message queue manager; should be created before
                                                        // kasem!!!
                                                        // #endif
#endif  // ifdef osMac

#ifndef osWin32

// #if !defined(osMac) || i386
# include "kaIPC.h"

# define NUM_SEMS 100

class kaBlockSem : public nskaIPC::IPCKey {
public:
    int semid;
    
    kaBlockSem() {
        semid = -1;
    }
    
    int semCreate(const char *fname) {
        keyCreate(fname, true);
        unlockReadFD();
        closeReadFD();
        
        semid = semget(getKey(), 1, IPC_CREAT | 0666);
        if (semid < 0)
            throw kaBaseException("Cannot get semaphore ID for file %s (%s errno %d)", fname, strerror(errno), errno);
        return semid;
    }
};

class kaBlockSemArray {
    kaBlockSem k[NUM_SEMS];
    
public:
    kaBlockSemArray() {
# if KALATTICEDEBUG == 2
        fprintf(stderr, "kaSemaphore::kaSemaphore\n");
# endif  // if KALATTICEDEBUG == 2
    }
    
    ~kaBlockSemArray() {
# if KALATTICEDEBUG == 2
        fprintf(stderr, "kaSemaphore::~kaSemaphore\n");
# endif  // if KALATTICEDEBUG == 2
        for (int i = 0; i < NUM_SEMS; i++) {
            semctl(k[i].semid, 0, IPC_RMID);
            remove(k[i].getName());
        }
    }
    
    int Init(long pid) {
# if KALATTICEDEBUG
        fprintf(stderr, "kaSemaphore::Init %ld\n", pid);
# endif  // if KALATTICEDEBUG
        char semname[256];
        sprintf(semname, "/tmp/kaSemaphore%ld", pid);
        
        int i;
        for (i = 0; i < NUM_SEMS; i++)
            if (k[i].getName())
                if (!strcmp(k[i].getName(), semname))
                    return k[i].semid;
        
        for (i = 0; i < NUM_SEMS; i++)
            if (k[i].semid == -1) {
                break;
            } else {
                struct sembuf semoperation;
                semoperation.sem_num = 0;
                semoperation.sem_op  = 0;
                semoperation.sem_flg = 0;
                int rc = semop(k[i].semid, &semoperation, 1);
                
                //      printf("del %d %d %d\n", i, rc, semID[i]);
                if (rc == -1)
                    break;
            }
        if (i == NUM_SEMS)
            throw kaBaseException("Not enough semaphores for block/unblock");
        
        return k[i].semCreate(semname);
    }  // Init
    
    int Block(long pid) {
# if KALATTICEDEBUG
        fprintf(stderr, "kaSemaphore::Block %ld\n", pid);
# endif  // if KALATTICEDEBUG
        int semid = Init(pid);
        
        struct sembuf semoperation;
        semoperation.sem_num = 0;
        semoperation.sem_op  = -1;
        semoperation.sem_flg = 0;
        int rc = semop(semid, &semoperation, 1);
        
# if KALATTICEDEBUG
        fprintf(stderr, "kaSemaphore::Block rc %ld\n", rc);
# endif  // if KALATTICEDEBUG
        
        return rc;
    }
    
    int Unblock(long pid) {
# if KALATTICEDEBUG
        fprintf(stderr, "kaSemaphore::Unblock %ld\n", pid);
# endif  // if KALATTICEDEBUG
        int semid = Init(pid);
        
        struct sembuf semoperation;
        semoperation.sem_num = 0;
        semoperation.sem_op  = 1;
        semoperation.sem_flg = 0;
        int rc = semop(semid, &semoperation, 1);
# if KALATTICEDEBUG
        fprintf(stderr, "kaSemaphore::Unblock rc %ld\n", rc);
# endif  // if KALATTICEDEBUG
        
        return rc;
    }
};  // class kaBlockSemArray

static kaBlockSemArray kasem;

int kaBlockProc(pid_t p) {
    return kasem.Block(p);
}

int kaUnblockProc(pid_t p) {
    usleep(100000);
    return kasem.Unblock(p);
}

// #endif
#endif  // ifndef osWin32

void kaLockFile(int fd, int LOCK_STYLE) {
#ifdef osWin32
    locking(fd, _LK_LOCK, Size());
#else  // ifdef osWin32
# if defined(osSolaris) || defined(osAIX)
    lseek(fd, 0, SEEK_SET);
    lockf(fd, F_TLOCK, 1);
# else  // if defined(osSolaris) || defined(osAIX)
#  if KADEBUG
    if (LOCK_STYLE && LOCK_SH)
        cerr << "kaLockFile: Shared locking for fd " << fd << "...\n";
#  endif  // if KADEBUG
    flock(fd, LOCK_STYLE);
# endif  // if defined(osSolaris) || defined(osAIX)
#endif  // ifdef osWin32
}

void kaUnlockFile(int fd) {
#ifdef osWin32
    locking(fd, _LK_UNLCK, Size());
#else  // ifdef osWin32
# if defined(osSolaris) || defined(osAIX)
    lseek(fd, 0, SEEK_SET);
    lockf(fd, F_ULOCK, 1);
# else  // if defined(osSolaris) || defined(osAIX)
    flock(fd, LOCK_UN);
# endif  // if defined(osSolaris) || defined(osAIX)
#endif  // ifdef osWin32
}

#ifdef osWin32

int strcasecmp(const char *string1, const char *string2) {
    return _stricmp(string1, string2);
}

int strncasecmp(const char *s1, const char *s2, size_t n) {
    return _strnicmp(s1, s2, n);
}

void CreateMappingID(const char *fname, char *cID) {
    char full[_MAX_PATH];
    
    _fullpath(full, fname, _MAX_PATH);
    
    //  fprintf(stdout,"%s\n",full);
    
    int len        = strlen(full);
    const char *pf = &full[0];
    char *pID      = &cID[0];
    
    for (; *pf; pf++)
        if ((*pf != '\\') && (*pf != ':'))
            *pID++ = *pf;
    
    *pID = '\0';
}

int flock(int fd, int operation) {
    return locking(fd, operation, 1);
}

key_t kaFtok(const char *fname, int id) {
    return 1;
}

#else  // ifdef osWin32

key_t kaFtok(const char *fname, int id) {
# ifndef NO_SHM
    key_t k = -1;
    struct stat s;
    if (stat(fname, &s) == 0) {
        key_t d;
        d  = ((0x000000FF & s.st_dev) << (sizeof(key_t)-1)*8);
        d |= ((0x0000FF00 & s.st_dev) << (sizeof(key_t)-3)*8);
        d |= ((0x00FF0000 & s.st_dev) >> (sizeof(key_t)-3)*8);
        d |= ((0xFF000000 & s.st_dev) >> (sizeof(key_t)-1)*8);
        ino_t &i = s.st_ino;
        k = (i & ~d)  | (~i & d);
        k = (k & ~id) | (~id & k);
        if (k < 0)
            k = -k;
#  ifdef osMac
        
        // #if (__GNUC__ < 4)
        msg_queue_mgr.set_new_filename(fname);
        
        // #endif
#  endif  // ifdef osMac
    }
    return k;
    
# else  // ifndef NO_SHM
    return 1;
    
# endif  // ifndef NO_SHM
}  // kaFtok

#endif  // ifdef osWin32

#if defined(osMac) || defined(osLinux)

/*
 implementation of mkdirp taken from
 http://communicator.sourceforge.net/sites/MITRE/distributions/GalaxyCommunicator/src/libGalaxy/util/mkdirp.c
 */
int mkdirp(const char *pathname, mode_t mode) {
    size_t len;
    int i = 0;
    
    len = strlen(pathname);
    
    char *dirpath = (char *)malloc(len+1);
    
    if (0 == dirpath) {
        return -1;
    }
    
    strcpy(dirpath, pathname);
    
    if (dirpath[len-1] == '/') {
        len--;
        dirpath[len] = '\0';
    }
    for (i = 1; i <= len; i++) {
        if (('/' == dirpath[i]) ||
            ('\0' == dirpath[i]) ) {
            dirpath[i] = '\0';
            
        RETRY:
            
            
            if (0 != mkdir(dirpath, mode)) {
                switch (errno) {
                    case EEXIST:
                        break;
                        
                    case EINTR:
                        goto RETRY;
                        
                    default:
                        free(dirpath);
                        return -1;
                }
            }
            dirpath[i] = '/';
        }
    }
    free(dirpath);
    return 0;
}  // mkdirp

#endif  // osMac || osLinux

#ifdef osMac

// #if __ppc__ || __ppc64__

// //DF's prctl implementation - only for PPC architecture

// int kaBlockProc(pid_t pid)
// {
//   mach_port_t task_id = 0;
//   if (task_for_pid(mach_task_self(), pid, &task_id) == KERN_SUCCESS)
//     {
//       if (task_suspend(task_id) == KERN_SUCCESS)
//      return 0;
//     }

//   return -1;
// }

// int kaUnblockProc(pid_t pid)
// {
//   return kaUnblockProc(pid, 100000);
// }

// int kaUnblockProc(pid_t pid, int uwait)
// {
//   mach_port_t task_id = 0;
//   kern_return_t res = KERN_SUCCESS;
//   usleep(uwait);             // !!! Schrecklich, but what to do?
//   if ((res = task_for_pid(mach_task_self(), pid, &task_id)) == KERN_SUCCESS)
//   {
//     if ((res = task_resume(task_id)) == KERN_SUCCESS)
//     {
//       return 0;
//     }
//     std::cerr << "Error unblocking process id = " << pid;
//     std::cerr << "; error code " << res << '\n';
//     return -1;
//   }
//      std::cerr << "Error getting task_id for process id = " << pid;
//      std::cerr << "; error code " << res << '\n';
//   return -1;
// }

// #endif // __ppc__

# ifndef NOPTHREADS

int kaBlockProc(pthread_t pid) {
    mach_port_t mthread = pthread_mach_thread_np(pid);
    
    return thread_suspend(mthread);
}

int kaUnblockProc(pthread_t pid) {
    mach_port_t mthread = pthread_mach_thread_np(pid);
    
    return thread_resume(mthread);
}

# endif  // NOPTHREADS
/*
 *      Here goes DF's implemenation of message queue mechanism,
 *      not implemented in Mach/Darwin. FIFO channels are used
 *      instead.
 *
 *      msg... functions are emulated. One instance of class
 *      msg_queue_mgr_t is used for it. It contains an indexed
 *      map of fifo channels, index denominates the msg queue key.
 */

/*********************  class msg_queue_mgr_t   ***************************/

// #if (__GNUC__ < 4)

/*      Constructor. Does nothing significant.  */
msg_queues::msg_queue_mgr_t::msg_queue_mgr_t() {
    m_error         = 0;
    m_last_filename = "msg_default";
}

/*      Destructor. Destroys the fifo channels map.     */
msg_queues::msg_queue_mgr_t::~msg_queue_mgr_t() {
    typedef map<int, fifo_channel_t *>::iterator I;
    
    for (I iter = m_fifo_channels.begin(); iter != m_fifo_channels.end(); iter++) {
        if (iter->second != NULL)
            delete iter->second;
        iter->second = NULL;
    }
    
    m_fifo_channels.erase(m_fifo_channels.begin(), m_fifo_channels.end());
}

/*      new_fifo_channel: if there's no channel with this id - creates one      */

int msg_queues::msg_queue_mgr_t::new_fifo_channel(int channel_id, const char *filename, int perm) {
    if (!fifo_channel_exists(channel_id)) {
        try {
            if ((filename != NULL) && (::strlen(filename) != 0) ) {
                m_last_filename = filename;
            }
            m_fifo_channels[channel_id] = new fifo_channel_t(m_last_filename, perm);
        } catch (...) {
            m_error = errno;
            delete_fifo_channel(channel_id);
            return -1;
        }
    }
    return channel_id;
}

/*      delete_fifo_channel - safe procedure to erase existing fifo channel     */

int msg_queues::msg_queue_mgr_t::delete_fifo_channel(int channel_id) {
    if (m_fifo_channels.find(channel_id) != m_fifo_channels.end()) {
        if (fifo_channel_exists(channel_id))
            delete m_fifo_channels[channel_id];
        
        m_fifo_channels.erase(channel_id);
        
        return 0;
    } else {return -1;}
}

int msg_queues::msg_queue_mgr_t::send_message(int channel_id, const char *message, int flags) {
    if (fifo_channel_exists(channel_id) and message != NULL) {
        return m_fifo_channels[channel_id]->put_to_channel(message, flags);
        
        // flags are transformed from IPC_.. to O_.. inside put_to_channel.
    }
    
    return -1;
}

int msg_queues::msg_queue_mgr_t::receive_message(int channel_id, char *message, int flags) {
    if (fifo_channel_exists(channel_id) and message != NULL) {
        return m_fifo_channels[channel_id]->take_from_channel(message, flags);
        
        // flags are transformed from IPC_.. to O_.. inside take_from_channel.
    }
    
    return -1;
}

/*      Called from kaFtok.     */

void msg_queues::msg_queue_mgr_t::set_new_filename(const char *filename) {
    if (filename == NULL)
        m_last_filename = "msg_default";
    else if (strlen(filename) <= 0)
        m_last_filename = "msg_default";
    else
        m_last_filename = filename;
    m_last_filename += '\0';
}

/************************************************************************************
 *                               Attn!!! This should be called before kaFtok!!!
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
 *                                                                                                                     *
 *       Utility function to create a fifo file (and delete non-fifo file on its place)  *
 *       Used in BlackBoardServer and replacing original creation of a file for msgq id  *
 ************************************************************************************/
int msg_queues::msg_queue_mgr_t::create_fifo_file(const char *filename, int perm) {
    /* check file name */
    if (filename == NULL)
        return -1;
    
    if (strlen(filename) == 0)
        return -1;
    
    /* check if a non-fifo file with the same name is present*/
    bool fifo_exists = false;
    
    struct stat filestat;
    
    if (stat(filename, &filestat) != -1) {                  // file present
        if ((filestat.st_mode & 0170000) != S_IFIFO) {        // not a fifo
            if ((unlink(filename) == -1) and (errno != ENOENT)) // erase not-fifo file
                return -1;
        } else {fifo_exists = true;}
    }
    
    if (!fifo_exists) {
        if (mkfifo(filename, perm))
            return -1;
    }
    
    return 0;
} // msg_queues::msg_queue_mgr_t

/*      check if there's a fifo channel with this id already allocated  */

bool msg_queues::msg_queue_mgr_t::fifo_channel_exists(int channel_id) {
    if (m_fifo_channels.find(channel_id) == m_fifo_channels.end())
        return false;
    
    if (m_fifo_channels[channel_id] == NULL)
        return false;
    
    return true;
}

/*********************************      class fifo_channel_t    ***********************************/

/*      Constructor. Creates fifo file and opens it for read & write.   */
msg_queues::msg_queue_mgr_t::fifo_channel_t::fifo_channel_t(string filename, int perm) {
    int err = 0;
    
    //    m_filename = dirname(filename.c_str());
    //    if (m_filename.size() > 0)
    //            m_filename += '/';
    //    m_filename += MSG_FIFO_DIR;
    //
    //            /*      Create a folder for a FIFO      */
    //    err = mkdir(m_filename.c_str(), 0777);                  //u_mask of the process
    //    if (err == -1 and errno != EEXIST)
    //            throw bad_folder_exception();
    //
    //    m_filename += '/';
    //    m_filename += basename(filename.c_str());
    m_filename = filename;
    
    /*    Create FIFO file        */
    err = mkfifo(m_filename.c_str(), perm);
    if (err == -1 and errno != EEXIST)
        throw file_error_exception();
    
    m_fd = open(m_filename.c_str(), O_RDWR);
    if (m_fd <= 0)
        throw file_error_exception();
}

/*      Destructor. Closes the fifo file if it was open.        */
msg_queues::msg_queue_mgr_t::fifo_channel_t::~fifo_channel_t() {
    if (m_fd > 0)
        close(m_fd);
}

/*Write to channel; emulates msgsnd() UNIX function. Data length is fixed and max!      */

int msg_queues::msg_queue_mgr_t::fifo_channel_t::put_to_channel(const char *data, int msg_flags) {
    ssize_t err;
    
    if (m_fd == 0) {
        int flags = O_RDWR;
        if (msg_flags & IPC_NOWAIT)
            flags |= O_NONBLOCK;
        if ((m_fd = open(m_filename.c_str(), flags)) == -1)
            return -1;
    }
    
    if (data == NULL)
        return -1;
    
    if ((err = write(m_fd, data, sizeof(message_t))) == -1)
        return -1;
    
    return /*err*/ 0;
}

/*Read from channel; emulates msgrcv() UNIX function. Data length is fixed and max!     */

int msg_queues::msg_queue_mgr_t::fifo_channel_t::take_from_channel(char *data, int msg_flags) {
    int err;
    bool blocking = (msg_flags & IPC_NOWAIT) == 0;
    
    if (m_fd == 0) {
        int flags = O_RDWR;
        if (msg_flags & IPC_NOWAIT)
            flags |= O_NONBLOCK;
        if ((m_fd = open(m_filename.c_str(), flags)) == -1)
            return -1;
    }
    
    if (data == NULL)
        return -1;
    
    if (blocking) {
        err = 0;
        while (err == 0) {
            err = (int)read(m_fd, data, sizeof(message_t));
        }
    } else if ((err = (int)read(m_fd, data, sizeof(message_t))) <= 0) {
        return -1;
    }
    
    return err;
}  // msg_queues::msg_queue_mgr_t::fifo_channel_t::take_from_channel

/*************************************  end of class fifo_channel_t     *************************************/

// #endif  //if (__GNUC__ < 4)
#endif  // ifdef osMac
