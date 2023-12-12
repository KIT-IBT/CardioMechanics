/*
 * File: PETScLSE.h
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


#ifndef PETSCLSE_H
#define PETSCLSE_H

#include <kaBasicIO.h>
#include <kaSharedMatrixN.h>
#include <kaLatticeHeader.h>
#include <petsc.h>
#include <string>
#include <iostream>

using namespace nskaGlobal;

#undef CHKERRQ
#define CHKERRQ(n) do {                                                                                                 \
    if (PetscUnlikely(n)) {                                                                                             \
      PetscError(PETSC_COMM_SELF, __LINE__, PETSC_FUNCTION_NAME, __FILE__, (PetscErrorCode)n, PETSC_ERROR_REPEAT, " "); \
      throw kaBaseException("PETSc error");                                                                             \
    }                                                                                                                   \
} while (0)
#undef SETERRQ1
#define SETERRQ1(n, s, ...) do {                                                                                   \
    PetscError(PETSC_COMM_SELF, __LINE__, PETSC_FUNCTION_NAME, __FILE__, (PetscErrorCode)n, PETSC_ERROR_REPEAT, s, \
               __VA_ARGS__);                                                                                       \
} while (0)

#if PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7)

inline PetscErrorCode PetscOptionsHasName(const char a[], const char b[], PetscBool *c) {
  return PetscOptionsHasName(PETSC_NULL, a, b, c);
}

inline PetscErrorCode PetscOptionsGetString(const char a[], const char b[], char c[], size_t d, PetscBool *e) {
  return PetscOptionsGetString(PETSC_NULL, a, b, c, d, e);
}

#endif  // if PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7)

#ifdef PETSC_USE_64BIT_INDICES
# define PETSCINT_FORMAT "%lld"
#else  // ifdef PETSC_USE_64BIT_INDICES
# define PETSCINT_FORMAT "%d"
#endif  // ifdef PETSC_USE_64BIT_INDICES


typedef enum {
  ot_bin,
  ot_ascii,
  ot_matlab
} IOType;


class BoundaryStorage {
 public:
  PetscInt index;
  PetscScalar value;
};

//! Get Header Filename From Viewer
inline std::string GetHeaderFilename(const char *file) {
  return std::string(file) + ".header";
}

//! Load Lattice Info
inline void LoadLatInfo(const char *file, kaSharedMatrixN<double> &m, IndexType &xLattice, IndexType &yLattice,
                        IndexType &zLattice) {
#if KADEBUG
  fprintf(stderr, "LoadLatInfo()\n");
#endif  // if KADEBUG

  FILE *fp;
  fp = fopen(GetHeaderFilename(file).c_str(), "r");

  if (!fp)
    throw kaBaseException("Reading of file .header failed");

  kaRead(&xLattice, 1, fp);
  kaRead(&yLattice, 1, fp);
  kaRead(&zLattice, 1, fp);
  m.Restore(fp);
  if (!feof(fp))
    kaBaseException(".header not completely read (maybe not a vector originating from kaLattice?)!");
  fclose(fp);

#if KADEBUG == 2
  fprintf(stderr, "LoadLatInfo finished\n");
#endif  // if KADEBUG == 2
}

//! Load Grid Info
inline void LoadGridInfo(const char *file, PetscInt &nPoints, PetscInt &nCells, std::string &origFilename) {
#if KADEBUG
  fprintf(stderr, "LoadGridInfo()\n");
#endif  // if KADEBUG

  FILE *fp;
  fp = fopen(GetHeaderFilename(file).c_str(), "r");

  if (!fp)
    kaBaseException("Reading of file .header failed");
  kaRead(&nPoints, 1, fp);
  kaRead(&nCells, 1, fp);
  size_t nameLength;
  kaRead(&nameLength, 1, fp);
  char *filename = new char[nameLength+1];
  kaRead(filename, sizeof(char), nameLength, fp);
  filename[nameLength] = '\0';
  origFilename         = filename;
  delete[] filename;
  fclose(fp);

#if KADEBUG == 2
  fprintf(stderr, "LoadGridInfo finished\n");
#endif  // if KADEBUG == 2
}

//! Save Lattice Info
inline void SaveLatInfo(const char *file, const kaSharedMatrixN<double> &m, const IndexType xLattice,
                        const IndexType yLattice, const IndexType zLattice) {
#if KADEBUG
  fprintf(stderr, "SaveLatInfo()\n");
#endif  // if KADEBUG

  int mpirank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank);
  if (mpirank) {
    return;
  }

  FILE *fp;
  fp = fopen(GetHeaderFilename(file).c_str(), "w");

  if (fp) {
    kaWrite(&xLattice, 1, fp);
    kaWrite(&yLattice, 1, fp);
    kaWrite(&zLattice, 1, fp);
    m.Save(fp);
    fclose(fp);
  }

#if KADEBUG == 2
  fprintf(stderr, "SaveLatInfo finished\n");
#endif  // if KADEBUG == 2
}  // SaveLatInfo

//! Save Grid Info
inline void SaveGridInfo(const char *file, const PetscInt nPoints, const PetscInt nCells) {
#if KADEBUG
  fprintf(stderr, "SaveGridInfo()\n");
#endif  // if KADEBUG

  FILE *fp;
  fp = fopen(GetHeaderFilename(file).c_str(), "w");
  if (fp) {
    kaWrite(&nPoints, 1, fp);
    kaWrite(&nCells, 1, fp);
    size_t nameLength = strlen(file);
    kaWrite(&nameLength, 1, fp);
    kaWrite(file, sizeof(char), nameLength, fp);
    fclose(fp);
  }

#if KADEBUG == 2
  fprintf(stderr, "SaveGridInfo finished\n");
#endif  // if KADEBUG == 2
}

//! Load vector of specified format
inline void LoadVec(const char *name, Vec &v, IOType ot) {
#if KADEBUG
  fprintf(stderr, "LoadVec('%s')\n", name);
#endif  // if KADEBUG

  PetscErrorCode ierr;
  PetscViewer fd;  /* viewer */

  switch (ot) {
    case ot_bin:
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_READ, &fd); CHKERRQ(ierr);
      break;

    case ot_ascii:
    case ot_matlab:
      throw kaBaseException("Format of input vector %s is not supported", name);
      break;

    default:
      throw kaBaseException("Format of input vector %s is missing", name);
      break;
  }
  VecCreate(PETSC_COMM_WORLD, &v);
  ierr = VecLoad(v, fd); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

#if KADEBUG == 2
  fprintf(stderr, "PETScLSE::LoadVec finished\n");
#endif  // if KADEBUG == 2
}  // LoadVec

//! Save vector in specified format
inline void SaveVec(const char *name, const Vec &v, IOType ot) {
#if KADEBUG
  fprintf(stderr, "SaveVec('%s')\n", name);
#endif  // if KADEBUG

  PetscErrorCode ierr;
  PetscViewer fd;  /* viewer */

  switch (ot) {
    case ot_bin:
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
      break;

    case ot_ascii:
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
      break;

    case ot_matlab:
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB);
      break;

    default:
      kaBaseException("Format of output matrix %s is missing", name);
      break;
  }
  ierr = VecView(v, fd); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

#if KADEBUG == 2
  fprintf(stderr, "SaveVec finished\n");
#endif  // if KADEBUG == 2
}  // SaveVec

inline void SaveVecMatInfo(const char *name, const Vec &v, IOType ot, const kaSharedMatrixN<double> &m,
                           const IndexType xLattice, const IndexType yLattice, const IndexType zLattice) {
#if KADEBUG
  fprintf(stderr, "SaveVec('%s')\n", name);
#endif  // if KADEBUG

  PetscErrorCode ierr;
  PetscViewer fd;  /* viewer */

  switch (ot) {
    case ot_bin:
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
      SaveLatInfo(name, m, xLattice, yLattice, zLattice);
      break;

    case ot_ascii:
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
      break;

    case ot_matlab:
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB);
      break;

    default:
      kaBaseException("Format of output matrix %s is missing", name);
      break;
  }
  ierr = VecView(v, fd); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

#if KADEBUG == 2
  fprintf(stderr, "SaveVec finished\n");
#endif  // if KADEBUG == 2
}  // SaveVecMatInfo

class PETScLSE {
 public:
  IndexType xLattice, yLattice, zLattice;  // dimensions of 3D problem
  PetscInt size;
  PetscMPIInt mpi_size, mpi_rank;
  PetscInt Istart, Iend;
  kaSharedMatrixN<double> m;  // Geometry matrix

  Mat A;       // system matrix from FD/FE discretization of Poisson equation
  IOType mot;  // io format of system matrix
  Vec X;       // vector of potentials
  IOType xot;  // io format of vectors of potentials
  Vec B;       // vector of current source densities
  IOType bot;  // io format of vectors of source densities

  KSP ksp;     // Krylov subspace method context
  PC pc;       // Preconditioner

  PETScLSE();
  virtual ~PETScLSE();

  void Set(const kaSharedMatrixN<double> &mArg) {
    m.New(4);
    m = mArg;
  }

  void Set(const IndexType xLatticeArg, const IndexType yLatticeArg, const IndexType zLatticeArg) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::Set(%d,%d,%d)\n", xLatticeArg, yLatticeArg, zLatticeArg);
#endif  // if KADEBUG
    xLattice = xLatticeArg;
    yLattice = yLatticeArg;
    zLattice = zLatticeArg;
    size     = xLattice*yLattice*zLattice;
  }

  void CreateA(PetscInt, PetscInt);
  void CreateA();
  void SetPreallocationA(const PetscInt[], const PetscInt[]);

  // add matrix elements with significant values
  inline void MatSetValues(PetscInt m, PetscInt n, PetscInt Indx[], PetscScalar v[], InsertMode addv) {
    /*PetscInt d=0;
       for(PetscInt i=0; i<n; i++)
       if (fabs(v[i])>1e-7) {
        Indx[d]=Indx[i];
        v[d]=v[i];
        d++;
       }
       PetscErrorCode ierr = ::MatSetValues(A, 1, &m, d, &Indx[0], v, addv);CHKERRQ(ierr); */
    PetscErrorCode ierr = ::MatSetValues(A, 1, &m, n, &Indx[0], v, addv); CHKERRQ(ierr);
  }

  void CreateX();
  void CreateB();
  void AddColumnVec(Vec *Stim, PetscInt row, double value);
  void MatZeroColumns(PetscInt numRows, const PetscInt *rows);
  void MatVecScale(Vec scale);
  void SetVoltage(PetscInt row, double val, bool set_vec = true);
  void SetCurrent(PetscInt row, double val);

  void LoadLatInfo(const char *file) {
    ::LoadLatInfo(file, this->m, this->xLattice, this->yLattice, this->zLattice);
    size = xLattice*yLattice*zLattice;
  }

  void SaveLatInfo(const char *file) {
    ::SaveLatInfo(file, this->m, this->xLattice, this->yLattice, this->zLattice);
  }

  void LoadA(const char *name);
  void SaveA(const char *name);

  void LoadX(const char *name) {LoadVec(name, X, this->xot);}

  void SaveX(const char *name) {SaveVec(name, X, this->xot);}

  void LoadB(const char *name) {LoadVec(name, B, this->bot);}

  void SaveB(const char *name) {SaveVec(name, B, this->bot);}

  void LoadVec(const char *name, Vec &v, IOType ot);
  void SaveVec(const char *name, const Vec &v, IOType ot);
  void SaveVecMatInfo(const char *name, const Vec &v, IOType ot, const kaSharedMatrixN<double> &m,
                      const IndexType xLattice, const IndexType yLattice, const IndexType zLattice);

  void InitSolve(KSPType ksptype  = KSPGMRES, PCType pctype = PCBJACOBI, PetscReal rtol = 1e-05,
                 PetscReal abstol = 1e-50, PetscReal dtol = 10000, PetscInt maxits = 10000) {
    PetscErrorCode ierr;

    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, ksptype); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, rtol, abstol, dtol, maxits); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, pctype); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  }

  //! Solve LSE
  inline void Solve(const Vec &B, Vec &X) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::Solve()\n");
#endif  // if KADEBUG

    PetscErrorCode ierr;
    ierr = KSPSolve(ksp, B, X); CHKERRQ(ierr);
  }
};  // class PETScLSE

#endif  // ifndef PETSCLSE_H
