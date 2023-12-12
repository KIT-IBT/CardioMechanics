/*
 * File: PETScLSE.cpp
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


#include <PETScLSE.h>


PETScLSE::PETScLSE() {
#if KADEBUG == 2
    fprintf(stderr, "PETScLSE::PETScLSE()\n");
#endif  // if KADEBUG == 2
    A   = NULL;
    X   = NULL;
    B   = NULL;
    mot = ot_bin;
    xot = ot_bin;
    bot = ot_bin;
    
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
}

PETScLSE::~PETScLSE() {
#if KADEBUG == 2
    fprintf(stderr, "PETScLSE::~PETScLSE()\n");
#endif  // if KADEBUG == 2
    
    PetscErrorCode ierr;
    
    if (ksp) {
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    }
    if (A) {
        ierr = MatDestroy(&A); CHKERRQ(ierr);
    }
    if (X) {
        ierr = VecDestroy(&X); CHKERRQ(ierr);
        X    = NULL;
    }
    if (B) {
        ierr = VecDestroy(&B); CHKERRQ(ierr);
    }
}

void PETScLSE::CreateA(PetscInt DiagonalElements, PetscInt NonDiagonalElements) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::CreateA(PetscInt, PetscInt)\n");
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, size, size, DiagonalElements, PETSC_NULL,
                        NonDiagonalElements, PETSC_NULL, &A); CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
}

void PETScLSE::CreateA() {
#if KADEBUG
    fprintf(stderr, "PETScLSE::CreateA()\n");
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
    
    // ierr = MatSetType(A, MATMPIAIJ); CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size); CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); CHKERRQ(ierr);
    ierr = MatSetUp(A); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
}

void PETScLSE::SetPreallocationA(const PetscInt DiagonalElements[], const PetscInt NonDiagonalElements[]) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::SetPreallocationA()\n");
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    ierr = MatMPIAIJSetPreallocation(A, 27, DiagonalElements, 12, NonDiagonalElements); CHKERRQ(ierr);
}

void PETScLSE::CreateX() {
#if KADEBUG
    fprintf(stderr, "PETScLSE::CreateX()\n");
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &X); CHKERRQ(ierr);
    ierr = VecSetSizes(X, PETSC_DECIDE, size); CHKERRQ(ierr);
    ierr = VecSetType(X, VECMPI); CHKERRQ(ierr);
    ierr = VecSet(X, 0.); CHKERRQ(ierr);
}

void PETScLSE::CreateB() {
#if KADEBUG
    fprintf(stderr, "PETScLSE::CreateB()\n");
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &B); CHKERRQ(ierr);
    ierr = VecSetSizes(B, PETSC_DECIDE, size); CHKERRQ(ierr);
    ierr = VecSetType(B, VECMPI); CHKERRQ(ierr);
    ierr = VecSet(B, 0.); CHKERRQ(ierr);
}

void PETScLSE::AddColumnVec(Vec *Stim, PetscInt row, PetscScalar value) {
    PetscErrorCode ierr;
    PetscInt Istart, Iend;
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
    
    if ((row >= Istart) && (row < Iend)) {
        ierr = MatGetRow(A, row, &ncols, &cols, &vals); CHKERRQ(ierr);
        vector<PetscInt> indizes(ncols);
        vector<PetscScalar> values(ncols);
        for (int j = 0; j < ncols; j++) {
            indizes[j] = cols[j];
            values[j]  = vals[j]*value;
        }
        ierr = MatRestoreRow(A, row, &ncols, &cols, &vals); CHKERRQ(ierr);
        
        ierr = VecSetValues(*Stim, ncols, indizes.data(), values.data(), ADD_VALUES); CHKERRQ(ierr);
    }
}

void PETScLSE::MatVecScale(Vec scale) {
    PetscErrorCode ierr;
    PetscInt Istart, Iend;
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
    PetscScalar *pscale;
    VecGetArray(scale, &pscale);
    PetscInt i = Istart;
    for (; i < Iend; i++) {
        ierr = MatGetRow(A, i, &ncols, &cols, &vals); CHKERRQ(ierr);
        vector<PetscInt> indizes(ncols);
        vector<PetscScalar> values(ncols);
        for (PetscInt j = 0; j < ncols; j++) {
            indizes[j] = cols[j];
            values[j]  = vals[j]*(*(pscale+i));
            if (cols[j] == i) {
                if (values[j] < 1e-9)
                    values[j] = 1e-9;
            }
        }
        ierr = MatRestoreRow(A, i, &ncols, &cols, &vals); CHKERRQ(ierr);
        
        ierr = ::MatSetValues(A, 1, &i, ncols, indizes.data(), values.data(), INSERT_VALUES);  CHKERRQ(ierr);
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }
    VecRestoreArray(scale, &pscale);
}  // PETScLSE::MatVecScale

void PETScLSE::MatZeroColumns(PetscInt numRows, const PetscInt *rows) {
    PetscErrorCode ierr;
    PetscInt Istart, Iend;
    
    vector<PetscInt> ncols(numRows);
    const PetscInt  *cols;
    
    vector<PetscInt *> indizes(numRows);
    vector<PetscScalar *> values(numRows);
    
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
    int i = 0;
    for (i = 0; i < numRows; i++)
        if ((rows[i] >= Istart) && (rows[i] < Iend)) {
            ierr       = MatGetRow(A, rows[i], &ncols[i], &cols, PETSC_NULL); CHKERRQ(ierr);
            indizes[i] = new PetscInt[ncols[i]];
            values[i]  = new PetscScalar[ncols[i]];
            for (int j = 0; j < ncols[i]; j++) {
                indizes[i][j] = cols[j];
                values[i][j]  = 0.0;
            }
            ierr = MatRestoreRow(A, rows[i], &ncols[i], &cols, PETSC_NULL); CHKERRQ(ierr);
        }
    
    for (i = 0; i < numRows; i++)
        if ((rows[i] >= Istart) && (rows[i] < Iend)) {
            ierr = ::MatSetValues(A, ncols[i], indizes[i], 1, &rows[i], values[i], INSERT_VALUES);  CHKERRQ(ierr);
            delete indizes[i];
            delete values[i];
        }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
}  // PETScLSE::MatZeroColumns

void PETScLSE::SetVoltage(PetscInt row, double V, bool set_vec) {
#if KADEBUG == 2
    fprintf(stderr, "PETScLSE::SetVoltage(%d, %lf)\n", row, V);
#endif  // if KADEBUG == 2
    
    PetscErrorCode ierr;
    
    ierr = MatZeroRows(A, 1, &row, 1.0, 0, 0); CHKERRQ(ierr);
    if (set_vec)
        ierr = VecSetValues(B, 1, &row, &V, INSERT_VALUES); CHKERRQ(ierr);
}

void PETScLSE::SetCurrent(PetscInt row, double I) {
#if KADEBUG == 2
    fprintf(stderr, "PETScLSE::SetCurrent(%d,%lf)\n", row, I);
#endif  // if KADEBUG == 2
    
    PetscErrorCode ierr;
    
    ierr = VecSetValues(B, 1, &row, &I, ADD_VALUES); CHKERRQ(ierr);
}

//! Load matrix of specified format
void PETScLSE::LoadA(const char *name) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::LoadA('%s')\n", name);
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    PetscViewer fd;  /* viewer */
    
    switch (this->mot) {
        case ot_bin:
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_READ, &fd); CHKERRQ(ierr);
            break;
            
        case ot_ascii:
        case ot_matlab:
            throw kaBaseException("Format of input matrix %s is not supported", name);
            break;
            
        default:
            throw kaBaseException("Format of input matrix %s is missing", name);
            break;
    }
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetType(A, MATMPIAIJ);
    ierr = MatLoad(A, fd); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
    
    try {
        LoadLatInfo(name);
    } catch (...) {
        PetscInt m;
        ierr     = MatGetSize(A, &m, NULL); CHKERRQ(ierr);
        zLattice = size = m;
        xLattice = yLattice = 1;
        this->m.New(4);
        this->m.Identity();
    }
}  // PETScLSE::LoadA

//! Save matrix in specified format
void PETScLSE::SaveA(const char *name) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::SaveA('%s')\n", name);
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    PetscViewer fd;  /* viewer */
    int rank = 0;
    
    switch (this->mot) {
        case ot_bin:
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
                SaveLatInfo(name);
            break;
            
        case ot_ascii:
            ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
            break;
            
        case ot_matlab:
            ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
            ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB);
            break;
            
        default:
            throw kaBaseException("Format of output matrix %s is missing", name);
            break;
    }
    ierr = MatView(A, fd); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
}  // PETScLSE::SaveA

//! Load vector of specified format
void PETScLSE::LoadVec(const char *name, Vec &v, IOType ot) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::LoadVec('%s')\n", name);
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    PetscViewer fd;  /* viewer */
    
    switch (ot) {
        case ot_bin:
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_READ, &fd); CHKERRQ(ierr);
            try {
                LoadLatInfo(name);
            } catch (...) {
                PetscInt m;
                ierr     = MatGetSize(A, &m, NULL); CHKERRQ(ierr);
                zLattice = size = m;
                xLattice = yLattice = 1;
                this->m.New(4);
                this->m.Identity();
            }
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
    
#if KADEBUG
    fprintf(stderr, "PETScLSE::LoadVec finished\n");
#endif  // if KADEBUG
}  // PETScLSE::LoadVec

//! Save vector in specified format
void PETScLSE::SaveVec(const char *name, const Vec &v, IOType ot) {
#if KADEBUG
    fprintf(stderr, "PETScLSE::SaveVec('%s')\n", name);
#endif  // if KADEBUG
    
    PetscErrorCode ierr;
    PetscViewer fd;  /* viewer */
    
    switch (ot) {
        case ot_bin:
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
            SaveLatInfo(name);
            break;
            
        case ot_ascii:
            ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
            break;
            
        case ot_matlab:
            ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &fd); CHKERRQ(ierr);
            ierr = PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB);
            break;
            
        default:
            throw kaBaseException("Format of output matrix %s is missing", name);
            break;
    }
    ierr = VecView(v, fd); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
    
#if KADEBUG == 2
    fprintf(stderr, "PETScLSE::SaveVec finished\n");
#endif  // if KADEBUG == 2
}  // PETScLSE::SaveVec
