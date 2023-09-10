/*! \file  PETScLSESolve.cpp
   \brief Solving stiffnes matrices for boundary conditions

   \version 1.0.0

   \date Created Gunnar Seemann (14.07.08)\n
   Last Modified Eike Wuelfers (28.04.11)

   \author Frank Sachse, CVRTI - University of Utah\n
   Institute of Biomedical Engineering\n
   Universitaet Karlsruhe (TH)\n
   http://www.ibt.uni-karlsruhe.de\n
   Copyright 2000-2008 - All rights reserved.

   \sa Synopsis \ref PETScLSESolve
 */

#include <PETScLSE.h>
#include <PETScLSEConditions.h>
#include <timer.h>

class PETScAIJSolve : public PETScLSEConditions, public PETScLSE {
  bool verbose;
  bool saveab;  /* save A and b after assembling */
  char savepre[1024];  /* filename prefix to save to */
  int mpirank;
  timer t;

 public:
  PETScAIJSolve(int argc, char **argv);
  ~PETScAIJSolve();
  void ApplyConditions();
};


PETScAIJSolve::PETScAIJSolve(int argc, char **argv) {
#if KADEBUG
  fprintf(stderr, "PETScAIJSolve<ValTyp>::PETScAIJSolve()\n");
#endif  // if KADEBUG

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  PetscErrorCode ierr;
  char *matname = NULL, *condname = NULL, *vecname = NULL;
  verbose = false;

  int is = 1;
  matname  = argv[is++];
  condname = argv[is++];
  vecname  = argv[is++];

  PetscBool flg;
  PetscOptionsHasName(PETSC_NULL, "-verbose", &flg);
  if (flg)
    this->verbose = true;
  PetscOptionsHasName(PETSC_NULL, "-bin", &flg);
  if (flg)
    this->mot = ot_bin;
  PetscOptionsHasName(PETSC_NULL, "-ascii", &flg);
  if (flg)
    this->mot = ot_ascii;
  PetscOptionsHasName(PETSC_NULL, "-matlab", &flg);
  if (flg)
    this->mot = ot_matlab;
  ierr = PetscOptionsGetString(PETSC_NULL, "-saveab", savepre, sizeof(savepre)/sizeof(char), &flg);
  if (flg) {
    this->saveab  = true;
    PETScLSE::bot = PETScLSE::mot;
  }

  PETScLSE::LoadA(matname);
  PETScLSE::CreateX();
  PETScLSE::CreateB();

  PETScLSEConditions::Load(condname);
  PETScLSEConditions::Parse(PETScLSE::size);

  ApplyConditions();

  if (saveab) {
    string fn(savepre);
    IOType tmpiot = PETScLSE::mot;
    PETScLSE::mot = ot_matlab;
    PETScLSE::SaveA((fn+".mtx").c_str());
    PETScLSE::mot = tmpiot;
    PETScLSE::SaveVec((fn+".rhs").c_str(), PETScLSE::B, ot_matlab);
  }
  t.start();
  PETScLSE::InitSolve();
  PETScLSE::Solve(PETScLSE::B, PETScLSE::X);
  t.stop();

  if (!mpirank)
    cout << t << endl;

  PETScLSE::SaveX(vecname);
}

PETScAIJSolve::~PETScAIJSolve() {
#if KADEBUG
  fprintf(stderr, "PETScAIJSolve<ValTyp>::~PETScAIJSolve()\n");
#endif  // if KADEBUG
}

void PETScAIJSolve::ApplyConditions() {
#if KADEBUG
  fprintf(stderr, "PETScAIJSolve<ValTyp>::ApplyConditions\n");
#endif  // if KADEBUG
  PetscErrorCode ierr;
  PetscInt xLattice = PETScLSE::xLattice;
  PetscInt xyLattice = xLattice*PETScLSE::yLattice;
  PetscInt xyzLattice = xyLattice*PETScLSE::zLattice;
  PetscInt Istart, Iend;

  ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
  if (verbose)
    printf("[%d] Processing range: %d %d\n", mpirank, Istart, Iend);

  vector<BoundaryStorage> BS;
  Vec Stim;
  VecDuplicate(B, &Stim);

  for (auto it = PETScLSEConditions::SC.begin();
       it != PETScLSEConditions::SC.end(); ++it) {
    PetscInt p = it->position;
    if ((p >= Istart) && (p < Iend)) {
      switch (it->cflag) {
        case CT_I:
          SetCurrent(p, it->val);
          break;
        case CT_U:
          BoundaryStorage lBS;
          lBS.index = p;
          lBS.value = it->val;
          BS.push_back(lBS);
          AddColumnVec(&Stim, p, lBS.value);
          break;
        default:
          kaBaseException("Unknown boundary condition type");
          break;
      }
    }
  }

  ierr = VecAssemblyBegin(Stim); CHKERRQ(ierr);  // Vec Stim has to be assembled because of preceded call VecSetValues()
                                                 // (in function SetCurrent)
  ierr = VecAssemblyEnd(Stim); CHKERRQ(ierr);

  ierr = VecAXPY(B, -1.0, Stim); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(B); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(B); CHKERRQ(ierr);

  vector<PetscInt> indizes(BS.size());
  vector<PetscScalar> values(BS.size());

  vector<BoundaryStorage>::iterator iBS, bBS = BS.begin(), eBS = BS.end();
  int i = 0;
  for (iBS = bBS; iBS != eBS; iBS++, i++) {
    indizes[i] = (*iBS).index;
    values[i]  = (*iBS).value;
  }

  ierr = VecSetValues(B, BS.size(), indizes.data(), values.data(), INSERT_VALUES); CHKERRQ(ierr);

  MatZeroColumns(BS.size(), indizes.data());

  ierr = MatZeroRows(A, BS.size(), indizes.data(), 1.0, 0, 0); CHKERRQ(ierr);

  ierr = VecDestroy(&Stim); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(B); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(B); CHKERRQ(ierr);
}  // PETScAIJSolve::ApplyConditions

int main(int argc, char **argv) {
#if KADEBUG
  fprintf(stderr, "main\n");
#endif  // if KADEBUG
  PetscErrorCode ierr;
  PetscInitialize(&argc, &argv, (char *)0, "");

  if (argc < 3) {
    cerr << argv[0] << " <stiffness matrix filename> <condition filename> <solution vector filename>" << endl;
    cerr << "\t[-verbose]" << endl;
    exit(-1);
  }
  try {
    PETScAIJSolve LFH2P(argc, argv);
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    ierr = PetscFinalize(); CHKERRQ(ierr);
    exit(-1);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

/*! \page PETScAIJSolve PETScLSESolve
   Solves LSE for boundary conditions.

   \section SYNOPSIS_PETScLSESolve SYNOPSIS
   PETScLSESolve \<stiffness matrix filename\> \<condition filename\> \<solution vector filename\>\n

   \section OPTIONS_PETScLSESolve OPTIONS
   \param "<stiffness matrix filename>" The stiffness matrix, generated e.g. by LatticeFECC2PETScAIJ
   \param "<condition filename>" A file containing the stimulation Points (x, y, z, type, amplitude) with type eithe "I"
      or "U"
   \param "<solution vector filename>" Solution vector in PETSc-format

   \section DESCRIPTION_PETScLSESolve DESCRIPTION
   To be done...

   \section SOURCE_PETScAIJSolve SOURCE
   PETScLSESolve.cpp PETScLSE.cpp PETScLSE.h

   \section SEEALSO_PETScLSESolve SEE ALSO
   \ref PETScLSE

   \section CHANGELOG_PETScLSESolve CHANGELOG
   V1.0.0 - 08.04.2009 (Meike Karl): Starting with changelog\n
   V1.0.1 - 29.04.2009 (Meike Karl): Bugfix: Number of processes had significant influence onto the solution of Ve
      (using extracellular stimuli I or U); now fixed (paralellization problem)
 */
