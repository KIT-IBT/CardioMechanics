/*
 * File: BidomainMatrixGenerator.cpp
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


#include "FEQuadrature.h"

#include <petsc.h>
#include <petscerror.h>

#include <Material.h>
#include <kaMask.h>

#include "VTKIOHelper.h"
#include "ProgressBar.h"
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkCell.h>

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <set>
#include <array>

#include <Matrix3.h>
#include <Matrix4.h>
#include <Vector4.h>


using namespace math_pack;

// TODO
// - connectivity (for reordering)
// - fix theta for coupling
// - sigma handling is bad, lots of duplication
// - Parallel assembly

#define S(a) std::string(a)

#define ERR(a, b) do {                           \
    std::cerr<<"ERROR: "<<a<<std::endl; exit(b); \
} while (0)
#define WARN(a) do {                      \
    std::cerr<<"WARNING: "<<a<<std::endl; \
} while (0)
#define CHKE(n) CHKERRABORT(PETSC_COMM_WORLD, n)

typedef uint8_t mat_t;
#define MAX_MAT 256

void saveMatrix(std::string const &, Mat);
void printUsage(const char *);
void printHelp(void);

enum {
  ERR_PARAMETER = 1,
  ERR_SANITY    = 2,
  ERR_MATERIAL  = 3,
  ERR_SIGMA     = 4,
  ERR_CELL      = 5,
  ERR_GAUSS     = 6,
};

struct MaterialWrapper {
  MaterialWrapper(std::string const &fn_) : fn(fn_) {
    try {
      ml.reset(new MaterialListe(fn.c_str()));
    } catch (kaBaseException &e) {
      ERR("Material list in file " << fn << ":\n" << e, ERR_MATERIAL);
    }
  }

  Material *Get(int i) {
    Material *m = GetIfExists(i);

    if (!m) {
      ERR("Unable to find material class " << i << " in " << fn, ERR_MATERIAL);
    }
    return m;
  }

  Material *GetIfExists(int i) {
    Material *m = nullptr;

    try {
      m = ml->Suchen(i);
    } catch (kaBaseException &e) {
      ERR("Material list in file " << fn << ":\n" << e, ERR_MATERIAL);
    }
    return m;
  }

  std::string fn;
  std::unique_ptr<MaterialListe> ml;
};

static bool gzip = false;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printUsage(argv[0]);
    return 1;
  } else if (argv[1] == S("-help")) {
    printUsage(argv[0]);
    printHelp();

    return 1;
  }

  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKE(ierr);

  std::string input  = argv[1];
  std::string prefix = input.substr(0, input.rfind('.') + 1);
  kaMask bathMask;

  std::array<int, MAX_MAT> tissuePrecedence({{MAX_MAT}});
  {
    int i = 0;
    for (int m = 0; m < MAX_MAT; m++) {
      tissuePrecedence[m] = i++;
    }
  }

  bool bath = false, mono = false, bi = false, combined = false, coupled = false;
  bool lumping = false, fullLumping = false;
  bool verbose = false, isotropic = false, ignoreMissing = false;
  double theta_mass = 1;
  std::string matfBath, matfIntra, matfExtra;
  double frequency = 0, betaCm_dt = .0, theta = .0;
  double betaDeltaS                       = .0;
  int bdsp                                = 2;
  int gauss                               = 0;
  std::unique_ptr<FETetQuadrature> gaussq = nullptr;

  /* Parameter parsing */
  for (int arg = 2; arg < argc; ++arg) {
    if (S("-o") == argv[arg]) {
      prefix  = argv[++arg];
      prefix += ".";
    } else if (S("-v") == argv[arg]) {
      verbose = true;
    } else if ((S("-bath") == argv[arg]) && (arg+1 < argc)) {
      bathMask.Set(argv[++arg]);
      if ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
        matfBath = argv[++arg];
      }
      bath = true;
    } else if ((S("-mono") == argv[arg]) && (arg+1 < argc)) {
      matfIntra = argv[++arg];
      if ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
        matfExtra = argv[++arg];
      }
      mono = true;
    } else if ((S("-bi") == argv[arg]) && (arg+2 < argc)) {
      matfIntra = argv[++arg];
      matfExtra = argv[++arg];
      bi        = true;
    } else if ((S("-tissue") == argv[arg]) && (arg+1 < argc)) {
      tissuePrecedence.fill(MAX_MAT);
      tissuePrecedence[0] = MAX_MAT + 1;
      std::stringstream ss(argv[++arg]);
      std::string item;
      int i = 0;
      while (std::getline(ss, item, ',')) {
        tissuePrecedence[(mat_t)std::stoi(item)] = i++;
      }
    } else if (S("-combined") == argv[arg]) {
      combined = true;
    } else if (S("-fullLumping") == argv[arg]) {
      fullLumping = true;
    } else if ((S("-massLumping") == argv[arg]) && (arg+1 < argc)) {
      theta_mass = atof(argv[++arg]);
      lumping    = true;
    } else if ((S("-coupled") == argv[arg]) && (arg+1 < argc)) {
      betaCm_dt = atof(argv[++arg]);
      coupled   = true;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
        theta = atof(argv[++arg]);
      }
    } else if ((S("-bds") == argv[arg]) && (arg+1 < argc)) {
      betaDeltaS = atof(argv[++arg]);
      if ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
        bdsp = atoi(argv[++arg]);
      }
    } else if ((S("-freq") == argv[arg]) && (arg+1 < argc)) {
      frequency = atof(argv[++arg]);
    } else if (S("-iso") == argv[arg]) {
      isotropic = true;
    } else if (S("-ignoreMissing") == argv[arg]) {
      ignoreMissing = true;
    } else if (S("-gzip") == argv[arg]) {
      gzip = true;
    } else if ((S("-gauss") == argv[arg]) && (arg+1 < argc)) {
      gauss = atoi(argv[++arg]);
    } else {
      ERR("Argument " << argv[arg] << " not understood or not enough parameters.", ERR_PARAMETER);
    }
  }

  /* Sanity checks */
  if (mono && bi) {
    ERR("Can only do either mono- or bidomain.", ERR_SANITY);
  } else if (coupled && combined) {
    ERR("Can only do either coupled or combined bidomain matrices.", ERR_SANITY);
  } else if (bath && !bi) {
    ERR("Bath requires bidomain.", ERR_SANITY);
  } else if (combined && !bi) {
    ERR("-coupled requires bidomain.", ERR_SANITY);
  } else if (combined && !bi) {
    ERR("-combined requires bidomain.", ERR_SANITY);
  } else if (coupled && fullLumping) {
    ERR("-coupled can not be used with -fullLumping. Use -massLumping 1 instead.", ERR_SANITY);
  } else if (!mono && !bi) {
    ERR("Either mono- or bidomain must be selected.", ERR_SANITY);
  } else if ((theta < .0) || (theta > 1.0)) {
    ERR("Theta must be between 0 and 1.", ERR_SANITY);
  }
  if ((theta_mass < .0) || (theta_mass > 1.0)) {
    WARN("Theta for mass lumping is not between 0 and 1.");
  }
  if (ignoreMissing && !bi) {
    WARN("-ignoreMissing should only be used with -bi.");
  }

  if (gauss) {
    try {
      gaussq = FETetQuadrature::GetDefaultRuleByPoints(gauss);
    } catch (std::exception &e) {
      ERR("Gauss quadrature rule for " << gauss << " points:\n" << e.what(), ERR_GAUSS);
    }
  }

  if (coupled) {
    /* When buildiung coupled matrix, conductivities in the overlap of both sets
     * have to have combined conductivities. */
    combined = true;
  }

  /* Load data */
  vtkSmartPointer<vtkDataSet> mesh = VTKIOHelper::Read(argv[1]);

  const vtkIdType numCells = mesh->GetNumberOfCells();
  PetscInt numIntraCells   = 0;
  const vtkIdType numNodes = mesh->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> matArray = mesh->GetCellData()->GetArray("Material");

  std::set<vtkIdType> bathNodes;
  std::set<vtkIdType> intraNodes;

  std::map<mat_t, std::array<double, 4>> intraSigmas;
  std::map<mat_t, std::array<double, 4>> extraSigmas;

  kaMask hasBath, hasIntra, hasExtra;

  std::vector<PetscInt> cellsPerNode(numNodes, 0);
  std::vector<PetscInt> nodeNeighbours(numNodes, 1);

  for (vtkIdType i = 0; i < numCells; ++i) {
    mat_t m                    = (mat_t)matArray->GetTuple1(i);
    vtkSmartPointer<vtkCell> c = mesh->GetCell(i);
    if (!(bath && bathMask.m[m])) {
      numIntraCells++;
    }

    for (int j = 0; j < c->GetNumberOfPoints(); ++j) {
      vtkIdType n = c->GetPointId(j);
      if (bath && bathMask.m[m]) {
        bathNodes.insert(n);
        hasBath.m[m] = true;
      } else {
        intraNodes.insert(n);
        hasIntra.m[m] = true;
      }
      cellsPerNode[n]++;
      nodeNeighbours[n] += c->GetNumberOfPoints() - 1;
    }
  }

  const vtkIdType numIntraNodes = intraNodes.size();
  IS intra_is                   = NULL;
  std::map<PetscInt, PetscInt> globalToIntra;
  std::vector<PetscInt> intraNeighbours;

  if (numIntraNodes != numNodes) {
    intraNeighbours.resize(numIntraNodes);
    PetscInt node = 0;
    PetscInt *intraIdcs;  // must use petscmalloc so that iscreate can takeownership
    ierr = PetscMalloc(numIntraNodes * sizeof(PetscInt), &intraIdcs); CHKE(ierr);
    for (const auto &in : intraNodes) {
      intraIdcs[node]       = in;
      globalToIntra[in]     = node;
      intraNeighbours[node] = nodeNeighbours[in];
      ++node;
    }
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, numIntraNodes, intraIdcs,
                           PETSC_OWN_POINTER, &intra_is); CHKE(ierr);
  } else {
    intraNeighbours = nodeNeighbours;
  }

  /* Read material conductivities */
  if (bath && (bathNodes.size() == 0)) {
    WARN("Bath mask does not match any cells. Bath domain is empty.");
  } else if (bath && matfBath.size()) {
    MaterialWrapper matl(matfBath);

    for (int i = 0; i < MAX_MAT; ++i) {
      if (hasBath.m[i]) {
        Material *material = matl.Get(i);
        double s           = material->LFkappa::Hole(frequency);
        if (s < 0) {
          ERR("Conductivity for class " << i << " is negative: " << s, ERR_SIGMA);
        }

        extraSigmas[i] = {{s,
          s *material->HoleAnisotropyX(),
          s *material->HoleAnisotropyY(),
          s *material->HoleAnisotropyZ()}
        };
      }
    }
  } else {
    hasExtra = hasBath;
  }

  if (intraNodes.size() == 0) {
    WARN("Bath mask matches all cells. Intracellular domain is empty.");
  }

  if (verbose) {
    std::cerr << "Domains:\n";
    std::cerr << "  Intra: " << numIntraNodes << std::endl;
    std::cerr << "  Extra: " << numNodes << std::endl;
    std::cerr << "  Bath:  " << numNodes - numIntraNodes << std::endl;
  }

  if (bi) {
    MaterialWrapper mati(matfIntra);
    MaterialWrapper mate(matfExtra);

    for (int i = 0; i < MAX_MAT; ++i) {
      if (hasIntra.m[i] || hasExtra.m[i]) {  /* all intra elements have extra as well */
        Material *material = mate.Get(i);
        double s           = material->LFkappa::Hole(frequency);
        if (s < 0) {
          ERR("Conductivity for class " << i << " is negative: " << s, ERR_SIGMA);
        }

        extraSigmas[i] = {{s,
          s *material->HoleAnisotropyX(),
          s *material->HoleAnisotropyY(),
          s *material->HoleAnisotropyZ()}
        };

        if (hasIntra.m[i]) {
          material = mati.GetIfExists(i);
          if (!material) {
            if (ignoreMissing) {
              s              = .0;
              intraSigmas[i] = {{.0, .0, .0, .0}
              };
              if (verbose) {
                std::cerr << "Ignoring missing intracellular conductivity for class " << i << "\n";
              }
            } else {
              ERR("Unable to find material class " << i << " in " << matfIntra, ERR_MATERIAL);
            }
          } else {
            s = material->LFkappa::Hole(frequency);
            if (s < 0) {
              ERR("Conductivity for class " << i << " is negative: " << s, ERR_SIGMA);
            }
            intraSigmas[i] = {{s,
              s *material->HoleAnisotropyX(),
              s *material->HoleAnisotropyY(),
              s *material->HoleAnisotropyZ()}
            };
          }
        }
      }
    }
  } else {  // mono
    MaterialWrapper mati(matfIntra);
    for (int i = 0; i < MAX_MAT; ++i) {
      if (hasIntra.m[i]) {
        Material *material = mati.Get(i);
        double s           = material->LFkappa::Hole(frequency);
        if (s < 0) {
          ERR("Conductivity for class " << i << " is negative: " << s, ERR_SIGMA);
        }

        intraSigmas[i] = {{s,
          s *material->HoleAnisotropyX(),
          s *material->HoleAnisotropyY(),
          s *material->HoleAnisotropyZ()}
        };
      }
    }

    /* compute monodomain conductivity from intra and extra (geom. mean) */
    /* assumes parallel circuit */
    if (matfExtra.size()) {
      MaterialListe mate(matfExtra.c_str());
      for (int i = 0; i < MAX_MAT; ++i) {
        if (hasIntra.m[i]) {
          Material *material = mate.Suchen(i);
          if (!material) {
            ERR("Unable to find material class " << i << " in " << matfExtra, ERR_MATERIAL);
          }
          double s = material->LFkappa::Hole(frequency);
          if (s < 0) {
            ERR("Conductivity for class " << i << " is negative: " << s, ERR_SIGMA);
          }
          double se[4] = {s,                              s *material->HoleAnisotropyX(),
                          s *material->HoleAnisotropyY(), s *material->HoleAnisotropyZ()};

          intraSigmas[i][0] *= se[0] / (intraSigmas[i][0] + se[0]);
          intraSigmas[i][1] *= se[1] / (intraSigmas[i][1] + se[1]);
          intraSigmas[i][2] *= se[2] / (intraSigmas[i][2] + se[2]);
          intraSigmas[i][3] *= se[3] / (intraSigmas[i][3] + se[3]);
        }
      }
    }
  }

  /* Preallocate matrices */
  Mat intraMatrix = NULL, extraMatrix = NULL, massMatrix = NULL;
  Mat gaussMatrix = NULL, interpolMatrix = NULL;

  if (bi) {
    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numNodes, numNodes, 0,
                           nodeNeighbours.data(), &extraMatrix); CHKE(ierr);
  }


  ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numIntraNodes, numIntraNodes, 0,
                         intraNeighbours.data(), &intraMatrix); CHKE(ierr);

  ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numIntraNodes, numIntraNodes, 0,
                         intraNeighbours.data(), &massMatrix); CHKE(ierr);

  if (gauss) {
    PetscInt n = numIntraCells * gauss;
    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, n, numIntraNodes, 4,  /* TODO tet specific... */
                           NULL, &interpolMatrix); CHKE(ierr);

    std::vector<PetscInt> gaussPtsPerIntraNode(numIntraNodes);
    if (intra_is) {
      PetscInt const *intraIdcs = NULL;
      ierr = ISGetIndices(intra_is, &intraIdcs); CHKE(ierr);
      for (PetscInt i = 0; i < numIntraNodes; ++i) {
        gaussPtsPerIntraNode[i] = cellsPerNode[intraIdcs[i]] * gauss;
      }
      ierr = ISRestoreIndices(intra_is, &intraIdcs); CHKE(ierr);
    } else {
      for (PetscInt i = 0; i < numIntraNodes; ++i) {
        gaussPtsPerIntraNode[i] = cellsPerNode[i] * gauss;
      }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numIntraNodes, n, 0,
                           gaussPtsPerIntraNode.data(), &gaussMatrix); CHKE(ierr);
  }

  std::vector<mat_t> nodeMaterial(numNodes, 0);

  vtkSmartPointer<vtkDataArray> phiArray   = NULL;
  vtkSmartPointer<vtkDataArray> thetaArray = NULL;

  if (!isotropic) {
    phiArray   = mesh->GetCellData()->GetArray("Phi");
    thetaArray = mesh->GetCellData()->GetArray("Theta");

    if (!phiArray || !thetaArray) {
      vtkSmartPointer<vtkDataArray> fiberArray = mesh->GetCellData()->GetArray("Fiber");
      if (!fiberArray) {
        fiberArray = mesh->GetCellData()->GetArray("Orientation");
      }

      if (fiberArray && (fiberArray->GetNumberOfComponents() == 3)) {
        phiArray   = vtkSmartPointer<vtkDoubleArray>::New();
        thetaArray = vtkSmartPointer<vtkDoubleArray>::New();
        for (vtkIdType i = 0; i < fiberArray->GetNumberOfTuples(); ++i) {
          double *ori = fiberArray->GetTuple3(i);
          if ((ori[0] == .0) && (ori[1] == .0) && (ori[2] == .0)) {
            phiArray->InsertNextTuple1(255);
            thetaArray->InsertNextTuple1(255);
          } else {
            double phi   = ori[0] ? atan(ori[1] / ori[0]) : .0;
            double theta = acos(ori[2] / sqrt(ori[0] * ori[0] + ori[1] * ori[1] + ori[2] * ori[2]));
            if (phi < 0) {
              phi = M_PI + phi;
            }
            phiArray->InsertNextTuple1(phi*254.0/M_PI);
            thetaArray->InsertNextTuple1(theta*254.0/M_PI);
          }
        }
      }
    }

    if (phiArray && thetaArray) {
      for (int c = 0; c < numCells; ++c) {
        double phi   = phiArray->GetTuple1(c);
        double theta = thetaArray->GetTuple1(c);
        if ((phi < 0) || (phi > 255) || (theta < 0) || (theta > 255)) {
          std::cerr << "Cell " << c << " has phi=" << phi << " and theta=" << theta << std::endl;
        }
      }
    }
  }

  if (verbose) {
    if (!phiArray || !thetaArray) {
      std::cerr << "Fiber orientation ignored or not defined." << std::endl;
    } else {
      std::cerr << "Fiber orientation defined." << std::endl;
    }
  }

  PetscInt intraCell = 0;
  ProgressBar pb("Assembling", numCells);
  for (vtkIdType c = 0; c < numCells; ++c) {
    pb.set(c);

    vtkSmartPointer<vtkCell> cl = mesh->GetCell(c);
    PetscInt np                 = cl->GetNumberOfPoints();

    mat_t m = (mat_t)matArray->GetTuple1(c);

    Matrix3<double> sigma_i = Matrix3<double>::Identity() * intraSigmas[m][0],
                    sigma_e = Matrix3<double>::Identity() * extraSigmas[m][0];

    /* TODO signs in 3rd column are opposite in FDYZ.h. Does that make a
     * difference? */

    /* TODO bds fÃ¼r isotrop... */

    if (phiArray && thetaArray) {
      double phi   = phiArray->GetTuple1(c);
      double theta = thetaArray->GetTuple1(c);

      if ((theta != 255) && ((phi != 255) || (theta == 0) || (theta == 127))) {
        phi   = M_PI * phi / 254.0;
        theta = M_PI * theta / 254.0;

        Vector3<double> sigma_if = {intraSigmas[m][1], intraSigmas[m][2], intraSigmas[m][3]};
        Vector3<double> sigma_ef = {extraSigmas[m][1], extraSigmas[m][2], extraSigmas[m][3]};

        /* Rotation matrices from IBT wiki [[Lattice]] */
        Matrix3<double> r = Matrix3<double>(
          cos(phi), -sin(phi), 0,
          sin(phi),  cos(phi), 0,
          0,         0, 1) *
          Matrix3<double>(
          sin(theta), 0, -cos(theta),
          0,          1,           0,
          cos(theta), 0,  sin(theta));

        if (betaDeltaS != .0) {
          Vector3<double> f = r * Vector3<double>{1, 0, 0};
          Vector3<double> s = r * Vector3<double>{0, 1, 0};
          Vector3<double> n = r * Vector3<double>{0, 0, 1};

          double h_f = .0, h_s = .0, h_n = .0, min_f = HUGE, min_s = HUGE, min_n = HUGE;
          for (vtkIdType v = 0; v < cl->GetNumberOfPoints(); ++v) {
            Vector3<double> p(cl->GetPoints()->GetPoint(v));
            p    /= 1000.0;
            h_f   = ::max(h_f, f * p);
            h_s   = ::max(h_s, s * p);
            h_n   = ::max(h_n, n * p);
            min_f = ::min(h_f, f * p);
            min_s = ::min(h_s, s * p);
            min_n = ::min(h_n, n * p);
          }

          h_f = pow(h_f - min_f, bdsp);
          h_s = pow(h_s - min_s, bdsp);
          h_n = pow(h_n - min_n, bdsp);

          // TODO pow h, p
          sigma_if(1) += betaDeltaS * sqrt(sigma_if(0) / sigma_if(1)) * h_s;
          sigma_if(2) += betaDeltaS * sqrt(sigma_if(0) / sigma_if(2)) * h_n;
          sigma_if(0) += betaDeltaS * h_f;

          if (bi) {
            sigma_ef(1) += betaDeltaS * sqrt(sigma_ef(0) / sigma_ef(1)) * h_s;
            sigma_ef(2) += betaDeltaS * sqrt(sigma_ef(0) / sigma_ef(2)) * h_n;
            sigma_ef(0) += betaDeltaS * h_f;
          }
        }

        if (hasIntra.m[m]) {
          sigma_i = r * Matrix3<double>(
            sigma_if(0), .0, .0,
            .0, sigma_if(1), .0,
            .0, .0, sigma_if(2)) * r.GetTranspose();
        }

        if (bi) {  /* bath or not doesn't matter, all sigmas are in extraSigmas */
          sigma_e = r * Matrix3<double>(
            sigma_ef(0), .0, .0,
            .0, sigma_ef(1), .0,
            .0, .0, sigma_ef(2)) * r.GetTranspose();
        }
      }
    } else if (betaDeltaS != .0) {
      throw "Not implemented :(\n";
    }

    if (bi && combined && hasIntra.m[m]) {
      sigma_e += sigma_i;
    }

    for (int i = 0; i < np; ++i) {
      vtkIdType n = cl->GetPointId(i);
      if (hasIntra.m[m] && !hasIntra.m[nodeMaterial[n]]) {
        nodeMaterial[n] = m;
      } else if ((hasIntra.m[m] == hasIntra.m[nodeMaterial[n]]) &&
                 (tissuePrecedence[m] < tissuePrecedence[nodeMaterial[n]]) ) {
        nodeMaterial[n] = m;
      }
    }

    std::vector<double> ki(np *np);
    std::vector<double> ke(np *np);
    std::vector<double> mm(np *np);

    std::vector<PetscInt> p(np);
    for (vtkIdType i = 0; i < np; ++i) {
      p[i] = cl->GetPointId(i);
    }

    switch (cl->GetCellType()) {
      // case VTK_LINE: {
      //  Vector3<double> x1 = {mesh->GetPoint(p[0])[0], mesh->GetPoint(p[0])[1], mesh->GetPoint(p[0])[2]};
      //  Vector3<double> x2 = {mesh->GetPoint(p[1])[0], mesh->GetPoint(p[1])[1], mesh->GetPoint(p[1])[2]};
      //  double l = sqrt(cl->GetLength2()) / 1000.0; // mm -> m
      //  auto zeta(x)[] { (x - x1) / l };
      //  Vector3<double>[2] dNdX = {{},{}};

      //  for (int i = 0; i < 2; ++i) {
      //    for (int j = 0; j < 2; ++j) {
      //      ki[2*i+j] = (i==j ? 1.0 : -1.0) * sigma_i / l; // TODO XXX WRONG
      //      if (bi) {
      //        ke[2*i+j] = (i==j ? 1.0 : -1.0) * sigma_e / l;
      //      }
      //      mm[2*i+j] = (i == j ? 1.0/3.0 : 1.0/6.0) * l;
      //    }
      //  }
      // }
      case VTK_TRIANGLE: {
        Vector3<double> p0(mesh->GetPoint(p[0]));
        Vector3<double> p1(mesh->GetPoint(p[1]));
        Vector3<double> p2(mesh->GetPoint(p[2]));

        p0 /= 1000.;
        p1 /= 1000.;
        p2 /= 1000.;

        Vector3<double> p3 = 1./3. * (p0 + p1 + p2) + CrossProduct(p1 - p2, p2 - p0);
        Matrix3<double> vertices(p0 - p3, p1 - p3, p2 - p3);
        Matrix3<double> dNdX = vertices.GetInverse();

        Vector3<double> n = CrossProduct(p0 - p1, p0 - p2);
        double a          = .5 * n.Norm();
        n = n / n.Norm();

        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            Vector3<double> dNdXi = dNdX.GetRow(j) - (dNdX.GetRow(j) * n) * n;
            Vector3<double> dNdXj = dNdX.GetRow(k) - (dNdX.GetRow(k) * n) * n;

            if (hasIntra.m[m]) {
              ki[3 * j + k] = a * dNdXj * (sigma_i * dNdXi);
              mm[3 * j + k] = a * (k == j ? 1./6. : 1./12.);
            }

            if (bi) {
              // ke[3 * j + k] = a * sigma_e * d * e;
              ki[3 * j + k] = a * dNdXj * (sigma_e * dNdXi);
            }
          }
        }
      }
      break;
      case VTK_TETRA: {
        Matrix4<double> vertices;

        for (int j = 0; j < 4; j++) {
          vertices(0, j) = 1;
          for (int k = 1; k < 4; k++) {
            vertices(k, j) = mesh->GetPoint(p[j])[k-1] / 1000.0;  // mm -> m
          }
        }

        double J             = vertices.Det() / 6;
        Matrix4<double> dNdX = vertices.GetInverse();


        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            Vector3<double> dNdXi = {dNdX(i, 1), dNdX(i, 2), dNdX(i, 3)};
            Vector3<double> dNdXj = {dNdX(j, 1), dNdX(j, 2), dNdX(j, 3)};

            /* dNdX is const over the linear element, and the sum of
             * the weights w will be one, so we can just say,
             * regardless of the integration rule: */
            if (hasIntra.m[m]) {
              ki[4*i+j] = J * dNdXi * (sigma_i * dNdXj);
              /* For the mass matrix, we also have the analytical
               * solution, regardless of integration rule: */
              mm[4*i+j] = J * (i == j ? 1.0/10.0 : 1.0/20.0);
            }

            if (bi) {
              ke[4*i+j] = J * dNdXi * (sigma_e * dNdXj);
            }
          }
        }
      }
      break;
      default:
        ERR("Unsupported cell type (" << cl->GetCellType() << ") at index " << c, ERR_CELL);
    }  // switch

    if (bi) {
      ierr = MatSetValues(extraMatrix, np, p.data(), np, p.data(), ke.data(), ADD_VALUES);
      CHKE(ierr);
    }

    if (hasIntra.m[m]) {
      if (intra_is) {
        for (vtkIdType i = 0; i < np; ++i) {
          p[i] = globalToIntra[p[i]];
        }
      }

      ierr = MatSetValues(intraMatrix, np, p.data(), np, p.data(), ki.data(), ADD_VALUES);
      CHKE(ierr);
      ierr = MatSetValues(massMatrix, np, p.data(), np, p.data(), mm.data(), ADD_VALUES);
      CHKE(ierr);

      if (gauss) {
        for (int j = 0; j < gauss; ++j) {
          PetscInt gaussidx = intraCell * gauss + j;

          if (intraCell >= numIntraCells) {
            std::cerr << numIntraCells << " " << intraCell << "\n";
          }

          ierr = MatSetValues(interpolMatrix, 1, &gaussidx, np, p.data(),
                              gaussq->GetWeights(), INSERT_VALUES);

          ierr = MatSetValues(gaussMatrix, np, p.data(), 1, &gaussidx,
                              gaussq->GetParametricCoordinates(j).data(), INSERT_VALUES);
        }
        intraCell++;
      }
    }
  }

  ierr = MatAssemblyBegin(intraMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  ierr = MatAssemblyEnd(intraMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  ierr = MatAssemblyBegin(massMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  ierr = MatAssemblyEnd(massMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  if (bi) {
    ierr = MatAssemblyBegin(extraMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
    ierr = MatAssemblyEnd(extraMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  }
  if (gauss) {
    ierr = MatAssemblyBegin(interpolMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
    ierr = MatAssemblyEnd(interpolMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
    ierr = MatAssemblyBegin(gaussMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
    ierr = MatAssemblyEnd(gaussMatrix, MAT_FINAL_ASSEMBLY); CHKE(ierr);
  }

  if (fullLumping || lumping) {
    Vec lumpedMassVec = NULL;
    ierr = VecCreateSeq(PETSC_COMM_WORLD, numIntraNodes, &lumpedMassVec); CHKE(ierr);
    ierr = MatGetRowSum(massMatrix, lumpedMassVec); CHKE(ierr);

    if (fullLumping) {
      /* Invert lumped mass matrix */
      PetscScalar *pv = NULL;
      ierr = VecGetArray(lumpedMassVec, &pv); CHKE(ierr);
      for (PetscInt i = 0; i < numIntraNodes; ++i)
        pv[i] =  1.0 / pv[i];
      ierr = VecRestoreArray(lumpedMassVec, &pv); CHKE(ierr);
      ierr = VecAssemblyBegin(lumpedMassVec); CHKE(ierr);
      ierr = VecAssemblyEnd(lumpedMassVec); CHKE(ierr);
      ierr = MatDiagonalScale(intraMatrix, lumpedMassVec, NULL); CHKE(ierr);
      if (bi) {
        /* see bidomain eq. if we scale K_i, we must scale K_i+e the same
         * way. */
        Vec extraLumpedMass = NULL, extraLMsub = NULL;
        ierr = VecCreateSeq(PETSC_COMM_WORLD, numNodes, &extraLumpedMass); CHKE(ierr);
        ierr = VecSet(extraLumpedMass, 1.0); CHKE(ierr);
        ierr = VecGetSubVector(extraLumpedMass, intra_is, &extraLMsub); CHKE(ierr);
        ierr = VecCopy(lumpedMassVec, extraLMsub); CHKERRQ(ierr);
        ierr = VecRestoreSubVector(extraLumpedMass, intra_is, &extraLMsub); CHKE(ierr);

        ierr = MatDiagonalScale(extraMatrix, extraLumpedMass, NULL); CHKE(ierr);
        ierr = VecDestroy(&extraLumpedMass); CHKE(ierr);
      }
    } else {
      ierr = VecScale(lumpedMassVec, (1.0 - theta_mass)); CHKE(ierr);
      ierr = MatScale(massMatrix, theta_mass); CHKE(ierr);
      ierr = MatDiagonalSet(massMatrix, lumpedMassVec, ADD_VALUES); CHKE(ierr);
    }

    ierr = VecDestroy(&lumpedMassVec); CHKE(ierr);
  }

  if (coupled) {
    ERR("Fix this: Coupled only works if nodes are sorted, intra first!", 1337);
    Mat fullMat;

    /* Create full list of neighbours, first intra only, then both */
    std::for_each(intraNeighbours.begin(), intraNeighbours.end(),
                  [](PetscInt &i) {i *= 2;});  /* TODO does this work?! */
    intraNeighbours.resize(numIntraNodes + numNodes);
    for (PetscInt i = 0; i < numIntraNodes; ++i) {
      intraNeighbours[i+numIntraNodes] = intraNeighbours[i] + nodeNeighbours[i];
    }
    for (PetscInt i = numIntraNodes; i < numNodes; ++i) {
      intraNeighbours[i+numIntraNodes] = nodeNeighbours[i];
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numNodes + numIntraNodes,
                           numNodes + numIntraNodes, 0, intraNeighbours.data(), &fullMat); CHKE(ierr);

    for (PetscInt i = 0; i < numNodes; ++i) {
      PetscInt ncols;
      const PetscInt *cols;
      const PetscScalar *vals;

      ierr = MatGetRow(extraMatrix, i, &ncols, &cols, &vals); CHKE(ierr);

      PetscInt *newcols = new PetscInt[ncols];
      for (int j = 0; j < ncols; ++j) {
        newcols[j] = cols[i] + numIntraNodes;
      }

      PetscInt newrow = i + numIntraNodes;
      ierr = MatSetValues(fullMat, 1, &newrow, ncols, newcols, vals, INSERT_VALUES); CHKE(ierr);

      delete[] newcols;

      ierr = MatRestoreRow(extraMatrix, i, &ncols, &cols, &vals); CHKE(ierr);

      if (i < numIntraNodes) {
        ierr = MatGetRow(intraMatrix, i, &ncols, &cols, &vals); CHKE(ierr);

        PetscInt *newcols = new PetscInt[ncols];
        for (int j = 0; j < ncols; ++j) {
          newcols[j] = cols[i] + numIntraNodes;
        }

        /* Full Matrix Structure:
         *
         * [ [K_i] [K_i]       [0]         ]
         * [ [K_i] [K_i+e]     [K_e\K_i+e] ]
         * [ [0]   [K_e\K_i+e] [K_e]       ]
         *
         * */

        ierr = MatSetValues(fullMat, 1, &i, ncols, cols, vals, INSERT_VALUES);
        CHKE(ierr);
        ierr = MatSetValues(fullMat, 1, &i, ncols, newcols, vals, INSERT_VALUES);
        CHKE(ierr);

        // ierr = MatSetValues(fullMat, ncols, newcols, 1, &i, vals, INSERT_VALUES);
        // CHKE(ierr);
        PetscInt newrow = i + numIntraNodes;
        ierr = MatSetValues(fullMat, 1, &newrow, ncols, cols, vals, INSERT_VALUES);
        CHKE(ierr);

        delete[] newcols;

        ierr = MatRestoreRow(intraMatrix, i, &ncols, &cols, &vals);
        CHKE(ierr);
      }
    }


    ierr = MatScale(massMatrix, betaCm_dt); CHKE(ierr);

    /* Add betaCm/dt*M to fullMatrix */
    auto addMM = [ = ](Mat mmm) {
        for (PetscInt i = 0; i < numIntraNodes; ++i) {
          PetscErrorCode ierr;
          PetscInt ncols;
          const PetscInt *cols;
          const PetscScalar *m;
          ierr = MatGetRow(massMatrix, i, &ncols, &cols, &m); CHKE(ierr);

          ierr = MatSetValues(mmm, 1, &i, ncols, cols, m, ADD_VALUES);
          CHKE(ierr);

          ierr = MatRestoreRow(massMatrix, i, &ncols, &cols, &m); CHKE(ierr);
        }
      };

    if ((theta == .0) || (theta == 1.0)) {
      if (theta == 0.0) {
        ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
        ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
        ierr = MatScale(fullMat, -1.0); CHKE(ierr);
        ierr = MatScale(massMatrix, -1.0); CHKE(ierr);
      }

      addMM(fullMat);
      ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);

      // ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      // ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      saveMatrix(prefix+(theta == .0 ? "rhs.mat" : "lhs.mat"), fullMat);
      ierr = MatDestroy(&fullMat); CHKE(ierr);

      std::memset(intraNeighbours.data() + numIntraNodes, 0, numNodes * sizeof(PetscInt));
      ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numNodes + numIntraNodes,
                             numNodes + numIntraNodes, 0, intraNeighbours.data(), &fullMat); CHKE(ierr);

      addMM(fullMat);
      ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      saveMatrix((theta == .0 ? "lhs.mat" : "rhs.mat"), fullMat);
      ierr = MatDestroy(&fullMat); CHKE(ierr);
    } else {
      ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      Mat fullRHS;
      ierr = MatDuplicate(fullMat, MAT_COPY_VALUES, &fullRHS); CHKE(ierr);

      ierr = MatScale(fullMat, theta); CHKE(ierr);
      ierr = MatScale(fullRHS, -(1.0 - theta)); CHKE(ierr);

      addMM(fullMat);
      addMM(fullRHS);

      ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyBegin(fullRHS, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      ierr = MatAssemblyEnd(fullRHS, MAT_FINAL_ASSEMBLY); CHKE(ierr);
      saveMatrix(prefix+"lhs.mat", fullMat);
      saveMatrix(prefix+"rhs.mat", fullRHS);
      ierr = MatDestroy(&fullMat); CHKE(ierr);
      ierr = MatDestroy(&fullRHS); CHKE(ierr);
    }
  } else {
    saveMatrix(prefix+(mono ? "mat" : "intra.mat"), intraMatrix);
    if (!fullLumping) {
      saveMatrix(prefix+"mass.mat", massMatrix);
    }

    if (bi) {
      saveMatrix(prefix+(combined ? "i+e.mat" : "extra.mat"), extraMatrix);
    }
  }

  if (extraMatrix) {
    ierr = MatDestroy(&extraMatrix); CHKE(ierr);
  }
  ierr = MatDestroy(&intraMatrix); CHKE(ierr);
  ierr = MatDestroy(&massMatrix); CHKE(ierr);

  if (gauss) {
    saveMatrix(prefix+"n2e.mat", interpolMatrix);
    saveMatrix(prefix+"e2n.mat", gaussMatrix);

    ierr = MatDestroy(&interpolMatrix); CHKE(ierr);
    ierr = MatDestroy(&gaussMatrix); CHKE(ierr);
  }

  /* Save material vector */
  Vec matV;
  ierr = VecCreateSeq(PETSC_COMM_WORLD, numNodes, &matV); CHKE(ierr);
  PetscScalar *pv;
  ierr = VecGetArray(matV, &pv); CHKE(ierr);
  for (PetscInt i = 0; i < numNodes; ++i) {
    pv[i] = nodeMaterial[i];
  }
  ierr = VecRestoreArray(matV, &pv); CHKE(ierr);
  Vec intraMatV = matV;
  if (intra_is) {
    ierr = VecGetSubVector(matV, intra_is, &intraMatV); CHKE(ierr);
  }
  PetscViewer fd;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, (prefix+"vec"+(gzip ? ".gz" : "")).c_str(),
                               FILE_MODE_WRITE, &fd); CHKE(ierr);
  ierr = VecView(intraMatV, fd); CHKE(ierr);
  ierr = PetscViewerDestroy(&fd); CHKE(ierr);

  if (intra_is) {
    ierr = VecRestoreSubVector(matV, intra_is, &intraMatV); CHKE(ierr);
  }
  ierr = VecDestroy(&matV); CHKE(ierr);

  if (intra_is) {
    PetscViewer fd;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, (prefix+"is"+(gzip ? ".gz" : "")).c_str(),
                                 FILE_MODE_WRITE, &fd); CHKE(ierr);
    ierr = ISView(intra_is, fd); CHKE(ierr);
    ierr = PetscViewerDestroy(&fd); CHKE(ierr);
    ierr = ISDestroy(&intra_is); CHKE(ierr);
  }

  ierr = PetscFinalize(); CHKE(ierr);

  return 0;
}  // main

void saveMatrix(std::string const &fn, Mat m) {
  PetscViewer fd;
  PetscErrorCode ierr;

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, (fn+(gzip ? ".gz" : "")).c_str(),
                               FILE_MODE_WRITE, &fd); CHKE(ierr);
  ierr = MatView(m, fd); CHKE(ierr);
  ierr = PetscViewerDestroy(&fd); CHKE(ierr);
}

void printUsage(const char *name) {
  std::cerr << name << " <.vt?> [-o <prefix>] [-v]" << std::endl;
  std::cerr << "    [-bath <mask> [mat.def]]" << std::endl;
  std::cerr << "    [-mono <mat.def> [mat2.def]]" << std::endl;
  std::cerr << "    [-bi <intra.def> <extra.def>]" << std::endl;
  std::cerr << "    [-tissue <list>] [-freq <freq>]" << std::endl;
  std::cerr << "    [-combined] [-coupled <betaCm_dt>] [-gauss <n>]" << std::endl;
  std::cerr << "    [-fullLumping] [-massLumping <theta>]" << std::endl;
  std::cerr << "    [-iso] [-gzip] [-ignoreMissing]" << std::endl;
  std::cerr << name << " -help" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Built " __DATE__ " " __TIME__ << std::endl;
}

void printHelp() {
  std::cerr << "-v --- Enable additional informational (verbose) output.\n";
  std::cerr << std::endl;
  std::cerr << "-bath <mask> [mat.def] --- All elements matching <mask> will only be part of the extracellular\n";
  std::cerr << "       domain. Requires [-bi ...]. If [mat.def] is not given, <extra.def> is used.\n";
  std::cerr << std::endl;
  std::cerr << "-mono <mat.def> [mat2.def] --- Creates monodomain matrix. If [mat2.def] is given,\n";
  std::cerr << "       the conductivity will be the harmonic mean of the two given conductivities\n";
  std::cerr << "       for each cell (parallel circuit).\n";
  std::cerr << std::endl;
  std::cerr << "-bi <intra.def> <extra.def> --- Creates bidomain matrices with the intra- and\n";
  std::cerr << "       extracellular conductivities from the given material files.\n";
  std::cerr << std::endl;
  std::cerr << "-combined --- Creates bidomain matrices where intracellular matrix\n";
  std::cerr << "       is already added to extracellular matrix.\n";
  std::cerr << std::endl;
  std::cerr << "-coupled <betaCm_dt> [theta] --- Create fully coupled bidomain matrices (LHS and RHS)\n";
  std::cerr << "       <betaCm_dt> is (beta*C_m)/dt. theta is the time discretization weight\n";
  std::cerr << "       (0=explicit, 1=implicit, 0.5=Crank-Nicolson, default = 0).\n";
  std::cerr << std::endl;
  std::cerr << "-gauss <n> --- Create Gauss interpolation and integration matrices.\n";
  std::cerr << "       <n> is the number of Gauss points. Default integration rules are provided.\n";
  std::cerr << std::endl;
  std::cerr << "-fullLumping --- Use full mass lumping and multiply intra/mono matrix with\n";
  std::cerr << "       inverse of lumped mass matrix. Will not save mass matrix separately.\n";
  std::cerr << "       Does not apply with -coupled. Use -massLumping 0 instead.\n";
  std::cerr << std::endl;
  std::cerr << "-massLumping <theta> --- Use mass lumping. Saved mass matrix will be weighted sum:\n";
  std::cerr << "       M* = (theta)M + (1-theta)*M_lumped\n";
  std::cerr << "       Usually, you'd choose theta in [0,1]. Different choices may work, however.\n";
  std::cerr << std::endl;
  std::cerr << "-tissue <list> --- Comma-separated list of tissue classes in order of their,\n";
  std::cerr << "       precedence when creating the tissue vector. Highest precedence first.\n";
  std::cerr << "       Default: smaller numbers have higher precedence.\n";
  std::cerr << std::endl;
  std::cerr << "-freq <freq> --- Evaluate conductivities at frequency <freq> if supported.\n";
  std::cerr << std::endl;
  std::cerr << "-iso --- Ignore fiber orientation (isotropic case).\n";
  std::cerr << std::endl;
  std::cerr << "-ignoreMissing --- Ignore missing intra conductivities, i.e., set them to 0.\n";
  std::cerr << "       Only recommended to use with -bi.\n";
  std::cerr << std::endl;
  std::cerr << "-bds <value> --- CV stabilization factor beta*Delta_s.\n";
  std::cerr << std::endl;
  std::cerr << "-gzip --- Compress output files.\n";
  std::cerr << std::endl;
  std::cerr << "-o <prefix> --- Prefix for saving the resulting files. Defaults to the name\n";
  std::cerr << "       of the input geometry sans extension. Depending on the given options,\n";
  std::cerr << "       the following files may be created:\n";
  std::cerr << "           <prefix>.vec - Material class vector.\n";
  std::cerr << "           <prefix>.mat - Monodomain matrix.\n";
  std::cerr << "           <prefix>.mass.mat - Mass matrix.\n";
  std::cerr << "           <prefix>.intra.mat - Bidomain intra matrix.\n";
  std::cerr << "           <prefix>.extra.mat - Bidomain extra matrix.\n";
  std::cerr << "           <prefix>.i+e.mat - Bidomain combined matrix.\n";
  std::cerr << "           <prefix>.n2e.mat - Node-to-element Gauss point (interpolation) matrix.\n";
  std::cerr << "           <prefix>.e2n.mat - Element-to-node Gauss point (integration) matrix.\n";
  std::cerr << "           <prefix>.lhs.mat - Coupled LHS matrix.\n";
  std::cerr << "           <prefix>.rhs.mat - Coupled RHS matrix.\n";
  std::cerr << "           <prefix>.is - Index set of intracellular domain.\n";
}  // printHelp
