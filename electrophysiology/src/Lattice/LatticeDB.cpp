/*
 * File: LatticeDB.cpp
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


#include <LatticeDB.h>

LatticeDB::LatticeDB(char *name) {
  for (int i = 0; i < MaxLatticesInDB; i++) {
    m[i].New(4);
    mi[i].New(4);
    latname[i] = NULL;
  }
  dbname = strdup(name);
}

LatticeDB::~LatticeDB() {
  DeleteEntries();

  delete dbname;
}

void LatticeDB::DeleteEntries() {
  for (int i = 0; i < MaxLatticesInDB; i++) {
    delete latname[i];
    latname[i] = NULL;
  }
}

int LatticeDB::FindEntry(char *name) {
  int i = 0;

  for (i = 0; i < MaxLatticesInDB; i++)
    if (!latname[i]) {
      if (!name)
        break;
    } else if (name) {
      if (!strcmp(name, latname[i]))
        break;
    }

  return i == MaxLatticesInDB ? -1 : i;
}

int LatticeDB::FindLattice(kaPointDouble &P, int offset) {
  cerr << P << endl;
  int i;
  for (i = offset+1; i < MaxLatticesInDB; i++)
    if (latname[i]) {
      kaPointDouble L = P;
      L.XMatrixN(mi[i]);
      cerr << L << endl;
      if ((L.x >= 0) && (L.y >= 0) && (L.z >= 0)) {
        if ((L.x < xLattice[i]) && (L.y < yLattice[i]) && (L.z < zLattice[i]))
          break;
      }
    }

  return i == MaxLatticesInDB ? -1 : i;
}

int LatticeDB::Save() {
  int err = 0;

  FILE *fp = fopen(dbname, "w");

  for (int i = 0; i < MaxLatticesInDB; i++)
    if (latname[i]) {
      fprintf(fp, "%s\n", latname[i]);
      m[i].Save(fp);
      fwrite(&xLattice[i], sizeof(xLattice[0]), 1, fp);
      fwrite(&yLattice[i], sizeof(yLattice[0]), 1, fp);
      fwrite(&zLattice[i], sizeof(zLattice[0]), 1, fp);
    } else {
      fputs("", fp);
    }

  fclose(fp);

  return err;
}

int LatticeDB::Restore() {
#if KALATTICEDEBUG
  fprintf(stderr, "LatticeDB::Restore()\n");
#endif  // if KALATTICEDEBUG
  int err = 0;

  DeleteEntries();

  FILE *fp = fopen(dbname, "r");
  if (!fp)
    return 1;

  for (int i = 0; i < MaxLatticesInDB; i++) {
    char buf[1024];
    fgets(buf, sizeof(buf), fp);
    if (feof(fp))
      return -1;

    if (buf[0]) {
      latname[i] = strdup(buf);
      m[i].Restore(fp);
      mi[i] = m[i];

      //                        m[i].Print();
      mi[i].Inverse();
      fread(&xLattice[i], sizeof(xLattice[0]), 1, fp);
      fread(&yLattice[i], sizeof(yLattice[0]), 1, fp);
      fread(&zLattice[i], sizeof(zLattice[0]), 1, fp);
    }
  }

  fclose(fp);

  return err;
}  // LatticeDB::Restore

int LatticeDB::TransformkaPointDouble(int i, IndexType x, IndexType y, IndexType z, kaPointDouble &P) {
  P.x = x; P.y = y; P.z = z;
  P.XMatrixN(m[i]);
  return 1;
}

int LatticeDB::InsertEntry(char *name) {
#if KALATTICEDEBUG
  fprintf(stderr, "LatticeDB::InsertEntry(%s)\n", name);
#endif  // if KALATTICEDEBUG
  int i = FindEntry(name);
  if (i < 0)
    i = FindEntry();

  InsertEntry(i, name);
  return 1;
}

int LatticeDB::DeleteEntry(char *name) {
  int i = FindEntry(name);

  delete latname[i];
  latname[i] = NULL;
  return 1;
}

int LatticeDB::InsertEntry(int i, char *name) {
#if KALATTICEDEBUG
  fprintf(stderr, "LatticeDB::InsertEntry(%d, %s)\n", i, name);
#endif  // if KALATTICEDEBUG
  assert(i >= 0 && i < MaxLatticesInDB);

  kaLatticeCreate cls;
  strcpy(cls.Name, name);
  delete latname[i];
  latname[i] = strdup(name);

  kaLatticeSharedHeader sh(cls);
  m[i]        = sh.lhMatrix();
  xLattice[i] = sh.lhWidth();
  yLattice[i] = sh.lhHeight();
  zLattice[i] = sh.lhDepth();

  m[i].Print();

  mi[i] = m[i];
  mi[i].Inverse();

  mi[i].Print();

  return 0;
}  // LatticeDB::InsertEntry

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "%s DataBase\n", argv[0]);
    fprintf(stderr, "\t\t[-ins kaLattice... ]\n");
    fprintf(stderr, "\t\t[-del kaLattice... ]\n");
    fprintf(stderr, "\t\t[-point x y z]\n");
    fprintf(stderr, "\t\t[-save]\n");
    exit(-1);
  }

  bool save = false;

  LatticeDB LDB(argv[1]);
  LDB.Restore();

  for (int i = 2; i < argc; i++) {
    if (((strcmp(argv[i], "-ins") == 0) || (strcmp(argv[i], "-insert") == 0)) && (i+1 < argc)) {
      while (++i < argc) {
        if (argv[i][0] == '-')
          break;
        LDB.InsertEntry(argv[i]);
      }
      i--;
    } else if (((strcmp(argv[i], "-del") == 0) || (strcmp(argv[i], "-delete") == 0)) && (i+1 < argc)) {
      while (++i < argc) {
        if (argv[i][0] == '-')
          break;
        LDB.DeleteEntry(argv[i]);
      }
      i--;
    } else if ((strcmp(argv[i], "-point") == 0) && (i+3 < argc)) {
      kaPointDouble P;
      P.x = atof(argv[++i]);
      P.y = atof(argv[++i]);
      P.z = atof(argv[++i]);
      for (int o = LDB.FindLattice(P); o >= 0; o = LDB.FindLattice(P, o))
        printf("%s\n", LDB.GetName(o));
    } else if (strcmp(argv[i], "-save") == 0) {
      save = true;
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[i]);
      exit(-1);
    }
  }

  if (save)
    LDB.Save();
}  // main
