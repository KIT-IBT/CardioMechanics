/*
**      Name
**              LatticeDB.h
**
**      Usage
**              LatticeDB.cpp
**
**      Remark
**              Verwaltung von zusammengehoerigen Lattices
**
**
**      History
**              01.02.98 -fs creation
**
**
**  created at IBT - Universitaet Karlsruhe
*/


#ifndef LATTICEDB_H
#define LATTICEDB_H

#include <kaLattice.h>


const int MaxLatticesInDB = 256;


class LatticeDB {
  char *dbname;

  char *latname[MaxLatticesInDB];
  kaSharedMatrixN<double> m[MaxLatticesInDB];
  kaSharedMatrixN<double> mi[MaxLatticesInDB];
  IndexType xLattice[MaxLatticesInDB], yLattice[MaxLatticesInDB], zLattice[MaxLatticesInDB];
  int  InsertEntry(int i, char *name);
  void DeleteEntries();
  int  FindEntry(char *name = NULL);

 public:
  LatticeDB(char *name);
  ~LatticeDB();
  int InsertEntry(char *name);
  int DeleteEntry(char *name);
  int Restore();
  int Save();
  int FindLattice(kaPointDouble &, int offset = -1);
  int TransformkaPointDouble(int i, IndexType x, IndexType y, IndexType z, kaPointDouble &P);

  char *GetName(int i) {return latname[i];}
};

#endif  // ifndef LATTICEDB_H
