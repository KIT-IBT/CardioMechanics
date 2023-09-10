/**@file CellModelValues.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef CELLMODELVALUES_H
#define CELLMODELVALUES_H

#include <kaMachineOS.h>
#include <kaExceptions.h>
#include <kaIPC.h>
#include <kaRootDir.h>
#include <sstream>
#include <vector>

template<class T>
class CellModelMembers {
 public:
  CellModelMembers(void) {}

  virtual ~CellModelMembers(void) {}

  virtual int GetSize(void) = 0;
  virtual T  *GetBase(void) = 0;

  inline void Write(FILE *writeto) {
    fwrite(GetBase(), GetSize(), 1, writeto);
  }

  inline void Read(FILE *readfrom) {
    fread(GetBase(), GetSize(), 1, readfrom);
  }

  inline void Set(CellModelMembers &from) {
    memcpy(GetBase(), from.GetBase(), from.GetSize());
  }
};


class CellModelConstant {
 public:
  char Name[64];
  double Value;

  CellModelConstant() {
    Name[0] = 0;
    Value   = 0.0;
  }

  CellModelConstant & operator=(const CellModelConstant &t) {
    strcpy(Name, t.Name);
    Value = t.Value;
    return *this;
  }

  void write(ostream &fp) {
    fp<<Name<<" "<<Value<<endl;
  }

  void read(istream &fp) {
    fp>>Name>>Value;
    while (!fp.eof()) // read over potential comment
      if (fp.get() == '\n')
        break;
  }
};  // class CellModelConstant


static const char *MAGIC_CELL_ID       = "CELLMODEL*IBT*KA";
static const int32_t MAGIC_CELL_ID_LEN = 16;
static const int32_t OVERLAP_LEN       = 16;

enum ModelType { MT_ELPHY, MT_FORCE };


class CellModelValues {
 public:
  char ProgramID[32];
  char ModelName[32];
  char Overlapstring[128];
  int numElements;
  CellModelConstant mc[256];
  char comment[500];


  CellModelValues() {
    memset(comment, 0, sizeof(comment));
  }

  void init(char *ID, char *Name, int Elements) {
    strcpy(ProgramID, ID);
    strcpy(ModelName, Name);
    numElements = Elements;
    memset(comment, 0, sizeof(comment));
  }

  void init(char *ID, char *Name, int Elements, CellModelConstant *mcon, char *text = NULL) {
    init(ID, Name, Elements);
    for (int i = 0; i < numElements; i++)
      mc[i] = mcon[i];
    if (text)
      strcpy(comment, text);
  }

  void setConstant(int num, const char *Name, double Value) {
    strcpy(mc[num].Name, Name);
    mc[num].Value = Value;
  }

  void setOLString(char *ol) {
    strcpy(Overlapstring, ol);
  }

  //  void setComment( std::string text){
  //      comment = text.c_str();
  //  }

  char *getProgramID() {return ProgramID;}

  char *getModelName() {return ModelName;}

  int getnumElements() {return numElements;}

  const char *getconstName(int i) {return mc[i].Name;}

  double getconstValue(int i) {return mc[i].Value;}

  char *getComment() {return comment;}

  char *getOLString() {return Overlapstring;}

  /*  ostream& operator<<(ostream& _o, this)    */
  /*     {  */
  /*       _o << ProgramID << endl; */
  /*       _o << ModelName << endl; */
  /*       _o << numElements << endl; */

  /*       for (int i=0; i<numElements; i++) */
  /*    _o << mc[i].Name <<" "<< mc[i].Value <<endl; */

  /*       return _o;  */
  /*     } ; */

  void HeaderCheck(const char *filename, const char *modeltype = 0) {
    ifstream fp(filename);

    if (!fp)
      throw kaBaseException("CellModelValues::HeaderCheck - Unable to open file %s!", filename);
    HeaderCheck(fp, modeltype);
    fp.close();
  }

  void HeaderCheck(ifstream &filename, const char *modeltype = 0) {
    filename.seekg(0);
    std::stringstream ss;
    std::string str;
    int i = 0;
    while (i < 3 && !filename.eof()) {
      std::getline(filename, str);
      std::size_t found = str.find_first_not_of(" \t");
      if (found != std::string::npos)
        if (str.at(found) != '#') {
          ss << str << " ";
          ++i;
        }
    }
    ss >> str;
    if (str.size() < 32)
      strcpy(ProgramID, str.c_str());
    else
      throw kaBaseException("CellModelValues::HeaderCheck - Magic id must not be longer than 32 characters");

    if (strncmp(ProgramID, MAGIC_CELL_ID, MAGIC_CELL_ID_LEN))
      throw kaBaseException("CellModelValues::HeaderCheck - Header is incorrect!\n%s != %s", ProgramID, MAGIC_CELL_ID);

    ss >> str;
    if (str.size() < 32)
      strcpy(ModelName, str.c_str());
    else
      throw kaBaseException("CellModelValues::HeaderCheck - Model name must not be longer than 32 characters");
    if (modeltype && strcmp(ModelName, modeltype))
      throw kaBaseException("CellModelValues::HeaderCheck - Model type is incorrect!\n%s != %s", ModelName, modeltype);

    ss >> numElements;
    if ((numElements < 0) || (numElements > 256) )
      throw kaBaseException("CellModelValues::HeaderCheck - Number of parameters must be between 0 and 256");
  }  // HeaderCheck

  void Load(const char *filename, ModelType Modtyp = MT_ELPHY) {
    ifstream fp(filename);

    if (!fp)
      throw kaBaseException("CellModelValues::Load - Unable to open file!");
    Load(fp, Modtyp);
    fp.close();
  }

  void Load(ifstream &filename, ModelType Modtyp = MT_ELPHY) {
    filename.seekg(0);
    int i = 0;
    std::string str;
    while (i < 3 && !filename.eof()) {
      std::getline(filename, str);
      std::size_t found = str.find_first_not_of(" \t");
      if (found != std::string::npos)
        if (str.at(found) != '#')
          i++;
    }

    if (Modtyp == MT_FORCE)
      filename.getline(Overlapstring, sizeof(Overlapstring));

    for (int i = 0; i < numElements; ++i) {
      std::streampos pos = filename.tellg();
      while (std::getline(filename, str)) {
        std::size_t found = str.find_first_not_of(" \t");
        if (found != std::string::npos) {
          if (str.at(found) == '#')
            pos = filename.tellg();
          else
            break;
        } else {
          pos = filename.tellg();
        }
      }
      filename.clear();
      filename.seekg(pos);

      mc[i].read(filename);
    }
    memset(comment, 0, sizeof(comment));
    if (!filename.eof())
      filename.getline(comment, sizeof(comment), '\f');
  }  // Load

  bool Save(const char *filename, const char *OLData) {
    ofstream fp(filename);

    if (!fp) {
      return false;

      throw kaBaseException("CellModelValues::Save - Unable to open file!");
    }
    return Save(fp, OLData);
  }

  bool Save(ofstream &fp, const char *OLData) {
    bool rc = true;

    fp<<ProgramID<<endl;
    fp<<ModelName<<endl;
    fp<<numElements<<endl;
    fp<<OLData<<endl;
    for (int i = 0; i < numElements; i++)
      mc[i].write(fp);
    if (!comment[0]) {
      // fp<<"*"<<endl;
    } else {
      fp<<comment;
      fp<<'\f';
    }
    fp.close();
    return rc;
  }

  bool Save(const char *filename) {
    ofstream fp(filename);

    if (!fp) {
      return false;

      throw kaBaseException("CellModelValues::Save - Unable to open file!");
    }
    return Save(fp);
  }

  bool Save(ofstream &fp) {
    bool rc = true;

    fp<<ProgramID<<endl;
    fp<<ModelName<<endl;
    fp<<numElements<<endl;
    for (int i = 0; i < numElements; i++)
      mc[i].write(fp);
    if (!comment[0]) {
      // fp<<"*"<<endl;
    } else {
      fp<<comment;
      fp<<'\f';
    }
    fp.close();
    return rc;
  }
};  // class CellModelValues


#endif  // ifndef CELLMODELVALUES_H
