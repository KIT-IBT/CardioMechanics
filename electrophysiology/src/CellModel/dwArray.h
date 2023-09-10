/*
 *  dwArray.h
 *  dwIniFile
 *
 *  Created by Daniel Weiss on 12.08.04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef DWARRAY_H
#define DWARRAY_H

#include <kaMachineOS.h>
#include <kaExceptions.h>

template<class T>  class dwArray {
 private:
  enum DataType {
    DT_FLOAT, DT_P_FLOAT, DT_INT, DT_P_INT, DT_UNSIGNED_CHAR, DT_INI_ENTRY, DT_STRING, DT_HLFInfo, DT_DUMMY
  };
  DataType dt;
  int start;
  int numEntries;
  T *v;
  T *pv;
  void allocArray(uint32_t start, uint32_t numEntries);

  inline void free() {
    // printf("freeing array of %s with %i Bytes ...\n",DTName(dt).c_str(),numEntries*(int)sizeof(T));
    delete[]v;
    pv = NULL;
  }

  string DTName(DataType dt);

 public:
  dwArray(int size) {allocArray(0, size);}

  dwArray(uint32_t from, uint32_t to) {allocArray(from, to-from+1);}

  dwArray() {allocArray(0, 0);}

  ~dwArray() {free();}

  inline int lBound() {return start;}

  inline int uBound() {return start+numEntries-1;}

  inline void fillWidth(T value) {for (int i = 0; i < numEntries; i++) *(pv+i) = value; }

  string arrayInformation();

  inline T & operator[](uint32_t index);
  int        addEntry(T vNewEntry, int afterIndex = -69);
  int        deleteEntry(int index);
  void       alloc(uint32_t start, uint32_t numEntries);
};  // class dwArray

template<class T> string dwArray<T>::DTName(DataType dt) {
  switch (dt) {
    case DT_FLOAT:
      return "<float>";

    case DT_P_FLOAT:
      return "<*float>";

    case DT_INT:
      return "<int>";

    case DT_P_INT:
      return "<*int>";

    case DT_UNSIGNED_CHAR:
      return "<unsigned char>";

    case DT_STRING:
      return "<string>";

    case DT_HLFInfo:
      return "<HLFInfo>";

    case DT_INI_ENTRY:
      return "<iniFile Entry>";

    default:
      return "<unknown type>";

      break;
  }  // switch
}  // >::DTName

template<class T> string dwArray<T>::arrayInformation() {
  char out[1024]    = "Array of ";
  char buf[1024/16] = "";

  sprintf(buf, "%s %s (%i ... %i) with %i entries:\n", buf, DTName(dt).c_str(), start, start+numEntries-1, numEntries);
  strcat(out, buf);

  for (int i = 0; i < numEntries; i++) {
    switch (dt) {
      case DT_FLOAT:
        sprintf(buf, "[%i]:\t%f\n", i+start, *(pv+i));
        break;
      case DT_INT:
        sprintf(buf, "[%i]:\t%i\n", i+start, *(pv+i));
        break;
      case DT_P_FLOAT:
      case DT_P_INT:
        sprintf(buf, "[%i]:\t%p\n", i+start, *(pv+i));
        break;
      case DT_UNSIGNED_CHAR:
        sprintf(buf, "[%i]:\t%c\n", i+start, *(pv+i));
        break;
      default:
        sprintf(buf, "unknown typeid '%s' - can't show formated values ...\n", typeid(*v).name());
    }
    strcat(out, buf);
  }
  string a = out;
  return a;
}  // >::arrayInformation

template<class T> inline T & dwArray<T>::operator[](uint32_t index) {
  if (((int)index >= start) && ((int)index <= (start + numEntries - 1)))
    return (T &)(*(pv+index-start));
  else
    throw kaBaseException("index %i out of index range (%i ... %i)\n", index, start, start+numEntries-1);
}

template<class T> void dwArray<T>::allocArray(uint32_t from, uint32_t size) {
  if (!(strcmp(typeid(T).name(), "i"))) {
    dt = DT_INT;
  } else if (!(strcmp(typeid(T).name(), "Pi"))) {
    dt = DT_P_INT;
  } else if (!(strcmp(typeid(T).name(), "f"))) {
    dt = DT_FLOAT;
  } else if (!(strcmp(typeid(T).name(), "Pf"))) {
    dt = DT_P_FLOAT;
  } else if (!(strcmp(typeid(T).name(), "h"))) {
    dt = DT_UNSIGNED_CHAR;
  } else if (!(strcmp(typeid(T).name(), "8iniEntry"))) {
    dt = DT_INI_ENTRY;
  } else if (!(strcmp(typeid(T).name(), "7HLFInfo"))) {
    dt = DT_HLFInfo;
  } else if (!(strcmp(typeid(T).name(), "Ss"))) {
    dt = DT_STRING;
  } else {
    // printf("id:(%s)\n",typeid(T).name());
    dt = DT_DUMMY;
  }

  // printf("allocating array of %s [%i ... %i] = %i * %i Bytes = %i Bytes
  // ...\n",DTName(dt).c_str(),from,from+size-1,size,(int)sizeof(T),size*(int)sizeof(T));
  start      = from;
  numEntries = size;
  v          = new T[size];
  pv         = v;
}  // >::allocArray

template<class T> int dwArray<T>::addEntry(T vNewEntry, int afterIndex) {
  T *tmp  = new T[numEntries+1];
  T *ptmp = tmp;

  if ((afterIndex > numEntries+start-1) || (afterIndex == -69))
    afterIndex = numEntries+start-1;
  else if (afterIndex < start)
    afterIndex = start-1;

  afterIndex = afterIndex-start+1;

  for (int i = 0; i < afterIndex; i++)
    *(ptmp+i) = *(pv+i);
  *(ptmp+afterIndex) = vNewEntry;
  for (int i = afterIndex+1; i < numEntries+1; i++)
    *(ptmp+i) = *(pv+i-1);

  numEntries++;
  v  = new T[numEntries];
  pv = v;
  for (int i = 0; i < numEntries; i++)
    *(pv+i) = *(ptmp+i);

  delete[] tmp;
  ptmp = NULL;
  return afterIndex+start;
}  // >::addEntry

template<class T> int dwArray<T>::deleteEntry(int index) {
  T *tmp  = new T[numEntries-1];
  T *ptmp = tmp;

  if (index > numEntries+start-1)
    index = numEntries+start-1;
  else if (index < start)
    index = start;

  index = index-start;

  printf("delete_index = %i, numEntries = %i\n", index, numEntries);

  for (int i = 0; i < index; i++)
    *(ptmp+i) = *(pv+i);
  for (int i = index; i < numEntries-1; i++)
    *(ptmp+i) = *(pv+i+1);

  numEntries--;
  v  = new T[numEntries];
  pv = v;
  for (int i = 0; i < numEntries; i++)
    *(pv+i) = *(ptmp+i);

  delete[] tmp;
  ptmp = NULL;
  return index+start;
}  // >::deleteEntry

template<class T>
void dwArray<T>::alloc(uint32_t start, uint32_t numEntries) {
  this->allocArray(start, numEntries);
}

#endif  // ifndef DWARRAY_H
