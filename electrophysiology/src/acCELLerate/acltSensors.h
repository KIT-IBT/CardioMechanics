/*
 * File: acltSensors.h
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



#ifndef ACLTSENSORS_H
#define ACLTSENSORS_H

#include "acltTime.h"
#include <list>
#include <string>

const int MaxSensors = 1024;

typedef kaBaseException ACLTSensorError;

enum SensorType { ST_II, ST_IE, ST_IF, ST_UI, ST_UE, ST_UEM, ST_UF, ST_PARA };


//! Base class for sensor in the bidomain case


class SingleSensor {
public:
    SingleSensor() {fp = NULL;}
    
    ~SingleSensor() {
        if (fp)
            fclose(fp);
    }
    
    PetscInt SensorIndex;
    std::string filename;
    FILE *fp;
    acltTime beginsave, dtsave;
    SensorType cflag;
    string SensorParameterType;
    double tmp;
    inline void ScanLine(const char *);
    
    void Init() {
        fp = fopen(filename.c_str(), "w");
        if (!fp)
            throw ACLTSensorError("Cannot open ACLTSSensor file '%s'", filename.c_str());
    }
    
    void Save(acltTime time, double val) {
        fprintf(fp, "%s %lf\n", time.toStr().c_str(), val);
        fflush(fp);
    }
};  // class SingleSensor

void SingleSensor::ScanLine(const char *buf) {
    char type[256], begin[64], dt[64];
    char fn[256];
    int  rc;
    
    rc              = sscanf(buf, PETSCINT_FORMAT " %255s %63s %63s %255s", &SensorIndex, fn, begin, dt, type);
    this->beginsave = begin;
    this->dtsave    = dt;
    
#if KADEBUG
    fprintf(stderr, "ScanLine '" PETSCINT_FORMAT "' '%s' '%s' '%s' '%s' \n",
            SensorIndex, fn, beginsave.toStr().c_str(), dtsave.toStr().c_str(), type);
#endif  // if KADEBUG
    
    if (rc != 5)
        throw ACLTSensorError("ACLTSSensor line includes not enough words");
    
    filename = fn;
    
    SensorParameterType = type;
    if (!strcmp(type, "Ie"))
        cflag = ST_IE;
    else if (!strcmp(type, "Ii"))
        cflag = ST_II;
    else if (!strcmp(type, "If"))
        cflag = ST_IF;
    else if (!strcmp(type, "Ve"))
        cflag = ST_UE;
    else if (!strcmp(type, "Vem"))
        cflag = ST_UEM;
    else if (!strcmp(type, "Vm"))
        cflag = ST_UI;
    else if (!strcmp(type, "Vf"))
        cflag = ST_UF;
    else
        cflag = ST_PARA;
}  // SingleSensor::ScanLine

class ACLTSensors {
protected:
    int numOfSensors;
    std::list<SingleSensor> SC;
    
public:
    inline ACLTSensors();
    virtual ~ACLTSensors() {}
    
    inline void Load(const char *name);
    inline void Parse(PetscInt VectorSize);
    inline void Init();
    
    inline int GetAnz() const {return numOfSensors;}
};


ACLTSensors::ACLTSensors() {
    numOfSensors = 0;
}

void ACLTSensors::Load(const char *name) {
#if KADEBUG
    fprintf(stderr, "ACLTSSensors::Load %s\n", name);
#endif  // if KADEBUG
    
    numOfSensors = 0;
    
    ifstream fp(name);
    if (fp.good()) {
        std::string line;
        while (std::getline(fp, line)) {
            SingleSensor sc;
            sc.ScanLine(line.c_str());
            SC.push_back(sc);
            numOfSensors++;
        }
    } else {
        throw ACLTSensorError("Unable to read sensors from file %s", name);
    }
}

/* This function used to parse the sensorfile that had been loaded into memory by Load(). Now,
 * the parsing is already done in Load() and Parse() more or less just verifies.
 * TODO: Merge Load() and Parse(). */
void ACLTSensors::Parse(PetscInt VectorSize) {
#if KADEBUG
    fprintf(stderr, "ACLTSensors::Parse(%ld)\n", (long)VectorSize);
#endif  // if KADEBUG
    
    for (auto sc = SC.cbegin(); sc != SC.cend(); ++sc) {
        if ((sc->SensorIndex < 0) || (sc->SensorIndex >= VectorSize))
            throw ACLTSensorError("Sensor at " PETSCINT_FORMAT " is outside of vector of length " PETSCINT_FORMAT,
                                  sc->SensorIndex, VectorSize);
    }
    
#if KADEBUG
    fprintf(stderr, "ACLTSSensors::Parse finished Sensors# %d\n", numOfSensors);
#endif  // if KADEBUG
}

void ACLTSensors::Init() {
#if KADEBUG
    fprintf(stderr, "ACLTSensors::Init()\n");
#endif  // if KADEBUG
    
    for (auto sc = SC.begin(); sc != SC.end(); ++sc)
        sc->Init();
}

#endif  // ifndef ACLTSENSORS_H
