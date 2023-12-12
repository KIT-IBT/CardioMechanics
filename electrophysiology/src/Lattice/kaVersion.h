/*
 * File: kaVersion.h
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


#ifndef KAVERSION_H
#define KAVERSION_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <kaExceptions.h>

// by Ouss
// I find a static Array is a bad programming techinque,
// a better alternative will be a std list of anOption.
const int MAXOPTIONS = 100;

enum optionTypeEnum {
  OT_optional = 1,
  OT_required = 2
};

class anOption {
 public:
  optionTypeEnum optionType;
  std::string key;
  std::string values;
  std::string defaultValue;
  std::string desc;
};

class kaVersion {
 public:
  kaVersion(int major, int minor, int revision, int argc, char *const argv[]);
  kaVersion(int major, int minor, int revision, int year, int month, int day, int argc, char *const argv[]);
  void init(int major, int minor, int revision, int year, int month, int day, int argc, char *const argv[]);
  void addOption(std::string key, std::string values = "", optionTypeEnum OT = OT_required,
                 std::string defaultValue            = "", std::string desc = "");

  void addOptionalOption(std::string key, std::string values = "", std::string defaultValue = "",
                         std::string desc                    = "") {
    addOption(key, values, OT_optional, defaultValue, desc);
  }

  void addRequiredOption(std::string key, std::string values = "", std::string defaultValue = "",
                         std::string desc                    = "") {
    addOption(key, values, OT_required, defaultValue, desc);
  }

  void addVersionDescription(std::string descriptions) {versiondescription = descriptions;}

  void printHelp(void);

 protected:
  anOption o[MAXOPTIONS];
  int major, minor, revision;
  int year, month, day;
  int lastOption;
  std::string toolName;
  std::string getVersionString(void);
  std::string getFullVersionString(void);
  std::string versiondescription;
}; // class kaVersion

#endif  // ifndef KAVERSION_H
