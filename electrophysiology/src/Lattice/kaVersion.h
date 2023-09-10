/**@file kaVersion.h
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

#ifndef KAVERSION_H
#define KAVERSION_H
/* -------------------------------------------------------

    kaVersion.h

    Ver. 1.1.0

    Created:       Daniel Weiss      (24.07.2008)
    Last modified: Frank M. Weber    (11.11.2008)

    Institute of Biomedical Engineering
    Universitaet Karlsruhe (TH)

    http://www.ibt.uni-karlsruhe.de

    Copyright 2000-2008 - All rights reserved.

        HISTORY:

        v1.1.0 (2008-11-11) (fmw): Added support for date information,
                                                           use constructor with (major, minor, rev, year, month, day)

        v1.0.1 (2008-10-29) (oj): <Document changes here>

   ------------------------------------------------------ */
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
