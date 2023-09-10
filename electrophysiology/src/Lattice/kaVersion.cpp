/* -------------------------------------------------------

    kaVersion.cpp

    Ver. 1.1.1

    Created:       Daniel Weiss      (24.07.2008)
    Last modified: Oussama Jarrousse    (05.12.2008)

    Institute of Biomedical Engineering
    Universitaet Karlsruhe (TH)

    http://www.ibt.uni-karlsruhe.de

    Copyright 2000-2008 - All rights reserved.

        HISTORY:

        v1.1.2 (2008-12-05) (oj116): Bugfix: in init(), the function doesn't write "done 1" to the error stream anymore.

        v1.1.1 (2008-11-17) (fmw): getVersionString now only returns version w/o date
                                                           getFullVersionString returns version + date + description

        v1.1.0 (2008-11-11) (fmw): Added support for date information,
                                                           use constructor with (major, minor, rev, year, month, day)

        v1.0.1 (2008-10-29) (oj116): Annotations!

   ------------------------------------------------------ */

#include "kaVersion.h"

kaVersion::kaVersion(int _major, int _minor, int _revision, int _year, int _month, int _day, int argc,
                     char *const argv[]) {
  init(_major, _minor, _revision, _year, _month, _day, argc, argv);
}

kaVersion::kaVersion(int _major, int _minor, int _revision, int argc, char *const argv[]) {
  init(_major, _minor, _revision, -1, -1, -1, argc, argv);
}

void kaVersion::init(int _major, int _minor, int _revision, int _year, int _month, int _day, int argc,
                     char *const argv[]) {
  major    = _major;
  minor    = _minor;
  revision = _revision;
  year     = _year;
  month    = _month;
  day      = _day;

  lastOption = 0;
  toolName   = argv[0];
  for (int i = 0; i < argc; i++) {
    // std::cout<<i<<": "<<argv[i]<<endl;
    if (strcmp(argv[i], "-version") == 0) {
      std::cout<<getVersionString();
      exit(0);
    }
  }

  // cerr << "done 1" << endl;
}

std::string kaVersion::getVersionString(void) {
  // What!!! 50 how about using the std::string << function :-) the std::string is a dynamic char array anyway.
  char ret[100];

  sprintf(ret, "%i.%i.%i", major, minor, revision);
  return ret;
}

std::string kaVersion::getFullVersionString(void) {
  // What!!! 50 how about using the std::string << function :-) the std::string is a dynamic char array anyway.
  char ret[1000];

  sprintf(ret, "%i.%i.%i", major, minor, revision);

  // Only if a date was set, add date information to the string
  if ( (year > -1) && (month > -1) && (day > -1) ) {
    sprintf(ret, "%s (%04d-%02d-%02d)", ret, year, month, day);
  }
  if (!versiondescription.empty()) {
    sprintf(ret, "%s %s", ret, versiondescription.c_str());
  }
  return ret;
}

void kaVersion::printHelp(void) {
  std::cerr<<"\n";
  std::cerr<<toolName;
  std::cerr<<" ";
  std::cerr<<"v"<<getFullVersionString().c_str();
  std::cerr<<"\n";
  std::cerr<<"usage:";
  std::cerr<<"\n";

  for (int co = 0; co < lastOption; co++) {
    if (o[co].key.length() > 0) {
      o[co].key = "-"+o[co].key;
      if (o[co].values.length() > 0)
        o[co].values = " "+o[co].values;
    }
    if (o[co].defaultValue.length() > 0)
      o[co].defaultValue = " (default: "+o[co].defaultValue+")";
    if (o[co].desc.length() > 0)
      o[co].desc = ": "+o[co].desc;
    if (o[co].optionType == OT_required)
      std::cerr<<"\t"<<"<"<<o[co].key.c_str()<<o[co].values.c_str()<<">"<<o[co].defaultValue<<o[co].desc<<"\n";
    else if (o[co].optionType == OT_optional)
      std::cerr<<"\t"<<"["<<o[co].key.c_str()<<o[co].values.c_str()<<"]"<<o[co].defaultValue<<o[co].desc<<"\n";
    else
      throw kaBaseException("unknown optionType!\n");
  }
  std::cerr<<"\t"<<"[-version]\n";
}  // kaVersion::printHelp

void kaVersion::addOption(std::string key, std::string values, optionTypeEnum OT, std::string defaultValue,
                          std::string desc) {
  // std::cerr<<"adding option '"<<key<<"' ...\n";
  if (lastOption == MAXOPTIONS)
    throw kaBaseException("cannot add more than %i options!", MAXOPTIONS);
  o[lastOption].key          = key;
  o[lastOption].values       = values;
  o[lastOption].optionType   = OT;
  o[lastOption].defaultValue = defaultValue;
  o[lastOption].desc         = desc;
  lastOption++;
}
