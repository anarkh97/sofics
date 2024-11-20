#ifndef _CONTROLINFO_H_
#define _CONTROLINFO_H_

#include <cstdio>
#include <cstring>

struct ControlInfo {
   FILE *checkfileptr;
   char *checkfile;
   char *nodeSetName;
   char *elemSetName;
   char *bcondSetName;
   char *decomposition;
   FILE *decPtr;
   char *currentRestartFile;
   char *lastRestartFile;
   char *outputExt;
   char *FlagRST;
   ControlInfo() { checkfile   = (char *) "check"; nodeSetName   = (char *) "nodeset";
                   elemSetName = (char *) "elemset"; bcondSetName = (char *) "bcondset";
                   decomposition = (char *) "DECOMPOSITION";
                   decPtr = 0; checkfileptr = 0;
                   currentRestartFile = 0;
                   lastRestartFile    = 0; outputExt = (char *) ""; FlagRST = (char *) "old"; }
  bool hotRestart() {
    return (lastRestartFile && strcmp(FlagRST,"new") != 0);
  }
};

#endif
