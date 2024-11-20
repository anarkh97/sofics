#ifndef _FILTEREDFILE_H_
#define _FILTEREDFILE_H_

#include <cstdio>

class FilteredFile {
    FILE *file;
    char buffer[200];
    char *epos;
  public:
    FilteredFile(const char *filename);
    int findToken(const char *token);
    char *getLine();
    char *getLineAfterToken();
};

#endif
