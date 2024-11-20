#include <Hetero.d/FilteredFile.h>
#include <cstdlib>
#include <string>
#include <strings.h>
#include <cstring>

using std::strlen;
using std::strstr;

FilteredFile::FilteredFile(const char *filename)
{
  file = fopen(filename,"r");
  if (file == 0) {
     fprintf(stderr,"File %s not found, giving up\n",filename); exit(1); 
  }
}

int
FilteredFile::findToken(const char *token)
{
  char*     theLine;

  while ((theLine = getLine()) != 0 && (theLine = strstr(theLine, token)) == 0)
 ;

  if(theLine == 0) { return -1; }
  int tlen = strlen(token);
  epos = theLine + tlen;
  return 1;
}

char *
FilteredFile::getLineAfterToken()
{
 return epos;
}

char *
FilteredFile::getLine()
{
 char *cl;
 do {
   cl = fgets(buffer, 200, file);
 }
 while(cl != 0 && (cl[0] == '*' || cl[0] == 0));
 return cl;
}
