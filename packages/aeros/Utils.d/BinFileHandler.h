#ifndef _BIN_FILE_HANDLER_H_
#define _BIN_FILE_HANDLER_H_

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>

#define MAXLINE 500

//------------------------------------------------------------------------------

class BinFileHandler {
public:
#if defined(__SGI) || defined(__SUNPRO_CC) || defined(SALINAS) || defined(WINDOWS)
  typedef long OffType; 
#elif defined(__APPLE__) && defined(__MACH__)
  typedef long OffType;
#elif defined(__LINUX)
  typedef long long OffType;
#else
#error Update the definition of OffType for your machine
#endif

private:

  OffType cpos;

  int headersize;

  double version;

  bool swapBytes;

  int fileid;

  FILE *file;

public:

  BinFileHandler(const char *, const char *, double = 0.0);
  ~BinFileHandler();

  template<class Scalar>
  void swapVector(Scalar *, int);

  template<class Scalar>
  void read(Scalar *, int);

  template<class Scalar>
  void write(Scalar *, int);

  void seek(OffType);

  OffType tell();

  double getVersion() const { return version; }

  int get_fileid() { return fileid; }

};

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::read(Scalar *p, int nobjs)
{
	size_t n1; ssize_t n2;
	if (file)
		n1 = fread(p, sizeof(Scalar), nobjs, file);
	else
		n2 = ::read(fileid, p, nobjs*sizeof(Scalar));

	cpos += nobjs*sizeof(Scalar);

	if (swapBytes) swapVector(p, nobjs);

}

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::write(Scalar *p, int nobjs)
{

  if (swapBytes) swapVector(p, nobjs);

  ssize_t n2;
  if (file) fwrite(p, sizeof(Scalar), nobjs, file);
  else n2 = ::write(fileid, p, nobjs*sizeof(Scalar)); 

  cpos += nobjs*sizeof(Scalar);

  if (swapBytes) swapVector(p, nobjs);

}

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::swapVector(Scalar *p, int nobjs)
{

  for (int obj = 0; obj < nobjs; ++obj) {

    Scalar x = p[obj];

    char *px = (char *) &x;
    char *pp = (char *) (p+obj);

    for (int c = 0; c < int(sizeof(Scalar)); ++c)
      pp[sizeof(Scalar)-1-c] = px[c];

  }

}


//------------------------------------------------------------------------------

inline
void BinFileHandler::seek(BinFileHandler::OffType size) 
{ 

  size += headersize;

  if (file) fseek(file, size, SEEK_SET);
  else lseek(fileid, size, SEEK_SET);

  cpos = size;

}


//------------------------------------------------------------------------------

inline
BinFileHandler::OffType BinFileHandler::tell() 
{ 

  OffType pos;

  if (file) 
    pos = ftell(file);
#ifdef __SGI
  else 
    pos = ::tell(fileid);
#else
  else
    pos = cpos;
#endif

  pos -= headersize;

  return pos;

}

//------------------------------------------------------------------------------

inline
BinFileHandler::BinFileHandler(const char *name, const char *flag, double ver) :
  version(ver), swapBytes(0),  fileid(0), file(0)
{

  int ierr = 0;

  if (std::strcmp(flag, "r") == 0) {
    fileid = open(name, O_RDONLY, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (std::strcmp(flag, "w") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (std::strcmp(flag, "ws") == 0) {
#ifdef WINDOWS
    fileid = open(name, O_WRONLY | O_CREAT | O_TRUNC , 0644);
#else
    fileid = open(name, O_WRONLY | O_CREAT | O_TRUNC | O_SYNC, 0644);
#endif
    if (fileid == -1) ierr = 1;
  }
  else if (std::strcmp(flag, "w+") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (std::strcmp(flag, "ws+") == 0) {
#ifdef WINDOWS
    fileid = open(name, O_WRONLY | O_CREAT, 0644);
#else
    fileid = open(name, O_WRONLY | O_CREAT | O_SYNC, 0644);
#endif 
    if (fileid == -1) ierr = 1;
  }
  else if (std::strcmp(flag, "rb") == 0) {
    file = fopen(name, "rb");
    if (!file) ierr = 1;
  }
  else if (std::strcmp(flag, "wb") == 0) {
    file = fopen(name, "wb");
    if (!file) ierr = 1;
  }
  else if (std::strcmp(flag, "rb+") == 0) {
    file = fopen(name, "rb+");
    if (!file) ierr = 1;
  }
  else {
    fprintf(stderr, "*** Error: wrong flag (%s) for \'%s\'\n", flag, name);
    exit(1);
  }

  if (ierr) {
    fprintf(stderr, "*** Error: unable to open \'%s\'\n", name);
    exit(1);
  }

  headersize = sizeof(int) + sizeof(double);
    
  int one = 1;
  if (std::strcmp(flag, "r") == 0 || std::strcmp(flag, "rb") == 0 || std::strcmp(flag, "rb+") == 0) {
    read(&one, 1);
    if (one != 1) swapBytes = 1;
    read(&version, 1);
  } 
  else if (flag[0] == 'w') {
    write(&one, 1);
    write(&version, 1);
  }

  cpos = headersize;

}

//------------------------------------------------------------------------------

inline
BinFileHandler::~BinFileHandler() 
{ 

  if (file) fclose(file); 
  else close(fileid);
}

//------------------------------------------------------------------------------

inline
int computeNumberOfDigits(int num)
{

  int digits = 1;

  while (num >= 10) {
    num /= 10;
    ++digits;
  }

  return digits;

}

//------------------------------------------------------------------------------

inline
char *computeClusterSuffix(int num, int maxNum)
{

  int numZeros = computeNumberOfDigits(maxNum) - computeNumberOfDigits(num);

  char zeros[100];
  char *suffix = new char[111];

  std::strcpy(zeros, "");
  for (int k=0; k<numZeros; ++k)
    std::strcat(zeros, "0");

  sprintf(suffix, "%s%d", zeros, num);

  return suffix;

}

//------------------------------------------------------------------------------

#endif
