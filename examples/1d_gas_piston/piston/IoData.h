/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <Utils.h>

using std::map;

/*********************************************************************
 * class IoData reads and processes the input data provided by the user
 *********************************************************************
*/
//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:
  map<int, DataType *> dataMap;

  void setup(const char *name, ClassAssigner *p) {
    SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);
    if (p) p->addSmb(name, smo);
    else addSysSymbol(name, smo);
  }

  ~ObjectMap() {
    for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      delete it->second;
  }
};

//------------------------------------------------------------------------------

struct InitialConditions {

  double disp;
  double velo;

  InitialConditions();
  ~InitialConditions() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SpringData {

  double stiffness;
  double mass;

  InitialConditions init_conds;

  SpringData();
  ~SpringData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TimeData {

  double dt;
  double tmax;

  TimeData();
  ~TimeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ProbeData {

  int node_id;
  const char* disp_file;
  const char* velo_file;
  const char* forc_file;

  int frequency;
  double frequency_dt;

  ProbeData();
  ~ProbeData() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct OutputData {

  ObjectMap<ProbeData> probes;

  enum VerbosityLevel {LOW = 0, MEDIUM = 1, HIGH = 2} verbose;

  int frequency;
  double frequency_dt;

  const char* prefix;
  const char* disp_file;
  const char* surf_file;

  OutputData();
  ~OutputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  SpringData spring_data;
  TimeData time_data;
  OutputData output_data;

public:

  IoData() {}
  IoData(int, char**);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void finalize();

};
#endif
