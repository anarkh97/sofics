/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Utils.h>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cmath>
#include <unistd.h>
#include <bits/stdc++.h> //INT_MAX
//#include <dlfcn.h>
using namespace std;

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------

InitialConditions::InitialConditions()
{

  disp = 0.0;
  velo = 0.0;

}

//------------------------------------------------------------------------------

void InitialConditions::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassDouble<InitialConditions>(ca, "Displacement", this, &InitialConditions::disp);
  new ClassDouble<InitialConditions>(ca, "Velocity", this, &InitialConditions::velo);

}

//------------------------------------------------------------------------------

SpringData::SpringData()
{

  stiffness = 0.0;
  mass = 1e-3;
  
}

//------------------------------------------------------------------------------

void SpringData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassDouble<SpringData>(ca, "Stiffness", this, &SpringData::stiffness);
  new ClassDouble<SpringData>(ca, "Mass", this, &SpringData::mass);

  init_conds.setup("InitialConditions", ca);

}

//------------------------------------------------------------------------------

TimeData::TimeData()
{

  dt = 0.0;
  tmax = 0.0;

}

//------------------------------------------------------------------------------

void TimeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassDouble<TimeData>(ca, "TimeStepSize", this, &TimeData::dt);
  new ClassDouble<TimeData>(ca, "MaxTime", this, &TimeData::tmax);

}

//------------------------------------------------------------------------------

ProbeData::ProbeData()
{

  node_id = -1;
  disp_file = "";
  velo_file = "";
  forc_file = "";

  frequency = -1;
  frequency_dt = -1;

}

//------------------------------------------------------------------------------

Assigner *ProbeData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 6, nullAssigner);

  new ClassInt<ProbeData>(ca, "NodeId", this, &ProbeData::node_id);
  new ClassStr<ProbeData>(ca, "Displacement", this, &ProbeData::disp_file);
  new ClassStr<ProbeData>(ca, "Force", this, &ProbeData::forc_file);
  new ClassStr<ProbeData>(ca, "Velocity", this, &ProbeData::velo_file);
  new ClassInt<ProbeData>(ca, "Frequency", this, &ProbeData::frequency);
  new ClassDouble<ProbeData>(ca, "TimeInterval", this, &ProbeData::frequency_dt);

  return ca;

}

//------------------------------------------------------------------------------

OutputData::OutputData() {

  verbose = LOW;

  prefix = "";
  disp_file = "";
  surf_file = "";

  frequency = -1;
  frequency_dt = -1;

}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassToken<OutputData>(ca, "VerboseScreenOutput", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::verbose), 3,
                             "Low", 0, "Medium", 1, "High", 2);
  
  new ClassStr<OutputData>(ca, "Prefix", this, &OutputData::prefix);
  new ClassStr<OutputData>(ca, "Displacement", this, &OutputData::disp_file);
  new ClassStr<OutputData>(ca, "Mesh", this, &OutputData::surf_file);
  new ClassInt<OutputData>(ca, "Frequency", this, &OutputData::frequency);
  new ClassDouble<OutputData>(ca, "TimeInterval", this, &OutputData::frequency_dt);

  probes.setup("Probe", ca);

}

//------------------------------------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  //Should NOT call functions in Utils (e.g., print(), exit_mpi()) because the
  //M2C communicator may have not been properly set up.
  readCmdLine(argc, argv);
  readCmdFile();
}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stdout,"\033[0;31m*** Error: Input file not provided!\n\033[0m");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stdout,"\033[0;31m*** Error: could not open \'%s\'\n\033[0m", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stdout,"\033[0;31m*** Error: command file contained parsing errors.\n\033[0m");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//------------------------------------------------------------------------------
// This function is supposed to be called after creating M2C communicator. So, 
// functions in Utils can be used.
void IoData::finalize()
{

}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{

  spring_data.setup("Spring");
  time_data.setup("Time");
  output_data.setup("Output");

}

//------------------------------------------------------------------------------
