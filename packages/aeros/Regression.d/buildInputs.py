#!/usr/bin/python
"""Tool to build input files used for regression testing of AERO codes"""

__author__ = "Mark A. Potts (mpotts@hpti.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2011/02/14"

import sys, os, re, glob, subprocess
#import argparse


def checkFilename(filename,MyList):
# for line in open("/lustre/home/mpotts/FEM/CMakeCache.txt.test"):
  for line in open("../../CMakeCache.txt"):
    if(("MUMPS_common" in line) & ("NOTFOUND" in line) & ("mumps" in filename)):
      return(0)
    if(("SPOOLES_spooles_LIBRARY" in line) & ("NOTFOUND" in line) & ("spooles" in filename)):
      return(0)
    if(("ARPACK_arpack_LIBRARY" in line) & ("NOTFOUND" in line) ):
      for list in MyList:
        for item in list:
          if("arpack" in item):
             return(0)
  return(1) 

def filterInputs(MyList):
  for i in range(len(MyList)):
    for j in range(len(MyList[i])):
      if(("arpack" in MyList[i][j])&(not("*arpack" in MyList[i][j]))):
        item = MyList[i][j]
        stext = "arpack"
        rtext = "*arpack"
        item = item.replace(stext,rtext)
        MyList[i][j] = item
        break

def buildInputs(params):
  mycwd = os.getcwd()
  print "at start, mycwd is %s \n"% mycwd
  if(mycwd.find("Regression.d") < 0):
    if(os.path.exists("Regression.d")==0):
      os.mkdir("Regression.d")  
    os.chdir("Regression.d")  

  ploc = -1
  bpath = re.compile("^\-c")
  i = 0
  copyInputs = 0

  for s in params:
    if(re.search(bpath,s)):
       ploc = i
       copyInputs = 1
    i=i+1

  if(ploc != -1):
    basepath = params[ploc+1]
    del params[ploc+1]
    del params[ploc]

  if(copyInputs == 1): # just copy the all the subdirectories and input files from the basepath parameter
    import shutil
    print "listing for " + basepath
    listdir = os.listdir(basepath)
    mycwd = os.getcwd()
    for infile in listdir:
      fullpath = os.path.join(basepath,infile)
      if (os.path.isdir(fullpath)==True)&(infile[0] != "."):
        if(os.path.exists(os.path.join(mycwd,infile)) != True):
          print "creating directory : " + infile
          newdir = os.mkdir(infile)
        
        print fullpath, os.path.exists(fullpath) 
        origfiles = fullpath + "/*.inp"
        print origfiles, infile
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/run.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/scp.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
  else:
    if(params[1] == 'ALL'):
      PROBLEM_NAMES=['statics','nlstatics','eigen','freqsweep','dynamics','nldynamics',\
                     'impe','tempstatics','tempdynamics','tempnldynamics',\
                     'tempnlstatics','dsvm1','dsvm11','dsvm31','dsvm13',\
                     'dsvm2','dsvm15','dsvm19','dsvm20','dsvm21','dsvm22',\
                     'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm29','dsvm30',\
                     'dsvm31','dsvm32','dsvm34','dsvm35a','dsvm35b','dsvm37','dsvm38',\
                     'dsvm39','dsvm40','vmmech003','vmmech063','vme1','vme2',\
                     'vme3','vme4','vme5','vme6','PreStressedMembrane','PlateUnderPressure']
    else:
      PROBLEM_NAMES = [params[1]]

    p = subprocess.Popen(["hostname"],stdout=subprocess.PIPE)
    retval = p.stdout.readline()
    if(retval.find("independence") != -1):
      host = "independence"
    elif(retval.find("su-ahpcrc") != -1):
      host = "su-ahpcrc"
    else:
      host = "other"
    print "host is %s " % retval
    for problem_type in PROBLEM_NAMES:
  
      if(os.path.exists(problem_type)==0):
        os.mkdir(problem_type)  

      os.chdir(problem_type)
      dirname = os.getcwd()
      command = "rm -f *"
      os.system(command)
      runfilename = "run."+problem_type 
      qsubfilename = "scp."+problem_type 
      RUNFILE = open(runfilename,"w")
      MPIFILE = open(qsubfilename,"w")
      if(host == "su-ahpcrc"):
        MPIFILE.write("#!/bin/bash\n#PBS -N test\n#PBS -V\n#PBS -l nodes=4:ppn=8,walltime=3:00:00\n\n")
      elif(host == "independence"):
        MPIFILE.write("#!/bin/bash\n#PBS -N test\n#PBS -V\n#PBS -l nodes=2:ppn=12,walltime=3:00:00\n\n")

      MPIFILE.write(". /opt/modules/Modules/3.2.6/init/bash\n module load intel openmpi\n")

      MPIFILE.write("cd %s\n" % dirname)
      MPIFILE.write("../create_mfiles.pl\n" )
     
      command = "chmod +x " + runfilename
      os.system(command)
      command = "chmod +x " + qsubfilename
      os.system(command)
#     command = "cp ../*.include ."
#     os.system(command)

      OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy", \
                "strainxz","inxforce","inyforce","inzforce","axmoment", \
                "aymoment","azmoment","energies","gvelocit","gacceler", \
                "displmod","rotatmod","gdispmod"]

      OUTPUT_EXTRAS = ["1"]
      OUTPUT2 = ""
      STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr","FETI",\
                 "FETI DP"]

      INCLUDE = ["\"../mesh.include\""]

      DYNAMICS = ["time\t0.0\t3.0e+0\t3.0e+0",\
                  "time\t0.0\t3.0e+1\t3.0e+1",\
                  "time\t0.0\t3.0e+2\t3.0e+2"]

      IMPE =      ["freq 10.0\nraydamp 1e-6 1.0",\
                   "freq 10.0\nraydamp 1e-7 1.0",\
                   "freq 10.0\nraydamp 1e-5 1.0"]

      NONLINEAR = ["maxitr 100\nnltol 1.0e-6\nrebuild 1",\
                   "maxitr 100\nnltol 1.0e-7\nrebuild 1",\
                   "maxitr 100\nnltol 1.0e-5\nrebuild 1"]

      EIGEN = ["arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-05",\
               "arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-04",\
               "arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-6"]

      SHIFT = ["0","0.1"]

      if(problem_type == "PlateUnderPressure"):
        OUTPUT = ["displacx","gvelocit"]
        OUTPUT2 = ["displacx"]
        OUTPUT_EXTRAS = [" 200 65"," 200 65"]
        DYNAMICS = [\
          "mech\t0.0\t0.5\ntime 0.0  6.4549E-06 3.0\nstable 0\nNONLINEAR\nMATLAW\n1 J2Plasticity 2.100E+11 3.000E-01 7.850E+03 2.500E+08 0.0 0.0 1.0e-10",\
          "mech\t0.0\t0.5\ntime 0.0  5.8234E-06 3.0\nstable 0\nNONLINEAR\nMATLAW\n1 HypoElastic 2.100E+11 3.000E-01 7.850E+03"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]

      if(problem_type == "vme6"):
        OUTPUT = ["displacx"]
        OUTPUT2 = ["displacy"]
        OUTPUT_EXTRAS = [" 10 1"," 10 1"]
        DYNAMICS = ["mech\t0.0\t0.5\ntime 0.0 0.0001 2.0\nraydamp 0 1.0"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\n*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]

      if(problem_type == "vme5"):
        OUTPUT = ["rotatioz"]
        OUTPUT_EXTRAS = [" 1 7"]
        DYNAMICS = ["mech\t0.0\t0.5\ntime 0.0 1.0e-6 0.4"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\n*","NONLINEAR\n*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]

      if(problem_type == "vme4"):
        OUTPUT = ["displacx"]
        OUTPUT2 = ["displacx"]
        OUTPUT_EXTRAS = [" 1 2"," 1 3"]
        DYNAMICS = ["mech\t0.0\t0.5\ntime 0.0 1.0e-3 2.0"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\n*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]

      if(problem_type == "vme3"):
        OUTPUT = ["displacx","gvelocit"]
        OUTPUT_EXTRAS = [" 1 2"]
        DYNAMICS = ["mech\t0.0\t0.5\ntime 0.0 0.001 2.0\nraydamp 0.0342905309 0"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\n*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INC_FILE_2 = INCLUDE_FILE + "\nMFTT\n0.0 0.0\n1.0 1.0\n2.0 2.0"
        INCLUDE = [INCLUDE_FILE,INC_FILE_2]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]

      if(problem_type == "vme2"):
        OUTPUT = ["displacy"]
        OUTPUT2 = ["gdisplac","normal_force_mag","normal_traction_mag","cdirnorx","cdirnory","cdirnorz","contact_area"]
        OUTPUT_EXTRAS = [" 500 94"," 50000"," 50000"," 50000"," 50000"," 50000"," 50000"," 50000"]
        DYNAMICS = ["time\t0.0\t2e-7\t1.0"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\nmech 0.0 0.5\nstable 0\n*","NONLINEAR\n*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]
 
      if(problem_type == "vme1"):
        OUTPUT = ["displacx"]
        OUTPUT_EXTRAS = [" 10 2"]
        DYNAMICS = ["time\t0.0\t0.001\t5.0"]
        NAMELIST = ["DYNAMICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["newmark\nmech 0.0 0.5 0.0 0.0\n*","*","*"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [DYNAMICS,OUTPUT,INCLUDE]
 
      if(problem_type == "PreStressedMembrane"):
        OUTPUT = ["gdisplac"]
        OUTPUT2 = ["stressxx","stressxy","stressyy"]
        OUTPUT_EXTRAS = [" 1"," 1 elemental"," 1 elemental"," 1 elemental"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps","sparse","spooles"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","NONLINEAR\nrebuild 1\nnltol 1e-6\nmaxit 10","*"]

      if(problem_type == "vmmech063"):
        OUTPUT = ["displacx"]
        OUTPUT2 = ["displacz","stressxx"]
        OUTPUT_EXTRAS = [" 1 1326"," 1 1326"," 1 11 upper"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","NONLINEAR\nmaxitr 10\nnltol 1.0e-6\nrebuild 1\ndlambda 0.08333333333 1.0\nunsymmetric","*"]

      if(problem_type == "vmmech003"):
        OUTPUT = ["geigenpa"]
        OUTPUT_EXTRAS = [" 10 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps", "mumps pivot","FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*"]

      if(problem_type == "dsvm40"):
        OUTPUT = ["displacy"]
        OUTPUT_EXTRAS = [" 1 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps","skyline","gmres",\
                   "mumps pivot","pcg","bcg","cr",\
                   "FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"\ncases 1 2","*","*","*"]

      if(problem_type == "dsvm39"):
        OUTPUT = ["displacx"]
        OUTPUT2 = ["displacx"]
        OUTPUT_EXTRAS = [" 10 371"," 10 1100"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        DYNAMICS = ["time\t0.0\t5.0e-3\t3.0"]
        STATICS = ["sparse","mumps","spooles",\
                 "spooles pivot","mumps pivot","FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","newmark\nmech 0.25 0.5 0.0 0.0\n*","*","*"]

      if(problem_type == "dsvm38"):
        OUTPUT = ["gvelocit"]
        OUTPUT2 = ["gvelocit"]
        OUTPUT_EXTRAS = [" 1 18"," 1 586"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        DYNAMICS = ["time\t0.0\t0.5e-5\t1.0"]
        STATICS = ["FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","newmark\nmech 0.8\n*","*","*"]

      if(problem_type == "dsvm37"):
        OUTPUT = ["gtempera"]
        OUTPUT_EXTRAS = ["1 42619"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        DYNAMICS = ["time\t21.6\t0\t21600"]
        STATICS = ["sparse","mumps","spooles",\
                 "spooles pivot","mumps pivot", "FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","newmark\nheat 0.5\n*","*","*"]

      if(problem_type == "dsvm35b"):
        OUTPUT = ["strainp1"]
        OUTPUT2 = ["displacz","stresszz"]
        OUTPUT_EXTRAS = [" 1"," 1"," 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps","spooles","gmres",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr",\
                 "FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*","*"]

      if(problem_type == "dsvm35a"):
        OUTPUT = ["gtempera"]
        OUTPUT2 = ["gtempera"]
        OUTPUT_EXTRAS = ["1 NG 1","1 NG 4"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps pivot"]
        NONLINEAR = ["maxitr 100\nnltol 1.0e-10\nrebuild 1"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"\n*\nTRBM\n1e-8\n*","*","*","*"]

      if(problem_type == "dsvm34"):
        OUTPUT = ["displmod"]
        OUTPUT_EXTRAS = [" 20"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse"]
        NONLINEAR = ["maxitr 20\nnltol 1.0e-6\ndlambda 0.05\t 1.0\npenalty 10 1.0e-6 1"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
        EXTRAS = ["CONSTRAINTS\naugmented 1e9\n*","*","*","*"]

      if(problem_type == "dsvm32"):
        OUTPUT = ["gtempera"]
        OUTPUT2 = ["gtempera"]
        OUTPUT_EXTRAS = ["1 NG 2","1 NG 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*","*"]

      if(problem_type == "dsvm30"):
        OUTPUT = ["displmod"]
        OUTPUT2 = ["stressxx","reaction"]
        OUTPUT_EXTRAS = [" 1","1 NG 2","1 NG 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm29"):
        OUTPUT = ["stressvm"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps pivot","spooles pivot"]
        NONLINEAR = ["maxitr 60\nnltol 1.0e-10\ndlambda 0.25 1.5"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm27b"):
        OUTPUT = ["displacz"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","skyline",\
                   "spooles","pcg","bcg","cr"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm27a"):
        OUTPUT = ["gtempera"]
        OUTPUT2 = ["heatflxz","heatreac"]
        OUTPUT_EXTRAS = [" 1"," 1"," 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","skyline",\
                   "spooles","pcg","bcg","cr"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*"]

      if(problem_type == "dsvm25"):
        OUTPUT = ["stressp1"]
        OUTPUT2 = ["stressp1"]
        OUTPUT_EXTRAS = ["1 NG 1","1 NG 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["tolfeti 1.0e-8\naux_coarse_solver solverhandle 1\nconstraints multipliers\n*\nSOLVERCNTL 1\nskyline\ntrbm 1e-5\n*\nGRBM\n*","*","*","*"]

      if(problem_type == "dsvm24"):
        OUTPUT = ["displacz"]
        OUTPUT_EXTRAS = [" 1 NG 1","1 NG 1 modphase"]
        IMPE = ["freq 500","freq 500","freq 500","freq 500\nraydamp 3.18309886e-5 0","freq 500\nraydamp 3.18309886e-5 0","freq 500\nraydamp 3.18309886e-5 0"]
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["spooles pivot","mumps pivot","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","include \"../fetidph.include\"","*","*"]

      if(problem_type == "dsvm23"):
        OUTPUT = ["displacy"]
        OUTPUT_EXTRAS = [" 1 11"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","pcg","bcg","cr",\
                   "FETI DP"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*"]

      if(problem_type == "dsvm22"):
        OUTPUT = ["gdisplac","displacx","stressvm","strainvm",\
                  "displacx","stressvm","strainvm","gdisplac"]
        OUTPUT2 = ["displacz","displacz","stresszz","stresszz"]
        OUTPUT_EXTRAS = [" 1"," 1 NG 1"," 1 NG 2"," 1 NG 1"," 1 NG 2"]
        NAMELIST = ["STATICS\n","","OUTPUT\n","INCLUDE ","INCLUDE "]
        STATICS = ["sparse","FETI DP"]
        STATICS_OPTS = ["constraints penalty 1e12","constraints multipliers"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        CONTACT = ["../dsvm22.contact1","../dsvm22.contact1","../dsvm22.contact2","../dsvm22.contact2",\
                   "../dsvm22.contact3","../dsvm22.contact3","../dsvm22.contact4","../dsvm22.contact4"]
        OPTIONSLIST = [STATICS,STATICS_OPTS,OUTPUT,INCLUDE,CONTACT]
#        EXTRAS = ["constraints penalty 1e12","*","*","*"]
        EXTRAS = ["*","*","*","*","*"]


      if(problem_type == "dsvm21"):
        OUTPUT = ["displacy"]
        NAMELIST = ["STATICS\n","EIGEN\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","pcg","bcg","cr"]
        EIGEN = ["arpack\nnsbspv 3\nneigpa 1","arpack\nnsbspv 3\nneigpa 1"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,EIGEN,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*\nGEPS\nBUCKLE\n*","*"]

      if(problem_type == "dsvm20"):
        OUTPUT = ["gdispmod"]
        NAMELIST = ["STATICS\n","EIGEN\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP"]
        EIGEN = ["arpack\nnsbspv 8\nneigpa 4"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,EIGEN,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*","*"]

      if(problem_type == "dsvm19"):
        OUTPUT = ["displacy"]
        OUTPUT_EXTRAS = ["1 27"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP"]
        INCLUDE = ["../dsvm19.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*"]

      if(problem_type == "dsvm15"):
        OUTPUT = ["displacy"]
        OUTPUT2 = ["displacy"]
        OUTPUT_EXTRAS = ["1 38046","1 37735"]
        IMPE =    ["freqsweep 0 500 11 50\nrecons pade 2 4 5"]
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["spooles pivot","mumps pivot","FETI DPH"]
        INCLUDE = ["../dsvm15.include"]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","include \"../fetidph.include\"","*","*"]

      if(problem_type == "dsvm2"):
        OUTPUT = ["stressxx","strainxx"]
        OUTPUT_EXTRAS = ["1 NG 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps pivot","mumps","spooles",\
                   "spooles pivot", "FETI DP"]
        INCLUDE = ["../dsvm2.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"","*","*"]

      if(problem_type == "dsvm13"):
        OUTPUT = ["geigenpa"]
        NAMELIST = ["STATICS\n","EIGEN\n","GEPS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps pivot", "spooles pivot"]
        GEPS = ["buckle"]
        EIGEN = ["arpack LA 4\nnsbspv 3\nneigpa 1\ntoleig 1.0e-6\nshift 100"]
        INCLUDE = ["../dsvm13.include"]
        OPTIONSLIST = [STATICS,EIGEN,GEPS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*","*"]

      if(problem_type == "dsvm31"):
        OUTPUT = ["stressvm","strainvm","strainxx",\
                "stressp3","stressxx","strainp3"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres",\
                   "spooles pivot","mumps pivot","pcg","bcg","cr",\
                   "FETI DP"]
        INCLUDE = ["../dsvm31.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["include \"../fetidp.include\"\n*","*","*"]

      if(problem_type == "dsvm11"):
        OUTPUT = ["displmod","gdisplac","displacz"]
        OUTPUT_EXTRAS = ["1 5"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        NONLINEAR = ["maxitr 40\nnltol 1.0e-6\nrebuild 1",\
                   "maxitr 10\nnltol 1.0e-5\nrebuild 1"]
        INCLUDE = ["../dsvm11.include"]
        EXTRAS = ["include \"../fetidp.include\"\n*","*","*","*"]
        STATICS = ["sparse","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]

        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]

      if(problem_type == "dsvm1"):
        OUTPUT = ["reaction"]
        OUTPUT2 = ["reaction"]
        OUTPUT_EXTRAS = [" 1 NG 1"," 1 NG 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres",\
                   "spooles pivot","mumps pivot","pcg","bcg","cr"]
        INCLUDE = ["../dsvm1.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "tempnlstatics"):
        STATICS = ["sparse","skyline","mumps","spooles", "spooles pivot","mumps pivot",\
                 "FETI DP"]
        EXTRAS = ["include \"../fetidp.include\"","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0 0.0 202.4 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempnldynamics"):
        DYNAMICS = ["heat\t0.5\ntime\t1.0\t1.0\t200.0"]
        EXTRAS = ["include \"../fetidp.include\"","newmark","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        STATICS = ["sparse","skyline","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]
        OUTPUT = ["gtempera","gtempvel"]
        OUTPUT_EXTRAS = [" 200"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempdynamics"):
        DYNAMICS = ["mech\t0.25\t0.5",\
                    "time\t1.0\t1.0\t200.0"]
        OUTPUT_EXTRAS = [" 200"]
        EXTRAS = ["include \"../fetidp.include\"","newmark","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        STATICS = ["sparse","skyline","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]
        OUTPUT = ["gtempera","gtempvel"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempstatics"):
        NAMELIST = ["STATICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr",\
                 "FETI DP"]
        EXTRAS = ["include \"../fetidp.include\"","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc..."]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0    0.0 202.4 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0"]
        OPTIONSLIST=[STATICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "dynamics"):
        EXTRAS = ["*","newmark\nmech\t0.25000\t0.5000\n*\ttime step\ttotal time","*","*","*"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]
        OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy","strainvm", \
                "strainxz","displmod","gdispmod","displacx","displacy","displacz",\
                "strainp1","strainp2","strainp3","stressp1","stressp2","stressp3",\
                "gvelocit","gacceler"]
        DYNAMICS = ["time\t0.0\t3.0e+0\t3.0e+2"]
        OUTPUT_EXTRAS = [" 100"]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
  
      if(problem_type == "nldynamics"):
        STATICS = ["sparse","spooles","FETI DP","mumps"]
        DYNAMICS = ["time\t0.0\t0.3e+0\t3.0e+1"]
        OUTPUT_EXTRAS = [" 100"]
        OUTPUT = ["gdisplac","stressvm","strainxx","strainxz",\
                  "stressxx","stressxz", "gvelocit","gacceler" ]
        EXTRAS = ["*","newmark\nmech\t0.8000\n*\ttime step\ttotal time","*","*","*","*","*"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,NONLINEAR,OUTPUT,INCLUDE]

      if(problem_type == "impe"):
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["spooles pivot","mumps pivot","FETI DPH"]
        OUTPUT = ["gdisplac","displmod","gdispmod","displacx","displacy","displacz"]
                
        INCLUDE = ["\"../mesh.include\""]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","include \"../fetidph.include\"","*","*","*"]
  
      if(problem_type == "freqsweep"):
        NAMELIST = ["STATICS\n","IMPE\n","OUTPUT\n","INCLUDE "]
        STATICS = ["spooles pivot","mumps pivot","FETI DPH"]
        IMPE =    ["freqsweep 1. 3. 3 10\nraydamp 1e-6 1.0",\
                   "freqsweep 1. 3. 3 10\nraydamp 1e-7 1.0",\
                   "freqsweep 1. 3. 3 10\nraydamp 1e-5 1.0"]
        OUTPUT = ["gdisplac","displmod","gdispmod","displacx","displacy","displacz"]
        OPTIONSLIST = [STATICS,IMPE,OUTPUT,INCLUDE]
        EXTRAS = ["*","recons pade 2 9 10","*","*","*"]

      if(problem_type == "nlstatics"):
        EXTRAS = ["include \"../fetidp.include\"\n*","dlambda 0.1 0.3\n*","*","*","*"]
        STATICS = ["sparse","skyline","mumps","spooles","spooles pivot","mumps pivot",\
                 "FETI DP"]
        OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy","strainvm", \
                "strainxz","displmod","gdispmod","displacx","displacy","displacz",\
                "strainp1","strainp2","strainp3","stressp1","stressp2","stressp3"]
  
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
      if(problem_type == "statics"):
        STATICS = ["sparse","skyline","mumps","spooles","gmres",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr",\
                 "FETI DP"]
        OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy","strainvm", \
                "strainxz","displmod","gdispmod","displacx","displacy","displacz",\
                "strainp1","strainp2","strainp3","stressp1","stressp2","stressp3"]
                
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["*","*","*"]
        OPTIONSLIST=[STATICS,OUTPUT,INCLUDE]

      if(problem_type == "eigen"):
        STATICS = ["spooles pivot","mumps pivot","FETI DPH"]
        EXTRAS = ["include \"../fetidph.include\"","*","*","*","*","*"]
        NAMELIST = ["STATICS\n","EIGEN\n","SHIFT ","OUTPUT\n","INCLUDE "]
        INCLUDE = ["\"../mesh_eigen.include\""]
        OUTPUT = ["geigenpa"]
        OPTIONSLIST=[STATICS,EIGEN,SHIFT,OUTPUT,INCLUDE]
 
 
      options = 0
      for i in range(len(OPTIONSLIST)):
        if(len(OPTIONSLIST[i]) >  options):
          options = len(OPTIONSLIST[i])
      for i in range(options):
        idname = problem_type
        for j in range(len(OPTIONSLIST)-1):
          if((NAMELIST[j] != "")&(OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find(" ") != 1)&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("-") == -1 )&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("\n") == -1 )&
                  (NAMELIST[j].find("INCLUDE") == -1 )&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("\t") == -1 )):
            idname = idname + "_" + OPTIONSLIST[j][i % len(OPTIONSLIST[j])]

        idname = idname.replace(" ","_",10)
  

#       if(problem_type == "dsvm1"): 
#         OUTPUT_FILENAME = idname+"_1.dat 1" + " NG 1\n" 
#         OUTPUT_FILENAME = OUTPUT_FILENAME + "reaction "+ idname+"_2.dat 1" + " NG 2" 
        OUTPUT_FILENAME = idname+".dat" 

#       filterInputs(OPTIONSLIST) 
        filename = idname+".inp"
        if(checkFilename(filename,OPTIONSLIST)): 
          FILE = open(filename,"w")
          FILE.write("CONTROL\n");
          FILE.write(idname);
          if(idname.find("temp") != -1):
            FILE.write("\n2\n\"nodeset\"\n\"elemset\"\n*\n");
          else:
            FILE.write("\n1\n\"nodeset\"\n\"elemset\"\n*\n");
          for j in range(len(OPTIONSLIST)):
            FILE.write(NAMELIST[j])
            FILE.write(OPTIONSLIST[j][i  % len(OPTIONSLIST[j])])
            if(NAMELIST[j].find("OUTPUT") != -1 ):
              # create an empty output file so we can then use its existence to activate the comparison even if the
              # code crashes before creating the file
              command = "touch " + OUTPUT_FILENAME
              os.system(command)
              FILE.write(" %s" % OUTPUT_FILENAME)
              FILE.write(" %s" % OUTPUT_EXTRAS[0])
              if(OUTPUT2 != ""):
                for jj in range(len(OUTPUT2)):
                  FILE.write("\n")
                  OUTPUT_FILENAME = idname + "_" + "%d.dat" % ( jj + 2)
                  # create an empty output file so we can then use its existence to activate the comparison
                  # code crashes before creating the file
                  command = "touch " + OUTPUT_FILENAME
                  os.system(command)
                  FILE.write("%s %s" % (OUTPUT2[jj],OUTPUT_FILENAME))
                  FILE.write(" %s" % OUTPUT_EXTRAS[jj+1])
            FILE.write("\n")
            if(NAMELIST[j].find("FETI DPH") != -1):
              FILE.write("include \"../fetidph.include\"\n*")
            elif(NAMELIST[j].find("FETI DP") != -1):
              FILE.write("include \"../fetidp.include\"\n*")
            else:
              FILE.write(EXTRAS[j])
            FILE.write("\n")
          FILE.write("END\n")
          FILE.close()
          command = "../../bin/aeros "
          re.compile("FETI");
          if(idname.find("FETI") != -1 ):
            command = command + "-n 2 --dec --nsub 4"
          print "Creating %s" % filename
          MPIFILE.write("echo mpirun -n 2 %s %s\n" % (command,filename.replace(" ","_")))
          MPIFILE.write("mpirun -n 2 --machinefile host.%d %s %s &\n" % (i,command,filename.replace(" ","_")))
          RUNFILE.write("echo %s %s\n" % (command,filename.replace(" ","_")))
          RUNFILE.write("%s %s\n" % (command,filename.replace(" ","_")))
      MPIFILE.write("wait ")
      for ii in range(options):
        MPIFILE.write("%")
        MPIFILE.write("%d " %(ii+1))
      MPIFILE.write("\n")
      os.chdir('../')

if __name__ == "__main__":

# parser = argparse.ArgumentParser(description='Build Regression Testing Inputs.')
# parser.add_argument('--c')
# basepath = parser.parse_args(sys.argv.split())

  params = sys.argv
  buildInputs(params)
