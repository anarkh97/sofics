#!/usr/bin/python
"""Tool designed to compare data files produced by regression testing and report discrepancies"""

__author__ = "Mark A. Potts (mpotts@drc.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2011/02/17"

import sys, os, re, md5, subprocess, math, glob, datetime, time
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


def directComp(basefile,file,SUMMARY_FILE,outstring):
  COMP = list(open(file).read().splitlines())
  BASE = list(open(basefile,"r").read().splitlines())
  if(file.find("eigen") != -1):  
    absComp = 1
  else:
    absComp = 0
  MaxDiff = 0.0
  RelDiff = 0.0
  TotDiff = 0.0
  nSample = 0
  avgSum = 0.0
  RelDiffVals = (0,0)
  MaxDiffLine = -1
  RelDiffLine = -1
  MaxDiffLoc = -1
  RelDiffLoc = -1
  MaxDiffVals = (-1,-1)
  if(len(BASE) != len(COMP)):
    print bcolors.WARNING + "WARNING-- %s and %s have different lengths!" % (file,basefile) + bcolors.ENDC
    print len(BASE)
    print len(COMP)
    outstring.append( "\tBaseline file length = %d, New file length = %d\n" %(len(BASE),len(COMP)))
    SUMMARY_FILE.write(outstring[0])
    result = 1
    return result
   
  for i in range(len(BASE)):
    if(i > 0):
      BASE[i].replace(","," ");
      COMP[i].replace(","," ");
      BASE[i].strip()
      COMP[i].strip()
      basewords = BASE[i].split()
      compwords = COMP[i].split()
      if(len(basewords) != len(compwords)):
#       print "WARNING--Different line lengths at line %d of %f!!" %(i,file)
        outstring.append( "\tWARNING--Different line lengths at line %d of %f!!" %(i,file) )
        result = 1
      else:
        for j in range(len(basewords)):
          if(not(basewords[j].isalpha())and not(compwords[j].isalpha()) and(basewords[j].find("__") == -1)):
            if(absComp):
              basewords[j] = math.fabs(float(basewords[j]))
              compwords[j] = math.fabs(float(compwords[j]))
            diff = float(basewords[j]) - float(compwords[j])
            avg = 0.5*(float(basewords[j]) + float(compwords[j]))
            TotDiff = math.fabs(diff) + TotDiff
            nSample = nSample + 1
            avgSum = avgSum + avg
            if((0.5*(float(basewords[j]) + float(compwords[j])))>= 1.0e-5):
              if(math.fabs(diff/avg) > RelDiff):
                RelDiff = math.fabs(diff/avg)
                RelDiffLine = i
                RelDiffLoc = j
                RelDiffVals = (basewords[j],compwords[j])
            if(math.fabs(diff) >= MaxDiff):
              MaxDiff = math.fabs(diff)
              MaxDiffLine = i
              MaxDiffLoc = j
              MaxDiffVals = (basewords[j],compwords[j])

  if(RelDiff > 0.0):
#   print bcolors.OKBLUE + "\t\t\t\tRel diff is %e on line %d at location %d" %(RelDiff,RelDiffLine,RelDiffLoc) + bcolors.ENDC
    outstring.append( "\tRel diff is %e on line %d at location %d\n" %(RelDiff,RelDiffLine,RelDiffLoc))
  else:
    outstring.append( "")
#  print bcolors.OKBLUE + "\t\t\t\tvals were %s " % (RelDiffVals,) + bcolors.ENDC
# outstring.append( "\tvals were %s \n" % (RelDiffVals,))
# SUMMARY_FILE.write(outstring[0])
  if((MaxDiff > 1.0e-8)&(RelDiff > 5.0e-2)):
#   TotDiff = TotDiff/nSample
#   print "Average Difference = %e %e %e %d \n" % (TotDiff/avgSum,TotDiff/nSample,avgSum,nSample)
#   print "Max diff was %e and RelDiff was %e \n" % (MaxDiff,RelDiff)
    result = 1
  else:
    result = 0
  return result

def dComp(params):
  import time, glob
  import subprocess

  summary_filename = "reg_test_summary"
  mycwd = os.getcwd()
  if(mycwd.find("Regression.d") == -1):
    os.chdir("Regression.d")
  if(os.path.exists(summary_filename)):
    last_mod = os.path.getmtime(summary_filename)
    time_now = time.mktime(time.localtime())
    print "time now = %d and last_mod = %d, diff = %d \n" % (time_now,last_mod,time_now - last_mod)
    if((time_now - last_mod) > 300):
      SUMMARY_FILE = open(summary_filename,"w")
    else:
      SUMMARY_FILE = open(summary_filename,"a")
  else:
    SUMMARY_FILE = open(summary_filename,"w")
  outstring = "\n\nSummary of regression test run as:\n%s\nfrom:%s" % (params,os.getcwd())
  SUMMARY_FILE.write(outstring)
  now = datetime.datetime.now()
  time_now = now.time()
  date_now = now.date()
  outstring = "\n\nTest started at %s on %s\n"% (time_now,date_now)
  SUMMARY_FILE.write(outstring)
  exactMatches = 0
  Total = 0
  loc = -1
  sloc = -1
  rloc = -1
  nloc = -1
  gloc = -1
  i = 0 
  pattern = re.compile("^\-r")
  runLocal = re.compile("^\-l")
  sendMail = re.compile("^\-s")
  newPlots = re.compile("^\-n")
  genBase = re.compile("^\-g")

  for s in params:
    if(pattern.match(s)):
       loc = i
    if(sendMail.match(s)):
       sloc = i
    if(runLocal.match(s)):
       rloc = i
    if(newPlots.match(s)):
       nloc = i
    if(genBase.match(s)):
       gloc = i
    i=i+1

  if(sloc != -1):
    sendMail = 1
    del params[sloc]
    if(nloc > sloc): nloc = nloc -1
    if(rloc > sloc): rloc = rloc -1
    if(loc > sloc): loc = loc -1
  else:
     sendMail = 0
 
  if(loc != -1):
    run = 1
    del params[loc]
    if(nloc > loc): nloc = nloc -1
    if(rloc > loc): rloc = rloc -1
    if(sloc > loc): sloc = sloc -1
  else:
    run = 0 

  if(rloc != -1):
    lrun = 1
    del params[rloc]
    if(nloc > rloc): nloc = nloc -1
    if(sloc > rloc): sloc = sloc -1
    if(loc > rloc): loc = loc -1
  else:
    lrun = 0 

  if(gloc != -1):
    genbase = 1
    lrun = 1
    del params[gloc]
  else:
    genbase = 0;

  if(nloc != -1):
    newP = 1
    del params[nloc]
    if(sloc > nloc): sloc = sloc -1
    if(rloc > nloc): rloc = rloc -1
    if(loc > nloc): loc = loc -1
    os.system("rm gnuplot_create")
  else:
    newP = 0


  batch = 0
  files = [] 
  plotList = []
  if((params[1] == 'ALL')|(params[1] == 'short')):
     
    if(params[1] == 'ALL'):
      PROBLEM_NAMES=['statics','nlstatics','eigen','dynamics','nldynamics',\
                   'freqsweep','impe','tempstatics','tempnlstatics',\
                   'tempdynamics','tempnldynamics','dsvm1','dsvm11','dsvm31',
                   'dsvm2','dsvm13','dsvm15','dsvm19','dsvm20','dsvm21','dsvm22',\
                   'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm29','dsvm30',\
                   'dsvm32','dsvm34','dsvm35a','dsvm35b','dsvm37','dsvm38','dsvm39',\
                   'dsvm40','vme1','vme2','vme3','vme4','vme5','vme6','vmmech003',\
                   'vmmech063','PreStressedMembrane','PlateUnderPressure']
    else:
      PROBLEM_NAMES=['nlstatics','freqsweep','impe','tempstatics','tempnlstatics',\
                   'tempdynamics','tempnldynamics','dsvm1','dsvm31',
                   'dsvm20','dsvm21','dsvm22',\
                   'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm30',\
                   'dsvm35b','dsvm40']
    del params[0:1]

    for names in PROBLEM_NAMES:
      indir = names+"/"
      os.chdir(indir)
      if(run == 1):
        os.system("rm *.dat test.* host.*")
        qsubScript = "scp."+names
        retcode = subprocess.Popen(["qsub", qsubScript ], stdout=subprocess.PIPE).communicate()[0]
        numbers = filter(lambda x: x.isdigit(), retcode)
        qsubRetfile = "test.o"+numbers
        while True:
          if(os.path.exists(qsubRetfile)):
            break
          time.sleep(10)
      if(batch == 1):
        os.system("rm -f *.dat test.* host.*")
        command = "./scp."+names+" >reg.out 2>&1"
        os.system(command)
      if(lrun == 1):
        os.system("rm -f *.dat test.* host.*")
        command = "./run."+names+" >reg.out 2>&1"
        os.system(command)
      os.chdir('../')
      for infile in glob.glob( os.path.join(indir, '*.dat') ):
        files.append(infile)
  else: 
    del params[0]
    indirs = params
    for indir in indirs:
      indirp = indir+"/"
      os.chdir(indir)
      if(run == 1):
        os.system("rm -f *.dat test.* host.*")
        qsubScript = "scp."+indir
        retcode = subprocess.Popen(["qsub", qsubScript ], stdout=subprocess.PIPE).communicate()[0]
        numbers = filter(lambda x: x.isdigit(), retcode)
        qsubRetfile = "test.o"+numbers
        while True:
          if(os.path.exists(qsubRetfile)):
            break
          time.sleep(10)
      if(batch == 1):
        os.system("rm -f *.dat test.* host.*")
        command = "./scp."+indir+" >reg.out 2>&1"
        os.system(command)
      if(lrun == 1):
        os.system("rm -f *.dat test.* host.*")
        command = "./run."+indir+" >reg.out 2>&1"
        os.system(command)
      os.chdir('../')
      for infile in glob.glob( os.path.join(indirp, '*.dat') ):
        files.append(infile)

  result = 0
  if(genbase != 1): 
    for file in files:
      basefile = "baseline/"+file
      file = mycwd+"/Regression.d/"+file
      p = subprocess.Popen("md5sum " + file, shell = True, stdout=subprocess.PIPE) # don't forget to "import subprocess"
      p.wait()
      Total += 1
      if(p.returncode == 0):
        newmd5 = p.stdout.readline().split() # now the md5 variable contains the MD5 sum
        p.wait() # some clean up
        p = subprocess.Popen("md5sum " + basefile, shell = True, stdout=subprocess.PIPE) # don't forget to "import subprocess"
        oldmd5 = p.stdout.readline().split() # now the md5 variable contains the MD5 sum
        p.wait() # some clean up
        if(newmd5[0] == oldmd5[0]):
          print bcolors.OKGREEN + " \tExact match " + bcolors.ENDC, file
          exactMatches += 1
          outstring = "\texact match " + file + "\n"
          SUMMARY_FILE.write(outstring)
        else:
          compstring = []
          dcResult = directComp(basefile,file,SUMMARY_FILE,compstring)
          if(dcResult == 1):
            result = 1
            print bcolors.FAIL + " \tDiscrepancy " + bcolors.ENDC, file
            outstring = "\tDiscrepancy " + file + "\n"
            SUMMARY_FILE.write(outstring)
            SUMMARY_FILE.write(compstring[0])
            gnuplotCalled = 1;
            plotList.append([file,basefile]);
          else:
            exactMatches += 1
            print bcolors.OKGREEN + " \tMatch" + bcolors.ENDC, file
            if(compstring != ""):
              print bcolors.OKBLUE    + "\t" +compstring[0] + bcolors.ENDC 
            outstring = "\tclose match " + file + "\n"
            SUMMARY_FILE.write(outstring)
            SUMMARY_FILE.write(compstring[0])
      else: 
        result = 1
        print bcolors.FAIL + " \tNo file" + bcolors.ENDC, file 
        outstring = "\tno match " + file + "\n"
  if(genbase == 1):
    outstring = "New baseline(s) generated\n"
  else:
    outstring = "%s--%d matches found out of %d total cases\n" %(indir,exactMatches,Total)
  print outstring
  SUMMARY_FILE.write(outstring)
  now = datetime.datetime.now()
  time = now.time()
  date = now.date()
  outstring = "Test completed at %s on %s\n"% (time,date)
  SUMMARY_FILE.write(outstring)

  if(len(plotList) > 0):
    if(os.path.exists("gnuplot_create")):
      PLOT_FILE = open("gnuplot_create","a")
    else:
      PLOT_FILE = open("gnuplot_create","w")
      PLOT_FILE.write("set terminal postscript color\n")
      PLOT_FILE.write("set output \"Discrepancies.ps\"\n")
      PLOT_FILE.write("set lmargin 10\nset rmargin 5\nset tmargin 5\nset bmargin 5\nset key bottom\nset key box\n")

    for files in plotList:
      words = files[0].split('/');
      title1 = words[1]
      title2 = "baseline--%s"%words[1]
      TEST_FILE = open(files[0],"r")
      linenum = 0
      count = 0
      for line in TEST_FILE:
        if (linenum == 2):
          words = line.split()
          count = len(words)
          break
        linenum = linenum+1
      if(title1.find("dsvm13") != -1):
        skip_num = 3
      elif(title1.find("dynamics") != -1):
        p = subprocess.Popen(["head","-2",files[0] ], stdout=subprocess.PIPE) 
        line_nums = p.stdout.readline().split() 
        line_nums = p.stdout.readline().split() 
        if(line_nums[0].find("___") == -1):
          skip_num = int(line_nums[0]) + 4
        else:
          skip_num = 4
      else:
        skip_num = 2
      if(count == 0):
        PLOT_FILE.write("set title \"%s--MISSING\"\n" % (files[0]))
        PLOT_FILE.write("plot \"%s\" every ::%d using 1 title \"%s\"\n" % (files[1],skip_num,title2))
      elif (count == 1):
          PLOT_FILE.write("set title \"%s comparison\"\n" % (files[0]))
          PLOT_FILE.write("plot \"%s\" every ::%d using 1 title \"%s\", \"%s\" every ::%d using 1 title \"%s\"\n" % (files[0],skip_num,title1,files[1],skip_num,title2));
      else:
        for col in range(2,count+1):
          PLOT_FILE.write("set title \"%s comparison of column %d\"\n" % (files[0],col))
          PLOT_FILE.write("plot \"%s\" every ::%d using 1:%d title \"%s\", \"%s\" every ::%d using 1:%d title \"%s\"\n" % (files[0],skip_num,col,title1,files[1],skip_num,col,title2));
    PLOT_FILE.close()
    os.system("rm Discrepancies.*")
    os.system('gnuplot gnuplot_create');
    os.system("ps2pdf Discrepancies.ps Discrepancies.pdf");
  if(sendMail == 1):
    # find out which host mail will be sent from
    p = subprocess.Popen(["hostname"],stdout=subprocess.PIPE)
    hostname = p.stdout.readline()

    msg = MIMEMultipart()
    hostname = hostname.rstrip('\n')
    msg['Subject'] = 'Regression Test results from %s' % hostname
    msg['From'] = 'mpotts@'+hostname
    msg['To'] = "pavery@stanford.edu"
    msg['Date'] = formatdate(localtime=True)
    f = 'Discrepancies.pdf'
    if(os.path.exists("Discrepancies.pdf")):
      msg.attach( MIMEText('These are the discrepancies from the last regression test'))
      part = MIMEBase('application', "octet-stream")
      part.set_payload( open(f,"rb").read() )
      Encoders.encode_base64(part)
      part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(f))
      msg.attach(part)
    else: 
      msg.attach( MIMEText('There were no discrepancies in the last regression test'))
    smtp = smtplib.SMTP("localhost")
    send_to = "pavery@stanford.edu"
    send_from = "mpotts@ahpcrcfe.stanford.edu"
    smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.close()

  sys.exit(result)

if __name__ == "__main__":
  params = sys.argv
  dComp(params)


