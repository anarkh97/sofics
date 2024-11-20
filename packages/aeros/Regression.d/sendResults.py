#!/usr/bin/python
"""Tool designed to compare data files produced by regression testing and report discrepancies"""

__author__ = "Mark A. Potts (mpotts@drc.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2012/02/2"

import sys, os, re, md5, subprocess, math, glob, datetime, time


def sResults(params):
  import time, glob
  import subprocess

  mail_file = open(".mail_command","w")
  next_build_num = list(open("/lustre/home/hudson/jobs/FEM build/nextBuildNumber","r").read().splitlines())
  build_num = int(next_build_num[0])-1
  log_file = "/lustre/home/hudson/jobs/FEM\ build/builds/" + str(build_num) +"/log"
  command = "tail -90 " + log_file + "| mail -s \"Regression Test Summary\" mpotts@hpti.com & "
  mail_file.write(command)
  mail_file.close()
  os.system(".mail_command")

if __name__ == "__main__":
  params = sys.argv
  sResults(params)


