#!/usr/bin/perl

  $hostfile = $ENV{'PBS_NODEFILE'};
  print "hostfile = $hostfile\n";
  system("rm host.*");
  @machs = split(/\n/,`cat $hostfile`);
  $len = @machs;
  $j = 0;
  $procs_per_job = 4;
  for($i = 0; $i < $len; $i++) {
    if(($i % $procs_per_job) == 0) {
      $filename = ">host.".$j;
      open(OFILE,$filename);
      $j++;
    }
    print OFILE "@machs[$i] \n";
    if(($i % $procs_per_job) == ($procs_per_job - 1)) { close OFILE; }
    
  }
