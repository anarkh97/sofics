#!/usr/bin/perl

  $src = @ARGV[0];
  $dest = @ARGV[1];
  if($src ne $dest) {
    $command = "cp $src/\*.include $dest";
    system "$command";
    $command = "cp $src/\*.py $dest";
    system "$command";
    $command = "ln -s $src/baseline $dest";
    system "$command";
    $command = "cp $src/\*.contact* $dest";
    system "$command";
    $command = "mkdir $dest/vme4";
    system "$command";
    $command = "cp $src/vme4/control.so $dest/vme4";
    system "$command";
  }  

