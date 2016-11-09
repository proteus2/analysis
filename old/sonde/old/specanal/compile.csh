#!/bin/csh -f
set srcdir = /data/sis/radiosonde/GWAP/specanal
set rundir = /data/sis/radiosonde/GWAP/run

mkdir -p $rundir
cd $srcdir
make clean
make
mv specanal $rundir
