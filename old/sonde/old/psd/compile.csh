#!/bin/csh -f
set srcdir = /data/sis/radiosonde/figure/FOT/radioanal/psd
set rundir = /data/sis/radiosonde/figure/FOT/radioanal/psd

mkdir -p $rundir
cd $srcdir
make clean
make
mv psdanal $rundir
