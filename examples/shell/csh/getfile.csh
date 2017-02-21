#!/bin/csh
set HOST = nomad2.ncep.noaa.gov
set USER = anonymous
set PASS = kimyh@yonsei.ac.kr
echo $HOST
echo $USER
echo $PASS
cd /data3/kyh/exam
ftp -in $HOST << EOF
user $USER $PASS
cd pub/reanalysis-2/month/dg3
bin 
lcd 2008
get dg3.ft06.200807.avrg.grib
bye
EOF
exit 0
