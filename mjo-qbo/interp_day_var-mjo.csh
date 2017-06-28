#!/bin/csh -f

#set exp_name = f2000_cntl
set exp_name = f2000_qbo_uni_test

set exp_root = /grace/s2/kyh/cesm_scr/$exp_name/run
set remapped = $DATD/mjo-qbo/$exp_name

set file_in	 = $exp_name.var-mjo.daily.nc
set file_out = $exp_name.var-mjo.remap.daily.nc

if (! -d $remapped) mkdir -p $remapped

if (! -f $remapped/$file_in) then
  cd $exp_root
  cdo select,name=U200,U850,FLUT *.cam.h1.*.nc $remapped/$file_in
endif

cd $remapped
if (! -f weights.nc) cdo gencon,r144x73 $file_in weights.nc

cdo remap,r144x73,weights.nc $file_in $file_out
exit

#cd /glade/scratch/cyoo/ncl/unicon1
#cdo -monavg wmi_index.f2000_cntl.nc wmi_mon_index.f2000_cntl.nc
#cdo -monavg wmi_index.f2000_qbo_cyclic.nc wmi_mon_index.f2000_qbo_cyclic.nc

