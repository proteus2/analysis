#!/bin/csh -f

#set exp_name = f2000_cntl
set exp_name = f2000_qbo_uni_test

set exp_root = /grace/s2/kyh/cesm_scr/$exp_name/run
set remapped = $DATD/mjo-qbo/$exp_name

set file_out = $exp_name.qbo.monthly.nc

if (! -d $remapped) mkdir -p $remapped

cd $exp_root
set files = ( *.cam.h0.*.nc )
foreach ff ( $files )
  if (! -f $remapped/umean-$ff ) then
    echo $ff
    ncwa -h -a lon -v U -d lat,-18.,18. -d lev,7.,90. $ff -O $remapped/umean-$ff
  endif
end
cd $remapped
ncrcat -h umean-*.nc -O $file_out

exit

