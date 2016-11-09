#! /bin/bash
# import : DIRTMP

fdir="/data18/GW_jet_wrf/Hist/wrfout_5min"

#== Parameter 1 ==================================
P101=( s VAR_NAMES  U  V  W  T  P  PH  -999 )
P102=( n ND  648  972  49  12 )
P103=( n IT  11  2 )         # the first time indice and number to interpolate
P104=( n NXP  432 )
#== Parameter 9 - I/O ============================
P990=( s FILE_I  $fdir/uvwTppPHP_region_2000-01-10_19:05:00_original.nc \
                 $fdir/uvwTppPHP_region_2000-01-10_20:05:00.nc )
P999=( s FILE_O  $fdir/uvwTppPHP_region_2000-01-10_19:05:00.nc \
                 $fdir/uvwTppPHP_region_2000-01-10_19:05:00.nc )
#=================================================

F_SOURCE='intp_rawdata'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

#== End ==========================================
exit 0

