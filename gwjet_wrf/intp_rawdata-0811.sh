#! /bin/bash
# import : DIRTMP

fdir="/prime0/kyh"

#== Parameter 1 ==================================
P101=( s VAR_NAMES  U  V  W  T  P  PH  -999 )
P102=( n ND  648  972  49  12 )
P103=( n IT  12  1 )         # the first time indice and number to interpolate
P104=( n NXP  432 )
#== Parameter 9 - I/O ============================
P990=( s FILE_I  $fdir/uvwTppPHP_region_2000-01-08_11:00:00_original.nc \
                 $fdir/uvwTppPHP_region_2000-01-08_12:00:00.nc )
P999=( s FILE_O  $fdir/uvwTppPHP_region_2000-01-08_11:00:00.nc )
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

