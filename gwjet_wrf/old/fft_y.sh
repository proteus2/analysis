#! /bin/bash
# import : TMPDIR

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1        ### CONTROLLED FOR SPLITTING - 2
D1=$EXTCTL1
D2=$EXTCTL2
vars=( U  V  W  T )
levs=( 200  300  500  600  850 )
iv=2
il=1
fidir="/data18/GW_jet_wrf/dat"
fitail="__08_00_05__12_00_00.nc"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( s VAR_NAME  ${vars[$iv]}_${levs[$il]} )
P102=( n NNN  1  1152  1 )   # the first and last time indices and interval
P103=( n JJ  2  649 )        # the first and last latitudinal indices
P104=( n II  1  432 )        # the first and last longitudinal indices
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='fft_y'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=$D1  ;  while [ $D -le $D2 ] ; do

  P990[2]=$fidir/x$D/${vars[$iv]}/${vars[$iv]}_${levs[$il]}__x$D$fitail
  ODIR=$fodir/x$D/fcoef_l
  P999[2]=$ODIR/fft_${vars[$iv]}_${levs[$il]}_l.nc
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

D=$(( D + 1 ))
done

#== End ==========================================
exit 0

