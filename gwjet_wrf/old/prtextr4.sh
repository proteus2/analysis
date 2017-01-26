#! /bin/bash
# import : TMPDIR

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1        ### CONTROLLED FOR SPLITTING - 2
D1=$EXTCTL1
D2=$EXTCTL2
vars=( U  V  W  T )
iv=0
fidir="/data18/GW_jet_wrf/dat"
fitail="__08_00_00__12_00_00__6hrly.nc"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( s VAR_NAME  ${vars[$iv]} )
P102=( s GRID_I  ${vars[$iv]} )
P103=( n NNN  1  17  1 )     # the first and last time indices and interval
P104=( n KK  1  49 )         # the first and last time indices and interval
P105=( n JJ  54  703 )       # the first and last latitudinal indices
P106=( n II  1  432 )        # the first and last longitudinal indices
P107=( n N_DOM  9 )          # No. of the zonally periodic domains
P108=( n D_CRIT  300 )       # great-circle distance [km] for running average
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P991=( s FILE_D  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='prtextr4'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=$D1  ;  while [ $D -le $D2 ] ; do

  P990[2]=$fidir/x$D/${vars[$iv]}/${vars[$iv]}__x$D$fitail
  P991[2]=$fidir/x$D/domain_x$D.nc
  ODIR=$fodir/x$D/prt_d300km
  P999[2]=$ODIR/prt_${vars[$iv]}__x$D$fitail
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P991[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

D=$(( D + 1 ))
done

#== End ==========================================
exit 0

