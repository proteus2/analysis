#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=49       ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
vars=( U  V  W )
iv=2
fidir="/data18/GW_jet_wrf/dat"
fistrlev="ml"
fitail="__08_00_00__12_00_00__5min.nc"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( s VAR_NAME  ${vars[$iv]} )
P102=( s GRID_I  ${vars[$iv]} )
P103=( n NNN  1  1153  1 )   # the first and last time indices and interval
P104=( n JJ  54  703 )       # the first and last latitudinal indices
P105=( n II  1  432 )        # the first and last longitudinal indices
P106=( n N_DOM  9 )          # No. of the zonally periodic domains
P107=( n D_CRIT  300 )       # great-circle distance [km] for running average
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P991=( s FILE_D  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='prtextr'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=6
K=$K1
if [ "${vars[$iv]}" = "W" -o "${vars[$iv]}" = "GEOP" ] ; then
  K=$(( K1 - 1 ))
  fistrlev="${fistrlev}s"
fi
while [ $K -le $K2 ] ; do

  KKK="$K"
  [ $K -lt 100 ] && KKK="0$K"
  [ $K -lt 10  ] && KKK="00$K"

  P990[2]=$fidir/x${D}_5min/${vars[$iv]}/${vars[$iv]}__${fistrlev}${KKK}__x$D$fitail
  P991[2]=$fidir/x${D}/domain_x$D.nc
  ODIR=$fodir/x${D}_5min/prt_d300km/${vars[$iv]}
  P999[2]=$ODIR/prt_${vars[$iv]}__${fistrlev}${KKK}__x$D$fitail
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P991[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

