#! /bin/bash
# import : TMPDIR

z0=-250
zi=250
#z9=25000
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=102      ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
vars=( d_nbe_maj  curv  nl_min  ageo  fvor  div  adv )
vars_f=( d_nbe-maj  d_nbe-curv  d_nbe-nl_min  d_nbe-ageo  d_nbe-fvor  div  adv_div )
iv=0
fidir="$DATD/gwjet_wrf/x6_5min/nbe_z/mean_d300km"
fitail="__08_00_00__12_00_00__5min.nc"

#== Parameter 1 ==================================
P101=( s VAR_NAME  ${vars[$iv]} )
P102=( n NNN  1  1153  1 )   # the first and last time indices and interval
P103=( n JJ  33  618 )       # the first and last latitudinal indices
P104=( n II  1  432 )        # the first and last longitudinal indices
P105=( n D_CRIT  300 )       # great-circle distance [km] for running average
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  ''  '' )
#=================================================

F_SOURCE='prtextr_z'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=6
K=$K1
while [ $K -le $K2 ] ; do

  Z=$(( ( ( K - 1 ) * zi ) + z0 ))
  ZZZZZ="$Z"
  [ $Z -lt 10000 ] && ZZZZZ="0$Z"
  [ $Z -lt 1000  ] && ZZZZZ="00$Z"
  [ $Z -lt 100   ] && ZZZZZ="000$Z"
  [ $Z -lt 10    ] && ZZZZZ="0000$Z"
  [ $Z -lt 0     ] && ZZZZZ="_-00${Z:1}"
  [ $Z -lt -10   ] && ZZZZZ="_-0${Z:1}"
  [ $Z -lt -100  ] && ZZZZZ="_-${Z:1}"

  P990[2]=$fidir/${vars_f[$iv]}__z${ZZZZZ}__x$D$fitail
  ODIR=$fidir/prt_d300km
  P999[2]=$ODIR/prt_${vars_f[$iv]}__z${ZZZZZ}__x$D$fitail
  P999[3]=$ODIR/mean_${vars_f[$iv]}__z${ZZZZZ}__x$D$fitail
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

