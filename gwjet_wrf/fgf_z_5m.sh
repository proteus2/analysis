#! /bin/bash
# import : TMPDIR

z0=0      # -250
z9=15000  # 25000
zi=250
zi0=250
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=61       ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
fidir="$DATD/gwjet_wrf"
fitail="__08_00_00__12_00_00__5min.nc"
vars_f=( U  V  T )  # Do not change the order !
D=6
ODIR=$fidir/x${D}_5min/fgf_z/mean_d300km

#== Parameter 1 ==================================
                   # Do not change the order !
P101=( s VAR_NAME  mean_U  mean_V  mean_T )
P102=( n NNN  1  1153  1 )   # the first and last time indices and interval
P103=( n JJ  1  650 )        # the first and last latitudinal indices
P104=( n II  1  432 )        # the first and last longitudinal indices
#== Parameter 9 - I/O ============================
P990=( s FILE_I  ''  ''  '' )
P998=( s FILE_O_D  $ODIR )
P999=( s FILE_O_T  '' )
#=================================================

F_SOURCE='fgf_wrf_z'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
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

  i=2
  for ff in ${vars_f[@]} ; do
    P990[$i]=$fidir/x${D}_5min/${ff}_z/prt_d300km/mean_${ff}__z${ZZZZZ}__x$D$fitail
    i=$(( i + 1 ))
  done
  P999[2]=__z${ZZZZZ}__x$D$fitail
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P998[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

