#! /bin/bash
# import : DIRTMP

z0=0      # -250
z9=25000  # 25000
zi=250
zi0=250
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=101       ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
vars=( U  V  W  T )
iv=2
fidir="$DATD/gwjet_wrf-2"
fitail="__08_00_00__12_00_00__5min.nc"
fodir="$DATD/gwjet_wrf-2"

#== Parameter 1 ==================================
P101=( s VAR_NAME  prt_${vars[$iv]} )
P102=( n NNN  1  1153  1 )   # the first and last time indices and interval
P103=( n JJ  1  650 )        # the first and last latitudinal indices
#P103=( n JJ  54  703 )       # the first and last latitudinal indices
P104=( n II  1  432 )        # the first and last longitudinal indices
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='fft_x'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=1
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

  P990[2]=$fidir/x${D}_5min/${vars[$iv]}_z/prt_d300km/prt_${vars[$iv]}__z${ZZZZZ}__x$D$fitail
  ODIR=$fodir/x${D}_5min/${vars[$iv]}_z/prt_d300km/fcoef_k
  P999[2]=$ODIR/fft_prt_${vars[$iv]}_k__z${ZZZZZ}.nc
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

