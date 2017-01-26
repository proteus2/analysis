#! /bin/bash
# import : TMPDIR

z0=8000   # 0
z9=8000   # 15000
zi=250
zi0=250
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1        ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
vars=( prt_U  prt_V  prt_W  prt_T )
varsd=( U  V  W  T )
iv=2
iwave=5
fotail="__11_00_00__12_00_00__1hrly"
fidir="$DATD/gwjet_wrf"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( n I_WAVE  $iwave )
P102=( s VAR_NAME  ${vars[$iv]} )
P103=( n NNN  289  577  12 )   # Day 7-8
P104=( n XY0  20.0463  21.06481 )
P105=( n LENGTHS  40.  25. )
P106=( n N_KL  217  270 )
P107=( n K_RNG  -999  -999 )    # the first and last zonal wavenumber
P108=( n L_RNG  -999  -999 )    # the first and last meridional wavenumber
P109=( n LAT_C  42. )
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='recon_xy'
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

  P990[2]=$fidir/x${D}_5min/${varsd[$iv]}_z/prt_d300km/fcoef_kl/fft_${vars[$iv]}_kl__z${ZZZZZ}-2a.nc
  ODIR=$fodir/x${D}_5min/${varsd[$iv]}_z/prt_d300km/recon_xy
  P999[2]=$ODIR/${vars[$iv]}__z${ZZZZZ}__x$D${fotail}-w${iwave}.nc
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}" "${P109[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

