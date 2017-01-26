#! /bin/bash
# import : TMPDIR

z0=8000   # -250
z9=8000   # 25000
zi=250
zi0=250
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1        ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
iwave=-999
fidir="$DATD/gwjet_wrf"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P100=( n I_WAVE  $iwave )
P101=( s VAR_NAME  fc_prt_W_r  fc_prt_W_i \
                   fc_prt_U_r  fc_prt_U_i \
                   fc_prt_V_r  fc_prt_V_i )
P102=( n NNN  289  577  1 )    # Day 7-8 (D0 = 6)
P103=( n K_RNG  0  216 )       # the first and last zonal wavenumber (0,)
P104=( n L_RNG  0  269 )       # the first and last meridional wavenumber (0,)
P105=( n LENGTHS  40.  25. )   # domain lengths in X, Y [deg]
P106=( n NPHI0  24 )           # number of phi, must be a multiple of 4
P107=( n INV_KH_ITV  18. )     # 1/(interval of Kh) [deg]
P108=( n LAT_C  42. )          # central latitude for k ([/deg] to [/km])
#== Parameter 9 - I/O ============================
P990=( s FILE_I  ''  ''  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='mfsa_wn'
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

  P990[2]=$fidir/x${D}_5min/W_z/prt_d300km/fcoef_kl/fft_prt_W_kl__z${ZZZZZ}-2a.nc
  P990[3]=$fidir/x${D}_5min/U_z/prt_d300km/fcoef_kl/fft_prt_U_kl__z${ZZZZZ}-2a.nc
  P990[4]=$fidir/x${D}_5min/V_z/prt_d300km/fcoef_kl/fft_prt_V_kl__z${ZZZZZ}-2a.nc
  ODIR=$fodir/x${D}_5min/mfs_z/prt_d300km
  P999[2]=$ODIR/mfsa_prt_wn__z${ZZZZZ}-w${iwave}.nc
  [ $iwave -eq 0 -o $iwave -eq -999 ] && P999[2]=$ODIR/mfsa_prt_wn__z${ZZZZZ}-2a.nc
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P100[*]}" "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

