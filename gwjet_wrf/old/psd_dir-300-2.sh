#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1        ### CONTROLLED FOR SPLITTING - 2
D1=$EXTCTL1
D2=$EXTCTL2
vars=( prt_U  prt_V  prt_W  prt_T )
levs=( 200  300  500  600  850 )
iv=2
il=1
fidir="$DATD/gwjet_wrf"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( s VAR_NAME  fc_${vars[$iv]}_${levs[$il]}_r \
                   fc_${vars[$iv]}_${levs[$il]}_i )
P102=( n K_RNG  0  216 )       # the first and last zonal wavenumber (0,)
P103=( n L_RNG  0  301 )       # the first and last meridional wavenumber (0,)
P104=( n O_RNG  0  431 )       # the first and last frequency (0,)
P105=( n LENGTHS  40.  28.  36. )  # domain lengths in X [deg], Y [deg], T [hr]
P106=( n NPHI   24 )           # number of phi, must be a multiple of 4
P107=( n C_ITV  2.  100. )     # interval and maximum of c [m/s]
P108=( n INV_KH_ITV  10. )     # 1/(interval of Kh) [deg]
P109=( n LAT_C  42. )          # central latitude for c ([deg/hr] to [m/s])
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='psd_dir'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=$D1  ;  while [ $D -le $D2 ] ; do

  P990[2]=$fidir/x$D/prt_d300km/fcoef_klo/fft_${vars[$iv]}_${levs[$il]}_klo-2.nc
  ODIR=$fodir/x$D/prt_d300km/psd
  P999[2]=$ODIR/psd_${vars[$iv]}_${levs[$il]}_phi_2d-2.nc
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}" "${P109[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

D=$(( D + 1 ))
done

#== End ==========================================
exit 0

