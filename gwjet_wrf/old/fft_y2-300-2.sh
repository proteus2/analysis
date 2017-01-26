#! /bin/bash
# import : TMPDIR

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
P102=( n NNN  1  1152  1 )   # the first and last time indices and interval
P103=( n JJ  2  303 )        # the first and last latitudinal indices
P104=( n K_RNG  0  216 )     # the first and last zonal wavenumber
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='fft_y2'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=$D1  ;  while [ $D -le $D2 ] ; do

  P990[2]=$fidir/x$D/prt_d300km/fcoef_k/fft_${vars[$iv]}_${levs[$il]}_k.nc
  ODIR=$fodir/x$D/prt_d300km/fcoef_kl
  P999[2]=$ODIR/fft_${vars[$iv]}_${levs[$il]}_kl-2.nc
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

