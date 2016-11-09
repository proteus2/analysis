#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=49       ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
vars=( U  V  W  T )
iv=2
fidir="$DATD/gwjet_wrf"
fistrlev="ml"
fodir="$DATD/gwjet_wrf"

#== Parameter 1 ==================================
P101=( s VAR_NAME  fc_${vars[$iv]}_r \
                   fc_${vars[$iv]}_i )
P102=( n NNN  1  1153  1 )   # the first and last time indices and interval
P103=( n JJ  304  605 )      # the first and last latitudinal indices
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

  P990[2]=$fidir/x${D}_5min/fcoef_k/fft_${vars[$iv]}_k__${fistrlev}${KKK}.nc
  ODIR=$fodir/x${D}_5min/fcoef_kl
  P999[2]=$ODIR/fft_${vars[$iv]}_kl__${fistrlev}${KKK}-1.nc
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

