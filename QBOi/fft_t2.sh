#! /bin/bash
# import : TMPDIR

EXTCTL1=11       ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=11       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1994
YYYY2=1994
M1=$EXTCTL1
M2=$EXTCTL2
#  u, v, w, T
iv=0

RA_SNAME=era-int_f
. ./predef_ra "$RA_SNAME"
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $RA_SNAME )
P103=( n MM  ''  1 )               # starting month and number of months
P104=( n HH  00  $RA_NHOUR )       # starting time [UTC] and frequency [/day]
#P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
P105=( n REFDATE  1950  1  1 )
#== Parameter 2 ==================================
P201=( s VAR_NAME  ${VAR_I_NAME[$iv]} )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n K_LARGEST  -999 )
P204=( n LAT_RNG  -999  -999 )
P205=( n P_RNG  -999  -999 )
P206=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
       #        H   D  M  Y
P901=( n NT_F4  $RA_NHOUR  30  1  1 )      # number of time series in one input file
#P901=( n NT_F4  4  30  1  1 )       # number of time series in one input file
P902=( n MISSV  1.0 )              # if no missing points, set 1.0
P990=( s FILE_I_HEAD  "$DATD/QBOi/fspec/fcoef_k" )
P991=( s FILE_I_FORM  XXXX/${P101[2]}.fft_${VAR_I[$iv]}_k.XXXX.XXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY                                YYYY  MM      -999 )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='fft_t2_ra'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh
. $UTIL/shell/ex_date_umstr.sh

#== Body =========================================

[ ${P206[2]} -ne 0 ] && P103[3]=$(( P103[3] + ( 2 * P206[2] ) ))

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  MM=$(( M ))  ;  if [ $MM -lt 10 ] ; then MM="0$MM" ; fi
  MM0=$M  ;  YYYY0=$YYYY
  if [ ${P206[2]} -ne 0 ] ; then
    MM0=$(( M - P206[2] ))
    if [ $MM0 -lt 1 ] ; then YYYY0=$(( YYYY - 1 )) ; MM0=$(( MM0 + 12 )) ; fi
  fi
  P102[2]=$YYYY0
  P103[2]=$MM0
  ODIR=$DATD/fspec/$RA_CODENAME/fcoef_ko/$YYYY
  P999[2]="$ODIR/${P101[2]}.fft_${VAR_I[$iv]}_ko.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}" "${P206[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun2 $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

M=$(( $M+1 ))
done
YYYY=$(( $YYYY+1 ))
done

#== End ==========================================
exit 0

