#! /bin/bash
# import : TMPDIR

EXTCTL1=9        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=11       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1980
YYYY2=1980
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  uanuj )
P103=( n MM  ''  3 )               # starting month and number of months
P104=( n HH  03  8 )               # starting time [UTC] and frequency [/day]
#P104=( n HH  03  8  1 )            # for time-averaged data
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( s VAR_NAME  u )
P202=( n PERIOD_SMALLEST  1 )
P203=( n K_LARGEST  -999 )
P204=( n LAT_RNG  -999  -999 )
P205=( n Z_RNG  -999  -999 )
P206=( n NMON_PATCH  0 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  -999 )              # the earlist date-of-file among input files
P902=( n NDAY_I  30 )              # number of days in one input file
P903=( n MISSV  1.e32 )            # if no missing points, set 1.0
P991=( s FILE_I_HEAD  '' )
P992=( s FILE_I_FORM  XXXX/${P101[2]}.fft_${P201[2]}_k.XXXX.XXXX.nc  -999 )
P993=( s FILE_I_XXXX  YYYY                             YYYY  MM      -999 )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='fft_t2_um'
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
  if [ ${P103[3]} -ne 1 ] ; then
    nn=$(( ( P103[3] - 1 )/2 ))
    MM0=$(( M - nn ))
    if [ $MM0 -lt 1 ] ; then YYYY0=$(( YYYY - 1 )) ; MM0=$(( MM0 + 12 )) ; fi
  fi
  P102[2]=$YYYY0
  P103[2]=$MM0
  IDIR=/genie0/kyh/dat/L60CGW-t/fcoef_k
  P991[2]="$IDIR"
  ODIR=/genie0/kyh/dat/L60CGW-t/fcoef_ko/$YYYY
  P999[2]="${P201[2]}"
  [ "${P201[3]}" == "-999" ] || P999[2]="${P201[2]}${P201[3]}"
  P999[2]="$ODIR/${P101[2]}.fft_${P999[2]}_ko.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.fft_${P999[2]}_ko.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}" "${P206[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

M=$(( $M+1 ))
done
YYYY=$(( $YYYY+1 ))
done

#== End ==========================================
exit 0

