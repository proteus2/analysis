#! /bin/bash
# import : TMPDIR

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1953
YYYY2=2006
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  uanuj )
P103=( n MM  ''  1 )
P104=( n HH  03  8 )
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
P107=( n NRES  96 )
#== Parameter 2 ==================================
P201=( s VAR_NAME  theta )
P202=( n K_LARGEST  20 )
P203=( n PERIOD_SMALLEST  1 )
P204=( n LAT_RNG  -30.0  30.0 )
P205=( n Z_RNG  8.4  -999 )
P206=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  -999 )              # the earlist date-of-file among input files
P902=( n NDAY_I  30 )              # number of days in one input file
P903=( n MISSV  1.e32 )
P991=( s FILE_I_HEAD  "/prime0/kyh/dat/L60CGW/fcoef_ko" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}.fft_${P201[2]}_ko.XXXX.XXXX.nc  -999 )
P993=( s FILE_I_XXXX  YYYY                              YYYY  MM      -999 )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='wavedcp_var'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh
. $UTIL/shell/ex_date_umstr.sh

#== Body =========================================

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  MM=$(( M ))  ;  if [ $M -lt 10 ] ; then MM="0$MM" ; fi
  P102[2]=$YYYY
  P103[2]=$MM
  ODIR=${P991[2]}/$YYYY
  P999[2]="$ODIR/${P101[2]}.fft_${P201[2]}_wave_ko.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && break

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}"
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

