#! /bin/bash
# import : TMPDIR

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1953
#1953
YYYY2=1955
#2006
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  uanuj )
P103=( n MM  ''  1 )               # starting month and number of months
P104=( n HH  03  8 )               # starting time [UTC] and frequency [/day]
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( s VAR_NAME  u )
P202=( s U_INTP  -999 )   # set "theta" when interpolating U to theta levels
                          ## set "zr" or "zt" when interpolating each model-
                          ##   level data to its pure altitudes.
P203=( n PERIOD_SMALLEST  -999 )
P204=( n K_LARGEST  -999 )
P205=( n LAT_RNG  -999  -999 )
P206=( n Z_RNG  -999  -999 )
P207=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  -999 )              # the earlist date-of-file among input files
P902=( n NDAY_I  30 )              # number of days in one input file
P903=( n MISSV  1.e32 )            # if no missing points, set 1.0
P991=( s FILE_I_HEAD  '' )
P992=( s FILE_I_FORM  XXXX/${P101[2]}.fft_${P201[2]}_k.XXXX.XXXX.nc  -999 )
P993=( s FILE_I_XXXX  YYYY                             YYYY  MM      -999 )
P994=( s VAR_I_NAME   fcr_${P201[2]}  fci_${P201[2]} )
P995=( s FILE_I_HEAD2  '' )
P996=( s FILE_I_FORM2  XXXX/${P101[2]}.fft_u_k.XXXX.XXXX.nc  -999 )
P997=( s FILE_I_XXXX2  YYYY                    YYYY  MM      -999 )
P998=( s VAR_I_NAME2   fcr_u )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='fft_t2i_um'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh
. $UTIL/shell/ex_date_umstr.sh

#== Body =========================================

[ ${P207[2]} -ne 0 ] && P103[3]=$(( P103[3] + ( 2 * P207[2] ) ))

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  MM=$(( M ))  ;  if [ $MM -lt 10 ] ; then MM="0$MM" ; fi
  MM0=$M  ;  YYYY0=$YYYY
  if [ ${P207[2]} -ne 0 ] ; then
    MM0=$(( M - P207[2] ))
    if [ $MM0 -lt 1 ] ; then YYYY0=$(( YYYY - 1 )) ; MM0=$(( MM0 + 12 )) ; fi
  fi
  P102[2]=$YYYY0
  P103[2]=$MM0
  IDIR=/prime0/kyh/dat/L60CGW/fcoef_k
  P991[2]="$IDIR"
  P995[2]="${P991[2]}"
  ODIR=/prime0/kyh/dat/L60CGW/fcoef_koi/$YYYY
  P999[2]="${P201[2]}"
  [ "${P201[3]}" == "-999" ] || P999[2]="${P201[2]}${P201[3]}"
  P999[2]="$ODIR/${P101[2]}.fft_${P999[2]}_koi.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.fft_${P999[2]}_koi.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}" "${P206[*]}" "${P207[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P998[*]}" "${P999[*]}"
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

