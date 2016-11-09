#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1979
YYYY2=2013
M1=$EXTCTL1
M2=$EXTCTL2
RA_N=5
RA_SNAMES=( era-int  merra  jra55  cfsr  merra2  era-int_f  jra55_f  merra2_f )
#  u, v, w, T
iv=0

RA_SNAME="${RA_SNAMES[$RA_N]}"
. ./predef_ra "$RA_SNAME"
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $RA_SNAME )
P104=( n HH  00  $RA_NHOUR )       # starting time [UTC] and frequency [/day]
#P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
P105=( n REFDATE  1950  1  1 )
#== Parameter 2 ==================================
P201=( s VAR_NAME  ${VAR_I_NAME[$iv]} )
P202=( n H_INTP  0 )
P203=( n K_LARGEST  -999 )
P204=( n LAT_RNG  -30.0  30.0 )
P205=( n P_RNG  500  0.5 )
#== Parameter 9 - I/O ============================
       #        H   D  M  Y
P901=( n NT_F4  ${NT_F4[@]} )
P902=( n P_PREDEF  ${P_PREDEF[@]} )
P903=( n MISSV  ${MISSV[0]} )
P990=( s FILE_I_HEAD  ${FILE_I_HEAD[0]} )
P991=( s FILE_I_FORM  ${FILE_I_FORM[@]} )
P992=( s FILE_I_XXXX  ${FILE_I_XXXX[@]} )
P993=( s VAR_I        ${VAR_I[$iv]}     )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year

F_SOURCE='fft_x_ra'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  MM=$(( M ))  ;  if [ $M -lt 10 ] ; then MM="0$MM" ; fi
  P102[2]=$YYYY
  P103[2]=$MM
  ODIR=$DATD/fspec/$RA_CODENAME/fcoef_k/$YYYY
  P999[2]="$ODIR/${P101[2]}.fft_${VAR_I[$iv]}_k.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P999[*]}"
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

