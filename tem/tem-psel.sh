#! /bin/bash

EXTCTL1=2009     ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=2009     ### CONTROLLED FOR SPLITTING - 2
YYYY1=$EXTCTL1
YYYY2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  era-int )
P102=( n YYYY  '' )
P103=( n MM  07  1 )               # starting month and number of months per year
P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
P105=( n P_SEL  1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10 5 3 2 1 -999 )
#== Parameter 2 ==================================
P201=( n OPT_AVRG  1 )             # daily-averaging option
#== Parameter 9 - I/O ============================
       #       H   D  M  Y
P901=( n NT_F4  1  30  1  1 )      # number of time series in one input file
P902=( n missv 1.0 )             # if no missing points, set 1.0
P990=( s FILE_I_HEAD  '/data11/data-arch/ERA-inter' )
P991=( s FILE_I_FORM  XXXX  XXXX  era-int.XXXX.anal.XXXX.pl.XXXXXXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY  MM            VAR_I     HH      YYYY MM      -999 )
P993=( s VAR_I       u  v  ome  t )  # Do not change this order.
#P994=( s VAR_I_NAME  u  v  omega  t )
P994=( s VAR_I_NAME  U_GDS0_ISBL  V_GDS0_ISBL  W_GDS0_ISBL  T_GDS0_ISBL )
  #-----
  MM2=$(( ${P103[2]}+${P103[3]}-1 ))
  if [ $MM2 -gt 12 ] ; then MM2=$(( $MM2-12 )) ; fi
  if [ $MM2 -lt 10 ] ; then MM2="0$MM2" ; fi
  MM2="${P103[2]}-$MM2"
  #-----
P999=( s FILE_O  '' )
#=================================================


F_SOURCE='tem-psel'
F_NAMELIST="`echo $PWD`/PROCESS/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1
while [ $YYYY -le $YYYY2 ] ; do

  P102[2]=$YYYY
  ODIR=/data17/kyh/dat/ERAinter
  P999[2]="$ODIR/${P101[2]}.tem-psel_ypt.$YYYY.${MM2}_daily.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST &> $F_LOG
  echo "pid : $$"

  YYYY=$(( $YYYY+1 ))

done

#== End ==========================================
exit 0

