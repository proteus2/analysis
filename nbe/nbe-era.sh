#! /bin/bash

EXTCTL1=7        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=7
#12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=2011
#1979
YYYY2=2011
#2012
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  era-int )
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year
P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
P105=( n PLEV  400  350  300  250  200  150  -999 )
#== Parameter 2 ==================================
P201=( n GP_DISCONT_LON  0 )       # option for smoothing GP in zonal direction (MERRA needs this.)
#== Parameter 9 - I/O ============================
       #        H  D  M  Y
P901=( n NT_F4  1  30 1  1 )      # number of time series in one input file
P902=( n missv 1.0 )              # if no missing points, set 1.0
P990=( s FILE_I_HEAD  '/data11/data-arch/ERA-inter' )
P991=( s FILE_I_FORM  XXXX  XXXX  ${P101[2]}.XXXX.anal.XXXX.pl.XXXXXXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY  MM               VAR_I      HH     YYYY MM      -999 )
#P993=( s VAR_I       u  v  h  omega )  # Do not change this order.
#P994=( s VAR_I_NAME  u  v  h  omega )
P993=( s VAR_I       u  v  gp  ome )  # Do not change this order.
P994=( s VAR_I_NAME  U_GDS0_ISBL  V_GDS0_ISBL  Z_GDS0_ISBL  W_GDS0_ISBL )
P999=( s FILE_O  '' )
#=================================================


F_SOURCE='nbe_ra'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  P102[2]=$YYYY
  P103[2]=$(( M ))
  MM=$M  ;  if [ $M -lt 10 ] ; then MM="0$MM" ; fi
  #-----
  MM2=$(( ${P103[2]}+${P103[3]}-1 ))
  if [ $MM2 -gt 12 ] ; then MM2=$(( $MM2-12 )) ; fi
  if [ $MM2 -lt 10 ] ; then MM2="0$MM2" ; fi
  MM2="$MM-$MM2"
  #-----
  ODIR=/data11/kyh/analy/nbe
#/prime0/kyh/dat/NBE/era-int
#  P999[2]="$ODIR/${P101[2]}.nbe_xyt.$YYYY.${MM2}.nc"
  P999[2]="$ODIR/${P101[2]}_nbe_$YYYY${MM}_daily.nc"

  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

M=$(( M+1 ))
done
YYYY=$(( $YYYY+1 ))
done

#== End ==========================================
exit 0

