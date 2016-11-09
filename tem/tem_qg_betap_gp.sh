#! /bin/bash

EXTCTL1=1983     ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1983     ### CONTROLLED FOR SPLITTING - 2
YYYY1=$EXTCTL1
YYYY2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  era-int )
P102=( n YYYY  '' )
P103=( n MM  11  5 )               # starting month and number of months per year
P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
#== Parameter 2 ==================================
P201=( n OPT_AVRG  1 )             # daily-averaging option
P202=( n LAT0  45 )
P203=( n LAT_RNG  10  80 )
#== Parameter 9 - I/O ============================
       #        H   D  M  Y
P901=( n NT_F4  1  30  1  1 )     # number of time series in one input file
P902=( n missv  1.0 )             # if no missing points, set 1.0
P990=( s FILE_I_HEAD  '/data11/data-arch/ERA-inter' )
P991=( s FILE_I_FORM  XXXX  XXXX  era-int.XXXX.anal.XXXX.pl.XXXXXXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY  MM            VAR_I     HH      YYYY MM      -999 )
P993=( s VAR_I       gp  v  ome )  # Do not change this order.
P994=( s VAR_I_NAME  Z_GDS0_ISBL  V_GDS0_ISBL  W_GDS0_ISBL )
P995=( l l_gpheight  False )
       #        H  D  M  Y
#P901=( n NT_F4  8  1  1  1 )      # number of time series in one input file
#P902=( n missv  1.e20 )           # if no missing points, set 1.0
#P990=( s FILE_I_HEAD  "/data11/data-arch/MERRA" )
#P991=( s FILE_I_FORM  XXXX/XXXX/${P101[2]}.XXXX.assm.pl.XXXXXXXXXXXX.nc  -999 )
#P992=( s FILE_I_XXXX  YYYY  MM             VAR_I        YYYY MM  DD      -999 )
#P993=( s VAR_I       h  v  omega )
#P994=( s VAR_I_NAME  h  v  omega )
#P995=( l l_gpheight  True )
  #-----
  MM2=$(( ${P103[2]}+${P103[3]}-1 ))
  if [ $MM2 -gt 12 ] ; then MM2=$(( $MM2-12 )) ; fi
  if [ $MM2 -lt 10 ] ; then MM2="0$MM2" ; fi
  MM2="${P103[2]}-$MM2"
  #-----
P999=( s FILE_O  '' )
#=================================================


F_SOURCE='tem_qg_betap_gp'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1
while [ $YYYY -le $YYYY2 ] ; do

  P102[2]=$YYYY
  ODIR=/data11/kyh/analy/tem
#/data17/kyh/dat/ERAinter
  P999[2]="$ODIR/${P101[2]}.temqgb_ypt.$YYYY.${MM2}_daily.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P995[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

  YYYY=$(( $YYYY+1 ))

done

#== End ==========================================
exit 0

