#! /bin/bash

EXTCTL1=1979     ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=2010
YYYY1=$EXTCTL1
YYYY2=$EXTCTL2
M1=1
M2=12
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  merra )
P104=( n HH  00  8 )               # starting time [UTC] and frequency [/day]
#== Parameter 2 ==================================
P201=( n OPT_AVRG  0 )             # daily-averaging option
P202=( n DZ_EQV_INTP  500.0 )        # [m]
#== Parameter 9 - I/O ============================
       #       H  D  M  Y
P901=( n NT_F4  8  1  1  1 )       # number of time series in one input file
P902=( n missv 1.e15 )             # if no missing points, set 1.0
P990=( s FILE_I_HEAD  '/data11/data-arch/MERRA' )
P991=( s FILE_I_FORM  XXXX  XXXX  merra.XXXX.assm.pl.XXXXXXXXXXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY  MM          VAR_I        YYYY MM  DD      -999 )
P993=( s VAR_I       u  v  omega  t )  # Do not change this order.
P994=( s VAR_I_NAME  u  v  omega  t )
#P994=( s VAR_I_NAME  U_GDS0_ISBL  V_GDS0_ISBL  W_GDS0_ISBL  T_GDS0_ISBL )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year

F_SOURCE='tem_reg'
F_NAMELIST="`echo $PWD`/PROCESS/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2 ] ; do
M=$M1        ;  while [ $M    -le $M2    ] ; do

  MM=$M  ;  if [ $M -lt 10 ] ; then MM="0$M" ; fi
  P102[2]=$YYYY
  P103[2]=$MM
  ODIR=/data17/kyh/dat/qbo/$YYYY/$MM
  P999[2]="$ODIR/${P101[2]}.tem_rg05_ypt.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST &> $F_LOG
  echo "pid : $$"

M=$(( $M+1 ))
done
YYYY=$(( $YYYY+1 ))
done

#== End ==========================================
exit 0

