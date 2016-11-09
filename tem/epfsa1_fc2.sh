#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1979
YYYY2=2012
M1=$EXTCTL1
M2=$EXTCTL2
RA_N=5
RA_SNAMES=( era-int  merra  jra55  cfsr  merra2  era-int_f  jra55_f  merra2_f )

RA_SNAME="${RA_SNAMES[$RA_N]}"
. ./predef_ra "$RA_SNAME"
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $RA_SNAME )
P103=( n MM  ''  1 )
P104=( n HH  00  $RA_NHOUR )
#P104=( n HH  00  4 )
P105=( n REFDATE  1950  1  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( n K_LARGEST  20 )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n LAT_RNG  -30.0  30.0 )
P204=( n P_RNG  500  0.5 )   # should be same with fft_*.sh
P205=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n NT_F4  ${NT_F4[@]} )
P902=( n P_PREDEF  ${P_PREDEF[@]} )
P903=( s WAVGRP  s  a  -999 )
P991=( s FILE_I_HEAD  ${FILE_I_HEAD[0]} )
P992=( s FILE_I_FORM  ${FILE_I_FORM[@]} )
P993=( s FILE_I_XXXX  ${FILE_I_XXXX[@]} )
       #             U                T
P994=( s VAR_I       ${VAR_I[0]}      ${VAR_I[3]}      )
P995=( s VAR_I_NAME  ${VAR_I_NAME[0]} ${VAR_I_NAME[3]} )
P996=( s FILE_I_HEAD2  "$DATD/fspec/$RA_CODENAME/fcoef_ko" )
P997=( s FILE_I_FORM2  XXXX/${P101[2]}.fft_XXXX_wave_ko.XXXX.XXXX.nc  -999 )
P998=( s FILE_I_XXXX2  YYYY                VAR_I        YYYY  MM      -999 )
       #              U                T                V                OMEGA
P981=( s VAR_I2       ${VAR_I[0]}      ${VAR_I[3]}      ${VAR_I[1]}      ${VAR_I[2]}      )
P982=( s VAR_I_NAME2  ${VAR_I_NAME[0]} ${VAR_I_NAME[3]} ${VAR_I_NAME[1]} ${VAR_I_NAME[2]} )
P999=( s FILE_O  '' )
P1000=( s FILE_O2  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='epfsa1_fc2_ra'
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
  ODIR=$DATD/tem/$RA_CODENAME/epf/$YYYY
  P999[2]="$ODIR/${P101[2]}.epfsa1_koyz.$YYYY.${MM}.nc"
  P1000[2]="$ODIR/${P101[2]}.epf2sa1_koyz.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P998[*]}" "${P981[*]}" "${P982[*]}" "${P999[*]}" "${P1000[*]}"
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

