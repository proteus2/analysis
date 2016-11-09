#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1951
YYYY2=2006
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  uantc )
P103=( n MM  ''  1 )
P104=( n HH  03  8 )
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( n K_LARGEST  -999 )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n LAT_RNG  -30.0  30.0 )
P204=( n Z_RNG  8.4  -999 )
#P203=( n LAT_RNG  -999  -999 )
#P204=( n Z_RNG  -999  -999 )
P205=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  1 )                 # the earlist date-of-file among input files
P902=( n NDAY_I  3 )               # number of days in one input file
P990=( s FID  pa  pa  pa )
P991=( s FILE_I_HEAD  "/hippo0/HG2CMIP/L60CGW" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}a.XXXX_XXXXXXXXXXXX00.nc  -999 )
P993=( s FILE_I_XXXX  FID              FID  YYYY MM  DD        -999 )
       #             U  PT     RHO
P994=( s VAR_I_NAME  u  theta  rho )  # Do not change this order.
P995=( s FILE_I_HEAD2  '' )
P996=( s FILE_I_FORM2  ${P101[2]}.fft_XXXXXXXX.XXXX.XXXX.nc  -999 )
P997=( s FILE_I_XXXX2              VAR_I_N AUX YYYY  MM      -999 )
P897=( s NL_AUX2  _ko  _ko_zr  _ko  -wc_ko )
P998=( s VAR_I_NAME2  u  theta  rhov  rhodz_dt )  # Do not change this order.
P999=( s FILE_O  '' )
P1000=( s FILE_O2  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='epf_z_fc2_um'
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
  IDIR2=$DATD/AOL60CGW/fcoef_ko/$YYYY
  P995[2]="$IDIR2"
  ODIR=$DATD/AOL60CGW/epf-wc/$YYYY
  P999[2]="$ODIR/${P101[2]}.epf_koyz.$YYYY.${MM}.nc"
  P1000[2]="$ODIR/${P101[2]}.epf2_koyz.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.epf_koyz.$YYYY.${MM}.nc
#  [ -e ${P1000[2]} ] && mv ${P1000[2]} $ODIR/old.${P101[2]}.epf2_koyz.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P897[*]}" "${P998[*]}" "${P999[*]}" "${P1000[*]}"
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

