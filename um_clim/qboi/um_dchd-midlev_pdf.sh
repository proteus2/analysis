#! /bin/bash
# import : TMPDIR

expname=L60CGW
expcode=uanuj
EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1979
YYYY2=2006
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $expcode )
P104=( n HH  00  8 ) # should be 0 # starting time [UTC] and frequency [/day]
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( n PDF_RNG  0.0  20.e3  81 )
P202=( n LAT_RNG  -10.0  10.0 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  1 )                 # the earlist date-of-file among input files
P902=( n NDAY_I  3 )               # number of days in one input file
P903=( n MISSV  0. )
P990=( s FID  pc  pc  pa  pa  pc )
P991=( s FILE_I_HEAD  "/hippo0/HG2AMIP/$expname" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}a.XXXX_XXXXXXXXXXXX00.nc  -999 )
P993=( s FILE_I_XXXX  FID              FID  YYYY MM  DD        -999 )
P994=( s VAR_I_NAME  dch_max  vars_cgw  u  v  indic_mc )
#P990=( s FID  pb )
#P992=( s FILE_I_FORM  XXXX/${P101[2]}a.XXXX_XXXXXXXXX.nc  -999 )
#P993=( s FILE_I_XXXX  FID              FID  YYYY MM       -999 )
#P994=( s VAR_I_NAME  precip_dc )
P996=( s FILE_ALT  "/hippo0/HG2AMIP/$expname/invariant/L60_z.nc" )
P997=( s VAR_ALT  ht  ht_1 )  # z_theta, z_rho
P999=( s FILE_O  '' )
P998=( s FILE_O2  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year

F_SOURCE='dchd-midlev_pdf_um'
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
  ODIR=$DATD/$expname/dchd_pdf/$YYYY
  P999[2]="$ODIR/${P101[2]}.dchd-midlev_pdf.$YYYY.${MM}.nc"
  P998[2]="$ODIR/${P101[2]}.dchd-nonmidlev_pdf.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.dchd-midlev_pdf.$YYYY.${MM}.nc
#  [ -e ${P998[2]} ] && mv ${P998[2]} $ODIR/old.${P101[2]}.dchd-nonmidlev_pdf.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P996[*]}" "${P997[*]}" "${P999[*]}" "${P998[*]}"
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

