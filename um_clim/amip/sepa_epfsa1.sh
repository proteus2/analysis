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
#== Parameter 2 ==================================
P201=( n K_LARGEST  -999 )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n LAT_RNG  -30.0  30.0 )
P204=( n Z_RNG  8.4  -999 )
#P203=( n LAT_RNG  -999  -999 )
#P204=( n Z_RNG  -999  -999 )
P205=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n DAY1  -999 )              # the earlist date-of-file among input files
P902=( n NDAY_I  30 )              # number of days in one input file
P990=( s FID  pj )
P991=( s FILE_I_HEAD  "/prime0/kyh/dat/L60CGW/epf-wc" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}.epXXXX_koyz.XXXX.XXXX.nc  -999 )
P993=( s FILE_I_XXXX  YYYY              AUX       YYYY  MM      -999 )
P893=( s NL_AUX  fsa0  fsa0  fsa0  fsa0 \
                 f2sa0  f2sa0  f2sa0  f2sa0 \
                 fsa0  fsa0  fsa0  fsa0 \
                 f2sa0  f2sa0  f2sa0  f2sa0 \
                 fsa1  fsa1  fsa1  fsa1 \
                 f2sa1  f2sa1  f2sa1  f2sa1 \
                 fsa1  fsa1  fsa1  fsa1 \
                 f2sa1  f2sa1  f2sa1  f2sa1 )
P994=( s VAR_I_NAME  f_y_s  f_z_s  epd_s  epd_z_s \
                     f_uv_s  f_uw_s  epd_uv_s  epd_uw_s \
                     f_y_a  f_z_a  epd_a  epd_z_a \
                     f_uv_a  f_uw_a  epd_uv_a  epd_uw_a \
                     f_y_s  f_z_s  epd_s  epd_z_s \
                     f_uv_s  f_uw_s  epd_uv_s  epd_uw_s \
                     f_y_a  f_z_a  epd_a  epd_z_a \
                     f_uv_a  f_uw_a  epd_uv_a  epd_uw_a )  # Do not change this order.
P995=( s FILE_I_HEAD2  "/hippo0/HG2AMIP/L60CGW" )
P996=( s FILE_I_FORM2  XXXX/${P101[2]}a.XXXX_XXXXXXXX.nc  -999 )
P997=( s FILE_I_XXXX2  FID              FID  YYYY MM      -999 )
P998=( s VAR_I_NAME2  u )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='sepa_epfsa1'
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
  P999[2]="$ODIR/${P101[2]}.epfsa_yz_recon1.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.epfsa0_yz_recon0.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P893[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P998[*]}" "${P999[*]}"
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

