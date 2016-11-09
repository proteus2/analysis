#! /bin/bash
# import : DIRTMP

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=1979
YYYY2=2012
M1=$EXTCTL1
M2=$EXTCTL2
RA_N=3
RA_SNAMES=( era-int  merra  jra55  era-int_f  cfsr )

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
P201=( n K_LARGEST  -999 )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n LAT_RNG  -30.0  30.0 )
P204=( n P_RNG  500  -999 )
#P203=( n LAT_RNG  -999  -999 )
#P204=( n P_RNG  -999  -999 )
P205=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n NT_F4  ${NT_F4[@]} )
P902=( n P_PREDEF  ${P_PREDEF[@]} )
P990=( s FILE_I_HEAD  "$DATD/tem/$RA_CODENAME/epf" )
P991=( s FILE_I_FORM  XXXX/${P101[2]}.epXXXX_koyz.XXXX.XXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY              AUX       YYYY  MM      -999 )
P892=( s NL_AUX  fsa0  fsa0  fsa0  fsa0 \
                 f2sa0  f2sa0  f2sa0  f2sa0 \
                 fsa0  fsa0  fsa0  fsa0 \
                 f2sa0  f2sa0  f2sa0  f2sa0 \
                 fsa1  fsa1  fsa1  fsa1 \
                 f2sa1  f2sa1  f2sa1  f2sa1 \
                 fsa1  fsa1  fsa1  fsa1 \
                 f2sa1  f2sa1  f2sa1  f2sa1 )
P993=( s VAR_I_NAME  f_y_s  f_z_s  epd_s  epd_z_s \
                     f_uv_s  f_uw_s  epd_uv_s  epd_uw_s \
                     f_y_a  f_z_a  epd_a  epd_z_a \
                     f_uv_a  f_uw_a  epd_uv_a  epd_uw_a \
                     f_y_s  f_z_s  epd_s  epd_z_s \
                     f_uv_s  f_uw_s  epd_uv_s  epd_uw_s \
                     f_y_a  f_z_a  epd_a  epd_z_a \
                     f_uv_a  f_uw_a  epd_uv_a  epd_uw_a )  # Do not change this order.
P994=( s FILE_I_HEAD2  "$DATD/fspec/$RA_CODENAME/fcoef_ko" )
P995=( s FILE_I_FORM2  XXXX/${P101[2]}.fft_XXXX_wave_ko.XXXX.XXXX.nc  -999 )
P996=( s FILE_I_XXXX2  YYYY                VAR_I        YYYY  MM      -999 )
P997=( s VAR_I2       ${VAR_I[0]}  ${VAR_I[3]}  ${VAR_I[1]}  ${VAR_I[2]} )
P998=( s VAR_I_NAME2  ${VAR_I_NAME[0]}  ${VAR_I_NAME[3]}  ${VAR_I_NAME[1]}  ${VAR_I_NAME[2]} )  # Do not change this order.
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='reconstr_varsa'
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
  ODIR=$DATD/fspec/$RA_CODENAME/fcoef_ko-wav4/$YYYY
  P999[2]="$ODIR/${P101[2]}.fft_var4_ko.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P892[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P998[*]}" "${P999[*]}"
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

