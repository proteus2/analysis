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
P201=( n K_LARGEST  -999 )
P202=( n PERIOD_SMALLEST  -999 )
P203=( n LAT_RNG  -30.0  30.0 )
P204=( n P_RNG  500  0.5 )   # should be same with epf_fc2.sh
P205=( n NMON_PATCH  1 )
#== Parameter 9 - I/O ============================
P901=( n NT_F4  ${NT_F4[@]} )
P902=( n P_PREDEF  ${P_PREDEF[@]} )
P990=( s FILE_I_HEAD  "$DATD/tem/$RA_CODENAME/epf" )
P991=( s FILE_I_FORM  XXXX/${P101[2]}.epf_koyz.XXXX.XXXX.nc  -999  \
                      XXXX/${P101[2]}.epf2_koyz.XXXX.XXXX.nc  -999 )
P992=( s FILE_I_XXXX  YYYY                     YYYY  MM      -999 )
P993=( s VAR_I_NAME  f_y  f_z  epd  epd_z  \
                     f_uv  f_uw  epd_uv  epd_uw )
P994=( s FILE_I_HEAD2  ${FILE_I_HEAD[0]} )
P995=( s FILE_I_FORM2  ${FILE_I_FORM[@]} )
P996=( s FILE_I_XXXX2  ${FILE_I_XXXX[@]} )
P997=( s VAR_I2       ${VAR_I[0]}      )  # u
P998=( s VAR_I_NAME2  ${VAR_I_NAME[0]} )
P999=( s FILE_O  ''  '' )
#=================================================
P102=( n YYYY  '' )

F_SOURCE='sepa_epf'
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
  ODIR=${P990[2]}/$YYYY
  P999[2]="$ODIR/${P101[2]}.epf_yz_recon0.$YYYY.${MM}.nc"
  P999[3]="$ODIR/${P101[2]}.epf2_yz_recon0.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  for i in {1..2} ; do

    # create namelist --------------------
    if [ $i -eq 1 ] ; then
      P991n="${P991[*]:0:4}"
      P993n="${P993[*]:0:6}"
      P999n="${P999[*]:0:3}"
    else
      P991n="${P991[*]:0:2} ${P991[*]:4:2}"
      P993n="${P993[*]:0:2} ${P993[*]:6:4}"
      P999n="${P999[*]:0:2} ${P999[3]}"
    fi
    cr_file $F_NAMELIST
    cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
    cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}" "${P204[*]}" "${P205[*]}"
    cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991n[*]}" "${P992[*]}" "${P993n[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P998[*]}" "${P999n[*]}"
    cat $F_NAMELIST

    # compile and run --------------------
    comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
    echo "pid : $$"

  done

M=$(( $M+1 ))
done
YYYY=$(( $YYYY+1 ))
done

#== End ==========================================
exit 0

