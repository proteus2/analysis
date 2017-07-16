#! /bin/bash

EXTCTL1=12       ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=2009
YYYY2=2009
#M12=( 1  12 )
M12=( $EXTCTL1  $EXTCTL2 )  # available only for NM=1
NM=1             # No. of months to calculate for at once
RA_N=0
RA_SNAMES=( era-int  merra  jra55  cfsr  merra2  era-int_f  jra55_f  merra2_f )

RA_SNAME="${RA_SNAMES[$RA_N]}"
. ./predef_ra "$RA_SNAME"
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $RA_SNAME )
P104=( n HH  00  $RA_NHOUR )          # starting time [UTC] and frequency [/day]
#== Parameter 2 ==================================
P201=( n LAT_RNG  -90.0  90.0 )
P202=( n P_RNG  925  10 )
P203=( n K_MAX  15 )
#== Parameter 9 - I/O ============================
P901=( n NT_F4     ${NT_F4[@]}    )
P902=( n MISSV     ${MISSV[0]}    )
#P903=( n P_PREDEF  ${P_PREDEF[@]} )
P990=( s FILE_I_HEAD  ${FILE_I_HEAD[0]} )
P991=( s FILE_I_FORM  ${FILE_I_FORM[@]} )
P992=( s FILE_I_XXXX  ${FILE_I_XXXX[@]} )
P993=( s VAR_I        ${VAR_I[@]:0:2}       ${VAR_I[@]:3:2}       \
                      ${VAR_I[2]} )
P994=( s VAR_I_NAME   ${VAR_I_NAME[@]:0:2}  ${VAR_I_NAME[@]:3:2}  \
                      ${VAR_I_NAME[2]} )
#P994=( s VAR_I_NAME  U  V  T  Z  W )
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  $NM )         # starting month and number of months per year

F_SOURCE='tem3d_qg_nonst-p86_ra'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

[ $NM -ne 1 ] && M12[1]=${M12[0]}

YYYY=$YYYY1  ;  while [ $YYYY -le $YYYY2    ] ; do
M=${M12[0]}  ;  while [ $M    -le ${M12[1]} ] ; do

  MM=$M  ;  if [ $M -lt 10 ] ; then MM="0$M" ; fi
  P102[2]=$YYYY
  P103[2]=$MM
  ODIR=$DATD/tem3d/$RA_CODENAME/$YYYY
  if [ $NM -eq 1 ] ; then
    ODIR=$ODIR/$MM
    MMM="$MM"
  else
    MMM=$(( ${P103[2]}+${P103[3]}-1 ))
    [ $MMM -gt 12 ] && MMM=$(( $MMM-12 ))
    [ $MMM -lt 10 ] && MMM="0$MMM"
    MMM="$MM-$MMM"
  fi
ODIR=/data11/kyh/analy/tem
  P999[2]="$ODIR/${P101[2]}.tem3d_qg_nonst-p86_xyp.$YYYY.$MMM.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
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

