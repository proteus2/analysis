#! /bin/bash

EXTCTL1=2009     ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=2009     ### CONTROLLED FOR SPLITTING - 2
YYYY1=$EXTCTL1
YYYY2=$EXTCTL2
RA_N=0
RA_SNAMES=( era-int  merra  jra55  cfsr  merra2  era-int_f  jra55_f  merra2_f )

RA_SNAME="${RA_SNAMES[$RA_N]}"
. ./predef_ra "$RA_SNAME"
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $RA_SNAME )
P102=( n YYYY  '' )
P103=( n MM  11  1 )               # starting month and number of months per year
P104=( n HH  00  $RA_NHOUR )       # starting time [UTC] and frequency [/day]
#== Parameter 2 ==================================
P201=( n OPT_AVRG  1 )             # daily-averaging option
P202=( n LAT_RNG  -30.0  90.0 )
P203=( n P_RNG  850  30 )
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
#P994=( s VAR_I_NAME  U_GDS0_ISBL  V_GDS0_ISBL  T_GDS0_ISBL  Z_GDS0_ISBL  \
#                     W_GDS0_ISBL )
P995=( s UNIT_H       $UNIT_H )
  #-----
  MM2=$(( ${P103[2]}+${P103[3]}-1 ))
  if [ $MM2 -gt 12 ] ; then MM2=$(( $MM2-12 )) ; fi
  if [ $MM2 -lt 10 ] ; then MM2="0$MM2" ; fi
  MM2="${P103[2]}-$MM2"
  #-----
P999=( s FILE_O  '' )
#=================================================


F_SOURCE='warc3d_qg_ra'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================

YYYY=$YYYY1
while [ $YYYY -le $YYYY2 ] ; do

  P102[2]=$YYYY
#  ODIR=$DATD/warc3d/$RA_CODENAME
  ODIR=/data11/kyh/analy/tem
  P999[2]="$ODIR/${P101[2]}.warc3dqg_ypt.$YYYY.${MM2}_daily.nc"
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

