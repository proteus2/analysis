#! /bin/bash
# import : TMPDIR

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=12       ### CONTROLLED FOR SPLITTING - 2
YYYY1=0001
YYYY2=0005
M1=$EXTCTL1
M2=$EXTCTL2
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  test1_0.9 )
P104=( n HH  00  1 ) # 8 ) # should be 0 # starting time [UTC] and frequency [/day]
P105=( n REFDATE  0001  1  1 )
P106=( n OPT_FEBLEAP  0 )
#== Parameter 2 ==================================
P201=( n DAYS_AVRG   0 )           # averaged days, > 1
#== Parameter 9 - I/O ============================
P901=( n DAY1  1 )                 # the earlist date-of-file among input files
P902=( n NDAY_I  1 )               # number of days in one input file
P903=( n MISSV  1.e32 )            # if no missing points, set 1.0
P990=( s FID  0  0  0  0 )
P991=( s FILE_I_HEAD  "/data11/kyh/analy/waccm_gwp/dat2" )
P992=( s FILE_I_FORM  ${P101[2]}.cam2.hXXXX.XXXX-XXXX.nc  -999 )
P993=( s FILE_I_XXXX                   FID  YYYY  MM      -999 )
       #             U  V  OMEGA  PT
P994=( s VAR_I_NAME  U  V  OMEGA  TH )  # Do not change this order.
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year

F_SOURCE='tem_cam'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh
. $UTIL/shell/ex_date_umstr.sh

#== Body =========================================

Y=$YYYY1  ;  while [ $Y -le $YYYY2 ] ; do
M=$M1     ;  while [ $M -le $M2    ] ; do

  MM=$(( M ))  ;  if [ $M -lt 10 ] ; then MM="0$MM" ; fi
  YYYY=$(( Y ))
  if [ $Y -lt 1000 ] ; then YYYY="0$YYYY" ; fi
  if [ $Y -lt 100  ] ; then YYYY="0$YYYY" ; fi
  if [ $Y -lt 10   ] ; then YYYY="0$YYYY" ; fi
  P102[2]=$YYYY
  P103[2]=$MM
  ODIR=/data11/kyh/analy/waccm_gwp/dat/$YYYY
  P999[2]="$ODIR/${P101[2]}.tem_yzt.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.tem_yzt.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

M=$(( $M+1 ))
done
Y=$(( $Y+1 ))
done

#== End ==========================================
exit 0

