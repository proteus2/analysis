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
#P104=( n HH  03  8 )               # starting time [UTC] and frequency [/day]
P104=( n HH  00  4 )               # starting time [UTC] and frequency [/day]
P105=( n REFDATE  1941  12  1 )
P106=( n OPT_30D  1 )
#== Parameter 2 ==================================
P201=( n P_INTP  104.914  85.7301  70.0952  57.5617  47.4395  39.2305  \
                 32.5301  27.0280  22.4926  18.7365  15.6013  12.9741  \
                 10.7673  8.91081  7.34516  6.02524  4.91298  3.97635  \
                 3.18840  2.52759  1.97616  1.51997  1.14672  0.843766 \
                 0.601832  0.409983  -999 )
#P201=( n P_INTP  115.818 95.0034 77.5195 63.5202 52.2569 43.1419  \
#                 35.7266 29.6562 24.6623 20.5366 17.1065 14.2379  \
#                 11.8312 9.80800 8.10368 6.66632 5.45450 4.43332  \
#                 3.57344 2.85078 2.24585 1.74286 1.32869 0.990873  \
#                 0.718601 0.501625 0.332068  -999 )
P202=( n LAT_RNG  -15.0  15.0 )
P203=( n DAYS_AVRG  0 )            # averaged days, non-negative integer
#== Parameter 9 - I/O ============================
P901=( n DAY1  1 )                 # the earlist date-of-file among input files
P902=( n NDAY_I  3 )               # number of days in one input file
P903=( n MISSV  1.e30 )            # if no missing points, set 1.0
P990=( s FID  pc  pa  pa  pa )
P991=( s FILE_I_HEAD  "/hippo0/HG2AMIP/$expname" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}a.XXXX_XXXXXXXXXXXX00.nc  -999 )
P993=( s FILE_I_XXXX  FID              FID  YYYY MM  DD        -999 )
                     # The first should be the pressure (at rho level).
                     # omega requires Temperature.
#
P994=( s VAR_I_NAME  p  u  v  -999 )  # u, v, or ht following p
ovar="uv"
P995=( s VAR_I_VERT_GRID  rho )    # rho / theta
#
#P994=( s VAR_I_NAME  p  theta  dz_dt  -999 )  # theta and/or dz_dt following p
#ovar="tw"
#P995=( s VAR_I_VERT_GRID  theta )    # rho / theta
#
P996=( s FILE_ALT  "/hippo0/HG2AMIP/$expname/invariant/L60_z.nc" )
P997=( s VAR_ALT  ht  ht_1 )  # z_theta, z_rho
P999=( s FILE_O  '' )
#=================================================
P102=( n YYYY  '' )
P103=( n MM  ''  1 )               # starting month and number of months per year

F_SOURCE='qboi_vintp_z2p_um'
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
  ODIR=/prime3/kyh/QBOi_upload/$expname/$YYYY
  P999[2]="$ODIR/${P101[2]}.${ovar}_xypt.$YYYY.${MM}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi
#  [ -e ${P999[2]} ] && mv ${P999[2]} $ODIR/old.${P101[2]}.tem_yzt.$YYYY.${MM}.nc

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}"
  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}" "${P203[*]}"
  cr_nl $F_NAMELIST FILEIO "${P901[*]}" "${P902[*]}" "${P903[*]}" "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P995[*]}" "${P996[*]}" "${P997[*]}" "${P999[*]}"
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

