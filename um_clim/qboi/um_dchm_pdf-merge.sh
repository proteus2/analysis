#! /bin/bash
# import : TMPDIR

expname=L60CGW
expcode=uanuj
#fmidname=dchm_pdf
#fmidname=dchm-midlev_pdf
fmidname=dchm-nonmidlev_pdf
#== Parameter 1 - CASE ===========================
P101=( s EXPNAME  $expcode )
P102=( n YYYY2  1979  2006 )
P103=( n MM  1  12 )               # starting month and number of months per year
#== Parameter 2 ==================================
#== Parameter 9 - I/O ============================
P991=( s FILE_I_HEAD  "$DATD/$expname/dchm_pdf" )
P992=( s FILE_I_FORM  XXXX/${P101[2]}.$fmidname.XXXX.XXXX.nc  -999 )
P993=( s FILE_I_XXXX  YYYY                      YYYY  MM      -999 )
P994=( s VAR_I_NAME  N_pop  \
                     dchmax  zcba    zcta    rho_ct  n_q     n_ct    \
                     t_ct    cq_x    cq_y    u_ct    v_ct    u_sfc   \
                     v_sfc   u_cb    v_cb    -999 )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='dchm_pdf-merge_um'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh
. $UTIL/shell/ex_date_umstr.sh

#== Body =========================================

  M2=${P103[2]}
  M3=$(( M2 + P103[3] - 1 ))  ;  [ $M3 -gt 12 ] && M3=$(( M3 - 12 ))
  [ $M2 -lt 10 ] && M2=0$M2  ;  [ $M3 -lt 10 ] && M3=0$M3

  ODIR=$DATD/$expname/dchm_pdf
  P999[2]="$ODIR/${P101[2]}.$fmidname.${P102[2]}-${P102[3]}.${M2}-${M3}.nc"
  if [ ! -d $ODIR ] ; then mkdir -p $ODIR ; fi

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST ANALCASE "${P101[*]}" "${P102[*]}" "${P103[*]}"
#  cr_nl $F_NAMELIST PARAM "${P201[*]}" "${P202[*]}"
  cr_nl $F_NAMELIST FILEIO "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

#== End ==========================================
exit 0

