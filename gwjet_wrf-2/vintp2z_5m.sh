#! /bin/bash
# import : TMPDIR

n0=1
n9=1153
ni=1
#n_com=$(( ( ( n9 - n0 ) / ni ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=1153     ### CONTROLLED FOR SPLITTING - 2
NC1=$EXTCTL1
NC2=$EXTCTL2
vars=( U  V  W  P  T )
vars_f=( U  V  W  PRS  T )
iv=2
fidir="/data18/GW_jet_wrf-2/dat"
fistrlev="ml"
fzstrlev="mls"
fitail="__08_00_00__12_00_00__5min.nc"
fodir="$DATD/gwjet_wrf-2"
fostrlev="z"

#== Parameter 1 ==================================
P101=( s VAR_NAME  ${vars[$iv]} )
P102=( s GRID_I  ${vars[$iv]} )
P103=( n N1  '' )                     # time index
P104=( n KK  1  49 )
P105=( s Z_NAME  PH )
P106=( n Z_DIVIDED_BY  9.81 )         # put g if the data is geopotential
P107=( n Z_OUT  -250  25000  250 )
#P106=( n Z_DIVIDED_BY  1. )           # put g if the data is geopotential
P108=( n H_SCALE  7.2 )               # [km]
#== Parameter 9 - I/O ============================
P990=( s FILE_I_H  '' )
P991=( s FILE_I_T  '' )
P992=( s FILE_Z_H  '' )
P993=( s FILE_Z_T  '' )
P994=( s FILE_D  '' )
P998=( s FILE_O_H  '' )
P999=( s FILE_O_T  '' )
#=================================================

F_SOURCE='vintp2z'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
D=1
if [ "${vars[$iv]}" = "W" ] ; then
  fistrlev="${fistrlev}s"
fi
NC=$NC1
while [ $NC -le $NC2 ] ; do

  NC4=$NC
  [ $NC -lt 1000 ] && NC4="0$NC"
  [ $NC -lt 100  ] && NC4="00$NC"
  [ $NC -lt 10   ] && NC4="000$NC"

  P103[2]=$(( ( ( NC - 1 ) * ni ) + n0 ))
  P990[2]=$fidir/x${D}_5min/${vars_f[$iv]}/${vars_f[$iv]}__${fistrlev}
  P991[2]=__x$D$fitail
  P992[2]=$fidir/x${D}_5min/GEOP/GEOP__${fzstrlev}
  P993[2]=__x$D$fitail
  P994[2]=$fidir/x${D}_5min/domain_x$D.nc
  ODIR=$fodir/x${D}_5min/${vars_f[$iv]}_z
  P998[2]=$ODIR/${vars_f[$iv]}__${fostrlev}
  P999[2]=__x$D$fitail"-$NC4"
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P993[*]}" "${P994[*]}" "${P998[*]}" "${P999[*]}"

  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

NC=$(( NC + 1 ))
done

#== End ==========================================
exit 0

