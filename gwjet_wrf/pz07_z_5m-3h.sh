#! /bin/bash
# import : TMPDIR

z0=250    # -250
z9=7000  # 25000
zi=250
zi0=250
#n_com=$(( ( ( z9 - z0 ) / zi ) + 1 ))  ;  echo $n_com  ;  exit

EXTCTL1=1        ### CONTROLLED FOR SPLITTING - 1
EXTCTL2=28       ### CONTROLLED FOR SPLITTING - 2
K1=$EXTCTL1
K2=$EXTCTL2
fidir="$DATD/gwjet_wrf"
fitail="__08_12_00__11_21_00__5min.nc"
fitail0="__08_00_00__12_00_00__5min.nc"
fotail="__09_00_00__11_18_00__3hrly.nc"
vars_f=( U  V  T )  # Do not change the order !
D=6

#== Parameter 1 ==================================
                   # Do not change the order !
P101=( s VAR_NAME  nbe  dpt_dt  dav_dt )
P108=( s VAR_NAME0  mean_U  mean_V  mean_T )
#P102=( n NNN  1  1153  1 )   # the first and last time indices and interval
P102=( n NNN  145  937  36 )     # the first and last time indices and interval
P109=( n NNN0  289  1081  36 )     # the first and last time indices and interval
P103=( n JJ  1  650 )        # the first and last latitudinal indices
P104=( n II  1  432 )        # the first and last longitudinal indices
P105=( n DZ_SHEAR  $(( zi0 * 2 )) )  # [m]
P106=( n DT_TEND   $(( 300 * 2 )) )  # [s]
P107=( n N_SMOOTH9  3 )      # iterating number of the 9-pt local smoothing
#== Parameter 9 - I/O ============================
P990=( s FILE_I  '' )
P991=( s FILE_V_I  ''  '' )
P992=( s FILE_I0  ''  ''  '' )
P999=( s FILE_O  '' )
#=================================================

F_SOURCE='pz07_wrf_z'
F_NAMELIST="$TMPDIR/namelist/namelist.$F_SOURCE-$$"
F_LOG="log/log.$F_SOURCE-$$"

#== Functions ====================================
. $UTIL/shell/cr_namelist.sh

#== Body =========================================
K=$K1
while [ $K -le $K2 ] ; do

  Z=$(( ( ( K - 1 ) * zi ) + z0 ))
  ZZZZZ="$Z"
  [ $Z -lt 10000 ] && ZZZZZ="0$Z"
  [ $Z -lt 1000  ] && ZZZZZ="00$Z"
  [ $Z -lt 100   ] && ZZZZZ="000$Z"
  [ $Z -lt 10    ] && ZZZZZ="0000$Z"
  [ $Z -lt 0     ] && ZZZZZ="_-00${Z:1}"
  [ $Z -lt -10   ] && ZZZZZ="_-0${Z:1}"
  [ $Z -lt -100  ] && ZZZZZ="_-${Z:1}"

  Z1=$(( Z - zi0 ))
  ZZZZZ1="$Z1"
  [ $Z1 -lt 10000 ] && ZZZZZ1="0$Z1"
  [ $Z1 -lt 1000  ] && ZZZZZ1="00$Z1"
  [ $Z1 -lt 100   ] && ZZZZZ1="000$Z1"
  [ $Z1 -lt 10    ] && ZZZZZ1="0000$Z1"
  [ $Z1 -lt 0     ] && ZZZZZ1="_-00${Z1:1}"
  [ $Z1 -lt -10   ] && ZZZZZ1="_-0${Z1:1}"
  [ $Z1 -lt -100  ] && ZZZZZ1="_-${Z1:1}"

  Z2=$(( Z + zi0 ))
  ZZZZZ2="$Z2"
  [ $Z2 -lt 10000 ] && ZZZZZ2="0$Z2"
  [ $Z2 -lt 1000  ] && ZZZZZ2="00$Z2"
  [ $Z2 -lt 100   ] && ZZZZZ2="000$Z2"
  [ $Z2 -lt 10    ] && ZZZZZ2="0000$Z2"
  [ $Z2 -lt 0     ] && ZZZZZ2="_-00${Z2:1}"
  [ $Z2 -lt -10   ] && ZZZZZ2="_-0${Z2:1}"
  [ $Z2 -lt -100  ] && ZZZZZ2="_-${Z2:1}"

  P990[2]=$fidir/x${D}_5min/pz07_z/mean_d300km/pz07terms3__z${ZZZZZ}__x$D$fitail
  P991[2]=$fidir/x${D}_5min/pz07_z/mean_d300km/pz07terms3__z${ZZZZZ1}__x$D$fitail
  P991[3]=$fidir/x${D}_5min/pz07_z/mean_d300km/pz07terms3__z${ZZZZZ2}__x$D$fitail
  i=2
  for ff in ${vars_f[@]} ; do
    P992[$i]=$fidir/x${D}_5min/${ff}_z/prt_d300km/mean_${ff}__z${ZZZZZ}__x$D$fitail0
    i=$(( i + 1 ))
  done
  ODIR=$fidir/x${D}_3hrly/pz07_z/mean_d300km
  P999[2]=$ODIR/pz07f__z${ZZZZZ}__x$D$fotail
  [ -d $ODIR ] || mkdir -p $ODIR

  # create namelist --------------------
  cr_file $F_NAMELIST
  cr_nl $F_NAMELIST PARAM "${P101[*]}" "${P102[*]}" "${P103[*]}" "${P104[*]}" "${P105[*]}" "${P106[*]}" "${P107[*]}" "${P108[*]}" "${P109[*]}"
  cr_nl $F_NAMELIST FILEIO "${P990[*]}" "${P991[*]}" "${P992[*]}" "${P999[*]}"
  cat $F_NAMELIST

  # compile and run --------------------
  comnrun $F_SOURCE $F_NAMELIST # &> $F_LOG
  echo "pid : $$"

K=$(( K + 1 ))
done

#== End ==========================================
exit 0

