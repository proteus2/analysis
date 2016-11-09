#! /bin/bash

#== Parameter 1 ==================================
P01=( s EXPNAME ctl )
P02=( n YYYYMM  201001 )
#== Parameter 2 ==================================
P03=( n DATE  01  31 )            # the first and the last date
P04=( n UTC  00  -999 )           # times specified in UTC
P05=( n FCT  0.0  10.0 )          # (forecast) time
P06=( n FCT_ITV  0.125 )          # sampling interval
P07=( s EXCEPTION  2009070000 )   # YYYYMMDDTT
P08=( n PLEV  1000.  925.  850.  750.  600.  500.  420.  350.  300.  250.   \
               200.  160.  125.  100.   80.   65.   55.   45.   37.   30.   \
                25.   20.   16.   12.5  10.    8.    6.5   5.5   4.5   3.7  \
                 3.    2.5   2.    1.6   1.25  1.    0.8   0.65  0.55  0.45 \
                 0.37  0.3   0.25  0.2   0.16  0.125 0.1   0.07  0.04  -999 )
P09=( n zh_int  12.0 )            # zh_theta(0) < zh_rho(1) < zh_int < zh_theta(1)

FILE_HEAD='/data4/kyh/UM_OPER/L70COLD'
FILE_HEAD2='/data8/kyh/UM_OPER/anal/L70COLD/ssw'
P00=( s FILE_I "" )
P99=( s FILE_O "" )
P98=( s FNAMET  'gwd_uob'  '' )
P97=( s VNAME_I  'field424'  'field425'  'field68'  'field69'  'field1596'  'field1597' )
P96=( s VNAME_O  'BGWD_x'  'BGWD_y'  'OGWD_x'  'OGWD_y'  'OBLD_x'  'OBLD_y' )
P10=( n fill_extp  0.0 )          # value for extrapolated points (1.0 for zero gradient)

FILE_HEAD_P='/data14/kyh/UM_OPER/BACKUP/L70COLD'
P89=( s FILE_P "" )
P88=( s FNAMET_P  'std' )
P87=( s VNAME_P  'p' )            # p_rho

F_SOURCE='ssw-zm_zr_ypft'
F_NAMELIST='namelist/namelist.zm_zr-'$$
F_LOG='log/log.zm_zr-'$$

#== Functions ====================================
function cr_file(){
  if [ -e $1 ] ; then rm -f $1 ; fi ; touch $1
}
function cr_nl(){
  echo >> $1 ; echo ' &'$2 >> $1
  for PI in "${@:3:$(($#-2))}" ; do
    PJ=( $PI )
    echo -n "   ${PJ[1]}  =  " >> $1
    sta=()
    for val in "${PJ[@]:2:$((${#PJ[@]}-2))}" ; do
      case "${PJ[0]}" in
        'n')  sta=( "${sta[@]}" "$val"',' ) ;;
        's')  sta=( "${sta[@]}" \'"$val"\'',' ) >> $1 ;;
        *  )  sta=' INPUT PARAMETER ERROR !!' ; exit 0
      esac
    done
    echo ${sta[@]} >> $1
  done
  echo ' /' >> $1
}

#== Body =========================================

for exp in "${P01[@]:2:$((${#P01[@]}-2))}" ; do
for mon in "${P02[@]:2:$((${#P02[@]}-2))}" ; do

# create namelist --------------------
cr_file $F_NAMELIST
cr_nl $F_NAMELIST ANALCASE "${P01[*]:0:2} $exp" "${P02[*]:0:2} $mon"
cr_nl $F_NAMELIST PARAM "${P03[*]}" "${P04[*]}" "${P05[*]}" "${P06[*]}" "${P07[*]}" \
                        "${P08[*]}" "${P09[*]}" "${P10[*]}"

EXP=`echo $exp | tr '[a-z]' '[A-Z]'`
P00[2]="$FILE_HEAD/$EXP/$mon"
P99[2]="$FILE_HEAD2/$exp.${P98[2]}_p$mon${P98[3]}.nc"
P89[2]="$FILE_HEAD_P/$EXP/$mon"

cr_nl $F_NAMELIST FILEIO "${P00[*]}" "${P99[*]}" "${P98[*]}" "${P97[*]}" "${P96[*]}" \
                         "${P89[*]}" "${P88[*]}" "${P87[*]}"
cat $F_NAMELIST

# compile and run --------------------
cp $F_NAMELIST namelist/nl.input
comrun $F_SOURCE #> $F_LOG
echo "pid : $$"

#-------------------------------------

done  # mon
done  # exp

#== End ==========================================
exit 0

