#! /bin/bash

#== Parameter 1 ==================================
P01=( s EXPNAME ctl )
P02=( n YYYYMM  201001 )
#== Parameter 2 ==================================
P03=( n DATE  01  31 )           # the first and the last date
P04=( n UTC  00  12  -999 )      # times specified in UTC
P05=( n FCT  0.5  5.0 )          # forecast time
P06=( n FCT_ITV  0.5 )           # interval

P07=( s EXCEPTION  2009070000 )  # YYYYMMDDTT
P08=( n XY_ITV  4  6 )   # 2.25 x 2.25 deg for N320

FILE_HEAD='/data10/kyh/UM_OPER/UM_OUT'
P00=( s FILE_I "" )
P99=( s FILE_O "" )

F_SOURCE='op-zme2'
F_NAMELIST='namelist/namelist.score-'$$
F_LOG='log/log.score-'$$

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
        'n')  sta=( ${sta[@]} $val',' ) ;;
        's')  sta=( ${sta[@]} \'$val\'',' ) >> $1 ;;
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
cr_nl $F_NAMELIST PARAM "${P03[*]}" "${P04[*]}" "${P05[*]}" "${P06[*]}" "${P07[*]}" "${P08[*]}"

P00[2]="$FILE_HEAD/$exp/$mon"
P99[2]="$FILE_HEAD/$exp/anal/$exp.zme2.$mon.nc"

cr_nl $F_NAMELIST FILEIO "${P00[*]}" "${P99[*]}"
cat $F_NAMELIST

# compile and run --------------------
cp $F_NAMELIST namelist/nl.input
comrun $F_SOURCE > $F_LOG
echo "pid : $$"

#-------------------------------------

done  # mon
done  # exp

#== End ==========================================
exit 0

