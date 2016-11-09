#! /bin/bash

#== Parameter 1 ==================================
P01=( s EXPNAME gwdc )
P02=( n YYYYMM  200907 201001 )
#== Parameter 2 ==================================
P03=( n LAT_AVG  2.0 )

FILE_HEAD='/data8/kyh/UM_OPER/anal/L50CYCLE'
FILE_HEAD_u='/data10/kyh/UM_OPER/L50CYCLE'
P00=( s FILE_I "" "" )
P99=( s FILE_O "" )

F_SOURCE='op-wr-dyn_ypft'
F_NAMELIST='namelist/namelist.wr_dyn-'$$
F_LOG='log/log.wr_dyn-'$$

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
cr_nl $F_NAMELIST PARAM "${P03[*]}"

EXP=`echo $exp | tr '[a-z]' '[A-Z]'`
P00[2]="$FILE_HEAD/$EXP/pp"
P00[3]="$FILE_HEAD_u/$EXP/$mon"
P99[2]="$FILE_HEAD/$EXP/pp/$exp.wr-dyn_ypft.$mon.nc"

cr_nl $F_NAMELIST FILEIO "${P00[*]}" "${P99[*]}"
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

