#! /bin/bash

#== Parameter 1 ==================================
P01=( s EXPNAME  ctl  gwdc )
P02=( n YYYYMM  201001 )
#== Parameter 2 ==================================
P03=( n DATE  01  31 )            # the first and the last date
P04=( n UTC  00  -999 )           # times specified in UTC
P05=( n FCT  1.0  10.0 )          # (forecast) time
P06=( n FCT_INTG  1.0 )           # integral interval
#P07=( n FCT_ITV  0.125 )         # sampling interval
P07=( n NF_INTG  9  9  9  9  9  9  9  9  9  9 )
                                  # sampling # for integral to every fcts
                                  #   = fct_intg/(sampling frequency)+1
P08=( n WN  0  3 )                # the first and the last wavenumber
                                  #   0, for the total perturbation
P09=( s EXCEPTION  2009070000 )   # YYYYMMDDTT

FILE_HEAD='/data4/kyh/UM_OPER/L70COLD'
FILE_HEAD2='/data8/kyh/UM_OPER/anal/L70COLD/ssw'
P00=( s FILE_I "" )
P99=( s FILE_O "" )
P98=( s FNAMET  'std_p'  '' )
P97=( s VNAME_I  'u'  'v'  'dz_dt'  'temp'  'ht' )  # u, v, w, T, gph, respectively

F_SOURCE='ssw-epf_ypft'
F_NAMELIST='namelist/namelist.epf-'$$
F_LOG='log/log.epf-'$$

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
cr_nl $F_NAMELIST PARAM "${P03[*]}" "${P04[*]}" "${P05[*]}" "${P06[*]}" "${P07[*]}" "${P08[*]}" "${P09[*]}"

EXP=`echo $exp | tr '[a-z]' '[A-Z]'`
P00[2]="$FILE_HEAD/$EXP/$mon"
P99[2]="$FILE_HEAD2"

cr_nl $F_NAMELIST FILEIO "${P00[*]}" "${P99[*]}" "${P98[*]}" "${P97[*]}"
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

