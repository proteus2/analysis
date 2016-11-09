#! /bin/bash

#== Functions ====================================
. $UTIL/shell/n10to36.sh
. $UTIL/shell/ex_mon_str.sh

ex_date_umstr(){
  if [ $# -eq 1 ] ; then
    umstr2date $1
  elif [ $# -eq 2 ] || [ $# -eq 4 ] ; then
    date2umstr $@
  else
    echo 'No. of arguments for FUNC:'$FUNCNAME' should be 1, 2 or 4'  ;  return 1
  fi
}

date2umstr(){
  local year month y10 y1 mdh umstr
  year=$(( $1 - 1800 ))  ;  month=$(( $2 ))
  y10=$(( $year / 10 ))
  y1=$(( $year - ( $y10 * 10 ) ))
  if [ $# -eq 2 ] ; then
    mdh="`ex_mon_str $month`"
  else
    local date hour
    date=$(( $3 ))  ;  hour=$(( $4 ))
    mdh="`n10to36 $month`""`n10to36 $date`""`n10to36 $hour`"
  fi
  umstr="`n10to36 $y10`""$y1$mdh"
  if [ ${#umstr} -ne 5 ] ; then return 2 ; fi
  echo "`n10to36 $y10`""$y1$mdh"
}

umstr2date(){
  local strin timstr y10 y1 year month st
  strin=$1
  if [ ${#strin} -lt 5 ] ; then return 3 ; fi
  timstr=${strin:$(( ${#strin} - 5 )):5}
  y10=$(( 36#${timstr:0:1} ))
  y1=${timstr:1:1}
  year=$(( $y10 * 10 + $y1 + 1800 ))
  month=$(( `ex_mon_str ${timstr:2:3}` ))
  st=$?
  if [ $st -eq 0 ] ; then
    if [ $month -lt 10 ] ; then month="0$month" ; fi
    echo $year $month
  else
    local date hour
    month="${timstr:2:1}"
    date="${timstr:3:1}"
    hour="${timstr:4:1}"
    if [[ "$month" > 'c' || "$date" > 'v' || "$hour" > 'o' ]] ; then return 4 ; fi
    month=$(( 13#$month ))
    date=$(( 32#$date ))
    hour=$(( 25#$hour ))
    if [ $month -lt 10 ] ; then month="0$month" ; fi
    if [ $date  -lt 10 ] ; then  date="0$date"  ; fi
    if [ $hour  -lt 10 ] ; then  hour="0$hour"  ; fi
    echo $year $month $date $hour
  fi
}

if [ "`basename $0`" == 'ex_date_umstr.sh' -a $# -ne 0 ] ; then ex_date_umstr $@ ; fi

