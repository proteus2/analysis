#! /bin/bash

#== Functions ====================================
ex_mon_str(){
  local strin cmon
  strin=$1
  cmon=( xxx jan feb mar apr may jun jul aug sep oct nov dec )
  if [ ${#strin} -lt 3 ] ; then
    if [ $strin -lt 1 -o $strin -gt 12 ] ; then return 1 ; fi
    echo ${cmon[$strin]}
  else
    local res imon mmm
    res=-1
    for imon in 1 2 3 4 5 6 7 8 9 10 11 12 ; do
      mmm="${strin:0:3}"
      if [[ "${mmm:0:1}" == [A-Z] ]] ; then mmm="`echo $mmm | tr '[A-Z]' '[a-z]'`" ; fi
      if [ "$mmm" == "${cmon[$imon]}" ] ; then res=$imon ; fi
    done
    if [ $res -eq -1 ] ; then
      return 2
    elif [ $res -lt 10 ] ; then
      res="0$res"
    fi
    echo $res
  fi
}

if [ "`basename $0`" == 'ex_mon_str.sh' -a $# -ne 0 ] ; then ex_mon_str $1 ; fi
