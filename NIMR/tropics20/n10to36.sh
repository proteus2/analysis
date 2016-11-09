#! /bin/bash

#== Functions ====================================
n10to36(){
  local dig res arg
  dig=( 0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z )
  res=( )
  for arg in $@ ; do
    if [ $arg -lt 0 -o $arg -ge 36 ] ; then return 1 ; fi
    res=( ${res[@]} "${dig[$arg]}" )
  done
  echo "${res[@]}"
}

if [ "`basename $0`" == 'n10to36.sh' -a $# -ne 0 ] ; then n10to36 $@ ; fi

