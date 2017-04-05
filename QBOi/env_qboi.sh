#!/bin/bash

##  ./env_qboi.sh 1 ua mon YPT r1i1p1
##  ./env_qboi.sh 5 ta 6hr XYPT r2i1p1 200011

# PARAMETER  :::::::::::::::::::::::::::::::::::::
dir_qboi="/grace/s3/kyh/dat/QBOi"   # archive directory
expno=${expno:=$1}       # 1 / 5
var=${var:=$2}           # ta / ua / ...
timeint=${timeint:=$3}   # mon / day / 6hr
diminfo=${diminfo:=$4}   # XYPT / XYT / YPT / ...
ripcode=${ripcode:=$5}   # r1i1p1 / ...

if [ "$expno" = "5" ] ; then
  if [ $# -ne 6 ] ; then echo "PARAMETER NOT SET:  exp5i" ; exit 5 ; fi
  exp5i=${exp5i:=$6}   # 200011 / ...
fi
#-------------------------------------------------

# STRING FOR THE EXP 5  ::::::::::::::::::::::::::
year5i=""  ;  mon5i=""  ;  mon5i_s=""
if [ "$expno" = "5" ] ; then
  year5i=${exp5i:0:4}
  mon5i=${exp5i:4:2}
  if [ "$mon5i" = "05" ] ; then
    mon5i_s=May
  else
    mon5i_s=Nov
  fi
fi
#-------------------------------------------------

# DIRECTORY  :::::::::::::::::::::::::::::::::::::
dir_cmam="$dir_qboi/CCCma/CMAM"
dir_miroc_a_ll="$dir_qboi/MIROC/MIROC-AGCM/MIROC-AGCM-LL"
dir_miroc_e="$dir_qboi/MIROC/MIROC-ESM"
dir_mri="$dir_qboi/MRI"
#-------------------------------------------------
dir_cmam="$dir_cmam/QBOiExp$expno$mon5i_s"
dir_miroc_a_ll="$dir_miroc_a_ll/QBOiExp$expno"
dir_miroc_e="$dir_miroc_e/QBOiExp$expno"
dir_mri="$dir_mri/QBOiExp${expno}"
if [ "$expno" = "5" ] ; then
  dir_mri="${dir_mri}_QBOhindcast"
elif [ "$expno" = "1" ] ; then
  dir_mri="${dir_mri}_amip"
fi
#-------------------------------------------------
dir_cmam="$dir_cmam/$timeint/atmos/$var/$ripcode"
dir_miroc_a_ll="$dir_miroc_a_ll/$timeint/atmos/$var/${mon5i_s}${year5i}*$ripcode"
dir_miroc_e="$dir_miroc_e/$timeint/atmos/$var/${mon5i_s}${year5i}*$ripcode"
dir_mri="$dir_mri/atm$diminfo$timeint"
#-------------------------------------------------

# FILE  ::::::::::::::::::::::::::::::::::::::::::
dimcode="$timeint"
[ "$timeint" = "mon" ] && dimcode="A$dimcode"
[ "${diminfo:0:1}" != "X" ] && dimcode="Z$dimcode"
[ "$timeint" = "6hr" ] && dimcode="${dimcode}P?ev"
#-------------------------------------------------
files_cmam=( $dir_cmam/${var}_${dimcode}_CMAM_QBOiExp${expno}${mon5i_s}_${ripcode}_$year5i${mon5i}*.nc )
files_miroc_a_ll=( $dir_miroc_a_ll/${var}_${dimcode}_MIROC-AGCM-LL_QBOiExp${expno}_${mon5i_s}${year5i}*${ripcode}_$year5i${mon5i}*.nc )
files_miroc_e=( $dir_miroc_e/${var}_${dimcode}_MIROC-ESM_QBOiExp${expno}_${mon5i_s}${year5i}*${ripcode}_$year5i${mon5i}*.nc )
if [ "$expno" = "5" ] ; then
  files_mri=( $dir_mri/mri_run-Fqboi04_QBOhindcast_${exp5i}01${ripcode:0:2}.atm$diminfo$timeint.*.nc )
elif [ "$expno" = "1" ] ; then
  files_mri=( $dir_mri/mri_run-Fqboi04_amip_003.atm$diminfo$timeint.*.nc )
fi
#-------------------------------------------------

st=0
for ff in ${files_cmam[@]} ${files_miroc_a_ll[@]} ${files_miroc_e[@]} ${files_mri[@]} ; do
  if [ ! -e $ff ] ; then
    echo "FILES DO NOT EXIST:  $ff"  ;  st=1
  fi
done
[ $st -eq 0 ] || ( exit $st )

for oo in  \
start_CMAM ${#files_cmam[@]} ${files_cmam[@]}  \
start_MIROC-AGCM-LL ${#files_miroc_a_ll[@]} ${files_miroc_a_ll[@]}  \
start_MIROC-ESM ${#files_miroc_e[@]} ${files_miroc_e[@]}  \
start_MRI ${#files_mri[@]} ${files_mri[@]}  \
; do
  echo $oo
done

