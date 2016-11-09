
rm -f log_ncl *.eps

## define functions

. ./ftn.sh

## download data

exe_wget="`command -v wget`"  ;  st_wget=$?
if [ $st_wget -ne 0 ] ; then exit_msg 1 ; fi
echo -e "\n Downloading $yyyy ...\n"
export year1 year2 yearf
./get-merra

## read data and draw plots using NCL

echo -e "\n Plotting with NCL...\n"
exe_ncl="`command -v ncl`"
if [ $? -ne 0 ] ; then exit_msg 2 ; fi
exe_ncl="`dirname $exe_ncl`"  ;  exe_ncl="`dirname $exe_ncl`"
export NCARG_ROOT=${NCARG_ROOT:-"$exe_ncl"}
if [ ! -d $NCARG_ROOT/lib/ncarg/nclscripts ] ; then exit_msg 3 ; fi
out1="easmi_${year1}-${year2}_${yearf}"
outfiles=( "$out1" )
export year1 year2 yearf
export out1
ncl fig_easmi.ncl > log_ncl
st_ncl=$?
if [ $st_ncl -ne 10 ] ; then exit_msg 4 ; fi

echo -e "\n Output files...\n"
if [ "$filetype" == "" ] ; then filetype="eps" ; fi
if [ "$filetype" != "eps" ] ; then
  exe_im="`command -v convert`"
  st_im=$?
  for o in ${outfiles[@]} ; do
    if [ $st_im -eq 0 ] ; then
      convert -density 300 $o.eps $o.$filetype
    else
      if   [ "$filetype" == "png" ] ; then dev_gs="png16m"
      elif [ "$filetype" == "jpg" ] ; then dev_gs="jpeg"
      fi
      gs -q -sDEVICE=$dev_gs -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=$o.$filetype $o.eps
    fi
    mv -f $o.$filetype $o.eps ../ ; echo "        $o.$filetype"
  done
else
  mv -f *.eps ../
  for o in ${outfiles[@]} ; do echo "        $o.eps" ; done
fi
echo

exit 0

