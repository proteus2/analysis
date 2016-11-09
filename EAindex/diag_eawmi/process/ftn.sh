#! /bin/bash

exit_msg(){
  case $1 in
  1 ) 
    echo "This program requires WGET for downloading data, which are"
    echo "not found. ( See www.gnu.org/software/wget/ )"
    ;;
  2 )
    echo "This program requires NCL for reading and plotting data,"
    echo "which is not found. NCL homepage: http://www.ncl.ucar.edu"
    ;;
  3 )
    echo "The environmental variable 'NCARG_ROOT' is not found. It"
    echo "should be set as the root directory of NCL. You can set it"
    echo "in the script 'run_wtem'."
    ;;
  4 )
    echo "The NCL script is exited with error. Please refer to"
    echo "'process/log_ncl'."
    ;;
  * )
    echo "Non-identified error"
  esac
  exit $1
}

