#! /bin/bash

exit_msg(){
  case $1 in
  1 ) 
    echo "This program requires WGET for downloading data, which are"
    echo "not found. ( See www.gnu.org/software/wget/ )"
    echo "Or, you can download the data manually using a web browser,"
    echo "e.g., Firefox by typing"
    echo "-----------------------------------------------------------"
    echo "firefox $FUB_datURL"
    echo "-----------------------------------------------------------"
    echo "and move the file '$fname' to the directory 'process'."
    echo "Finally, you can re-run this after setting 'opt_update=no'."
    ;;
  2 )
    echo "This program requires NCL for reading and plotting data,"
    echo "which is not found. NCL homepage: http://www.ncl.ucar.edu"
    ;;
  3 )
    echo "The environmental variable 'NCARG_ROOT' is not found. It"
    echo "should be set as the root directory of NCL. You can set it"
    echo "in the script 'run_qbo'."
    ;;
  4 )
    echo "The NCL script is exited with error. Please refer to"
    echo "'process/log_ncl'."
    ;;
  5 )
    echo "The input pressure is not available. The available pressures"
    echo "are listed in 'process/log_ncl'."
    ;;
  * )
    echo "Non-identified error"
  esac
  exit $1
}

