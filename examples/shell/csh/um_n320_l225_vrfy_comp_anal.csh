#!/bin/csh -xvf
###################################################################
# This verification is based on the anal-time (or, observation), so 
#      verification time 24 means that performed between ANLTIM anal 
#      and P1DTIM fcst.
#                                                   -- by jhson 
###################################################################
### DEFINE DIRECTORY
   source ~/UTIL/ENVI/util_src
   source ~/UTIL/ENVI/vrfy_ukcold_src

   set ANLTIM = `cat $UMGLLOGO/vrfy.tim`
   set P1DTIM = `$UTILEXET/maketim  $ANLTIM -24`
   set P2DTIM = `$UTILEXET/maketim  $P1DTIM -24`
   set P3DTIM = `$UTILEXET/maketim  $P2DTIM -24`
   set P4DTIM = `$UTILEXET/maketim  $P3DTIM -24`
   set P5DTIM = `$UTILEXET/maketim  $P4DTIM -24`

   set LOGFILE  = $UMGLLOGO/um_n320_vrfy_all_log.$ANLTIM

#############################################################################
### MAKE INFORMATION INPUT FILE (input files)

    echo $ANLTIM >  $UMGLLOGO/um_n320_l225_vrfy_infiles.dat
    foreach tim   ( $ANLTIM $P1DTIM $P2DTIM $P3DTIM $P4DTIM $P5DTIM )
       if ( $tim == $ANLTIM ) then
            set day = 0
       else if ( $tim == $P1DTIM ) then
            set day = 1
       else if ( $tim == $P2DTIM ) then
            set day = 2
       else if ( $tim == $P3DTIM ) then
            set day = 3
       else if ( $tim == $P4DTIM ) then
            set day = 4
       else if ( $tim == $P5DTIM ) then
            set day = 5
       endif
    foreach vari  ( h t u v )
    foreach level ( 850 500 200 )
       if ( ! -e  $UMGLDAIO/um_glob_out_bin_${vari}${level}.${tim} ||                    \
            ! -e  $UMGLDAIO/um_glob_out_bin_${vari}${level}.${tim} || ) then
            echo  '0' $day $vari $level $UMGLDAIO/um_glob_out_bin_${vari}${level}.${tim}  \
                      >> $UMGLLOGO/um_n320_l225_vrfy_infiles.dat
       else
            echo  '1' $day $vari $level $UMGLDAIO/um_glob_out_bin_${vari}${level}.${tim}  \
                      >> $UMGLLOGO/um_n320_l225_vrfy_infiles.dat
       endif
    end
    end
    end

#############################################################################
### DEFINE INPUT/OUTPUUT FILES
## input file
    setenv  F_FF11  $UMGLLOGO/um_n320_l225_vrfy_infiles.dat
#
## output file
    setenv  F_FF99  $UMGLDAOU/um_n320_l225_vrfy_anal.$ANLTIM
# 
#############################################################################
### set SOURCE

    set  SOURCE = um_n320_l225_vrfy_comp_anal

    if ( ! -e $UMGLEXET/$SOURCE.e ) then
        cd $UMGL_SRC
        ifort -static $UMGL_SRC/$SOURCE.f90 -o $UMGLEXET/$SOURCE.e
    endif

### run program
  ($UMGLEXET/$SOURCE.e > $UMGLLOGO/$SOURCE.lis ) >& $UMGLLOGO/$SOURCE.log

#############################################################################
###ERROR CHECK
   set     PCHK=$status
   echo    $PCHK : $SOURCE
   set     LOGDATE="DATE:"`date +%Y/%m/%d/%H:%M:%S`
   if (  "$PCHK" !=  0 )  then
   echo     $LOGDATE ' ERROR FOUNDED          '$SOURCE            >> $LOGFILE
   echo     '<<<WARNING>>>  CHECK '$VKWFLOGO'/'$SOURCE'.log   '   >> $LOGFILE
   exit     1
   endif
   echo     $LOGDATE ' SUCESSFUL FINISHED     '$SOURCE            >> $LOGFILE
#############################################################################
