###############################################################################
# Offline wrapper PROGRAM and support modules to populate the runtime variables
# namelist and offline input fields required to run UM gw_ussp_core.
#
# Created by Andrew Bushell - April 2017.
# 
# Current code owner(s) of ~/utils/off_gw_ussp - andrewbushell
#
###############################################################################
#
# src/ contains F90 code for offline testing of UM subroutine gw_ussp_core.
#
# To compile gw_ussp_core_mod [UM vn10.7 r34978] in offline mode:
# 
ifort -c gw_ussp_core_mod.F90
ifort -c gw_ussp_offline_mod.F90
ifort -c runvnamelist_off_mod.F90

ifort runvnamelist_off_mod.o gw_ussp_offline_mod.o gw_ussp_core_mod.o gw_ussp_core_offline.F90 -o gw_ussp_offline.exe

#
# data/ contains example files from QBOi project core paper
#
### runvar_QBOiP0_GA7 specifies runtime variables to gw_ussp_offline.exe
#    Note: adjust maximum input file dimensions in runvnamelist_off_mod
#          choose ussp_launch_factor = 1.3, cgw_scale_factor = 0.9600
#
### QBOiP0inp_1993-05-01.txt are pre-constructed input files based on ERA-I data
### QBOiP0inp_2005-11-01.txt
#    Advice on construction of input files (IDL script) can be offered
#
### fptot_UMGA7_P0_1993-05-01.txt are example output files
### fptot_UMGA7_P0_2005-11-01.txt
#
# To repeat example run:
#
ln -f -s runvar_QBOiP0_GA7 fort.16
ln -f -s QBOiP0inp_2005-11-01.txt fort.17

./gw_ussp_offline.exe > fptot_UMGA7_P0_2005-11-01.txt
mv fort.25 fptot_UMGA7_P0_2005-11-01.out

