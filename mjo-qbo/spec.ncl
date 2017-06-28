;==============================================================================
; load ncl libraries
;------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

begin
;==============================================================================
; set directory root
;------------------------------------------------------------------------------
;CASE	= "f2000_cntl"
CASE	= "f2000_qbo_uni_test"
MODEL   = "$DATD/mjo-qbo/"
OUTPUT	= "$DATD/mjo-qbo/"+CASE+"/"
FIGDIR	= "./fig/"
SCRIPT	= get_script_prefix_name ()
print ("run "+SCRIPT)
;------------------------------------------------------------------------------
Fout	= OUTPUT+"wmi_index."+CASE+".nc"
Fig	= FIGDIR+SCRIPT+"_"+CASE
;------------------------------------------------------------------------------
yr_sp = 6  ; spin-up years 
nyr = 24
if (CASE .eq. "f2000_cntl") then
  yr_sp = 0
  nyr = 10
end if
;------------------------------------------------------------------------------
; read ucomp
;------------------------------------------------------------------------------
files	= systemfunc ("ls "+MODEL+CASE+"/"+CASE+".var-mjo.remap.daily.nc")
fuo  	= addfile (files, "r")
;------------------------------------------------------------------------------
n1 = 365*yr_sp
n2 = (n1-1) + 365*nyr
u200	= fuo->U200(n1:n2,:,:)
printVarSummary (u200)
u850	= fuo->U850(n1:n2,:,:)
olr	= fuo->FLUT(n1:n2,:,:)
;------------------------------------------------------------------------------
dims	= dimsizes (u200)
ntim	= dims(0)
nlat	= dims(1)
nlon	= dims(2)
;------------------------------------------------------------------------------
; compute daily climatologies and anomalies
;------------------------------------------------------------------------------
time	= u200&time
ymd	= ut_calendar (time, -2)
ymdhms	= ut_calendar (time, 0)
year	= floattointeger (ymdhms(:,0))
month	= floattointeger (ymdhms(:,1))
day	= floattointeger (ymdhms(:,2))
ddd	= day_of_year (year, month, day)
yyyyddd	= year*1000 + ddd
yrfrac	= yyyymmdd_to_yyyyfrac (ymd, 0.)
;------------------------------------------------------------------------------
nhar	= 3
u200ClmDay	= clmDayTLL (u200, yyyyddd)
u200ClmDay_sm	= smthClmDayTLL (u200ClmDay, nhar)
u850ClmDay	= clmDayTLL (u850, yyyyddd)
u850ClmDay_sm	= smthClmDayTLL (u850ClmDay, nhar)
olrClmDay	= clmDayTLL (olr, yyyyddd)
olrClmDay_sm	= smthClmDayTLL (olrClmDay, nhar)
;------------------------------------------------------------------------------
u200Anom	= calcDayAnomTLL (u200, yyyyddd, u200ClmDay_sm)
u850Anom	= calcDayAnomTLL (u850, yyyyddd, u850ClmDay_sm)
olrAnom		= calcDayAnomTLL (olr, yyyyddd, olrClmDay_sm)
printVarSummary (olrAnom)
;------------------------------------------------------------------------------
; take average over latS to latN
;------------------------------------------------------------------------------
latS    = -10
latN    =  10
tmp	= u200Anom(time|:,lon|:,lat|:)
u2002D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
tmp	= u850Anom(time|:,lon|:,lat|:)
u8502D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
tmp	= olrAnom(time|:,lon|:,lat|:)
olr2D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
delete (tmp)

;;
nameseason = (/"winter","summer"/)
;nameseason = (/"annual"/)
optWFP = True
optWFP@cnLevelSelectionMode = "ManualLevels"
optWFP@cnMaxLevelValF = 3.0
optWFP@cnMinLevelValF = 0.25
optWFP@cnLevelSpacingF = 0.25

do ns=0, dimsizes(nameseason)-1
  wf = mjo_wavenum_freq_season (olr2D, ymd,nameseason(ns), False)
  mjo_wavenum_freq_season_plot(wf, nameseason(ns), FIGDIR, "eps", Fig, optWFP) 
end do
;;

end

