;==============================================================================
; load ncl libraries
;------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
load "MJO_composite.ncl"

begin
;==============================================================================
; set directory root
;------------------------------------------------------------------------------
;CASE	= "f2000_cntl"
CASE	= "f2000_qbo_uni_test"
MODEL   = "$DATD/mjo-qbo/"
OUTPUT  = "$DATD/mjo-qbo/"+CASE+"/"
FIGDIR	= "./fig/"
SCRIPT	= get_script_prefix_name ()
print ("run "+SCRIPT)
;------------------------------------------------------------------------------
READ	= False 
Fout	= OUTPUT+SCRIPT+"."+CASE+".nc"
Fig     = FIGDIR+SCRIPT+"_"+CASE
;------------------------------------------------------------------------------
LagS	= -25
LagE	= 25
nLag	= 51
nPhase	= 8
;---  The below should be the same as in indx_wmi.ncl  ------------------------
yr_sp = 6  ; spin-up years 
nyr = 24
if (CASE .eq. "f2000_cntl") then
  yr_sp = 0
  nyr = 10
end if
;------------------------------------------------------------------------------
; read data
;------------------------------------------------------------------------------
if (.not.READ) then
print ("read the rmm indices")
;findx	= addfile (MODEL+"rmm_index.nc", "r")
findx	= addfile (OUTPUT+"wmi_index."+CASE+".nc", "r")
indx	= findx->indx
rmm1	= findx->rmm1
rmm2	= findx->rmm2
phase	= findx->phase
;------------------------------------------------------------------------------
files   = systemfunc ("ls "+MODEL+CASE+"/"+CASE+".var-mjo.remap.daily.nc")
fin     = addfile (files, "r")
;------------------------------------------------------------------------------
print ("read the data")
n1 = 365*yr_sp
n2 = (n1-1) + 365*nyr
olr     = fin->FLUT(n1:n2,{-30:30},:)
printVarSummary (olr)
;------------------------------------------------------------------------------
dims    = dimsizes (olr)
ntim    = dims(0)
nlat    = dims(1)
nlon    = dims(2)
;------------------------------------------------------------------------------
time	= olr&time
ymd	= ut_calendar (time, -2)
ymdhms  = ut_calendar (time, 0)
year    = floattointeger (ymdhms(:,0))
month   = floattointeger (ymdhms(:,1))
day     = floattointeger (ymdhms(:,2))
ddd	= day_of_year (year, month, day)
yyyyddd	= year*1000 + ddd
print(ddd)
;------------------------------------------------------------------------------
; calculate anomaly and do filtering
;------------------------------------------------------------------------------
print ("extract the intraseasonal variability")
nhar    = 3
olrClmDay	= clmDayTLL (olr, yyyyddd)
olrClmDay_sm	= smthClmDayTLL (olrClmDay, nhar)
olr		= calcDayAnomTLL (olr, yyyyddd, olrClmDay_sm)
;------------------------------------------------------------------------------
ihp     = 2                             ; bpf=>band pass filter
nWgt    = 201
sigma   = 1.0                           ; Lanczos sigma
fca     = 1./100.
fcb     = 1./20.
wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
olr	= wgt_runave_n_Wrap (olr, wgt, 0, 0)
;------------------------------------------------------------------------------
iLH			= LH_mjo_composite (indx, phase)
iLH(0:100)		= 0
iLH(ntim-101:ntim-1)	= 0
;------------------------------------------------------------------------------
print ("do lag composite for "+nLag+" days")
colr	= new ((/nPhase,nLag,nlat,nlon/), typeof(olr), getFillValue(olr))
colr	= 0.
nDays	= new (nPhase, "integer")
do np=1, nPhase
  print ("phase "+np)
  iMJO	= ind ((month.le.2.or.month.ge.12).and.iLH.eq.1.and.phase.eq.np)
  if (.not.all(ismissing(iMJO))) then
    do nl=LagS, LagE
      colr(np-1,nl-LagS,:,:) = dim_avg_n_Wrap (olr(iMJO+nl,:,:), 0)
    end do
  end if
  nDays(np-1) = dimsizes (iMJO)
  delete (iMJO)
end do
colr!0		= "phase"
colr&phase	= ispan (1,nPhase,1)
colr!1		= "lag"
colr&lag	= ispan (LagS,LagE,1)
;==============================================================================
system ("rm "+Fout)
nCDF	= addfile (Fout, "c")
nCDF->colr=colr
nCDF->nDays=nDays
end if
;------------------------------------------------------------------------------
if (READ) then
nCDF	= addfile (Fout, "r")
colr	= nCDF->colr
nDays	= nCDF->nDays
end if
printVarSummary (colr)
colr&lat@units	= "degrees_north"
colr&lon@units	= "degrees_east"
;==============================================================================
; plot
;------------------------------------------------------------------------------
print ("plot the data")
wks = gsn_open_wks ("eps", Fig)
gsn_define_colormap (wks, "ViBlGrWhYeOrRe")
plot= new (8, graphic)
;==============================================================================
res				= True           
res@gsnDraw			= False
res@gsnFrame			= False
res@mpFillOn			= False
res@mpMinLatF			= -30.
res@mpMaxLatF			= 30.
res@mpCenterLonF		= 210.
res@gsnMajorLatSpacing		= 20.
res@gsnMajorLonSpacing		= 60.
res@gsnSpreadColors		= True
res@cnFillOn			= True
res@cnLinesOn			= False
res@cnLineLabelsOn		= False
res@lbLabelBarOn		= False
res@gsnContourZeroLineThicknessF= 2
res@gsnContourNegLineDashPattern= 1
res@cnLevelSelectionMode	= "ManualLevels"
res@cnLevelSpacingF		= 2.
res@cnMaxLevelValF		= 10.*res@cnLevelSpacingF
res@cnMinLevelValF		= -1.*res@cnMaxLevelValF
;------------------------------------------------------------------------------
resP				= True
resP@gsnFrame			= False
resP@gsnMaximize		= True
resP@gsnPanelLabelBar		= True
resP@gsnPaperOrientation	= "portrait"
;------------------------------------------------------------------------------
;res@gsnLeftString	= "Z at 300 mb"
;res@gsnRightString	= "m"
;iLag	= ispan (-15, 15, 5)
;resP@gsnPanelFigureStrings = "Lag "+iLag
;resP@gsnPanelLeft	= 0.
;resP@gsnPanelRight	= 0.5
;do nl=0, dimsizes (iLag)-1 
;  plot(nl)	= gsn_csm_contour_map_ce (wks,colr(2,{iLag(nl)},:,:),res)
;end do
;gsn_panel (wks, plot, (/7,1/), resP)
;------------------------------------------------------------------------------
;resP@gsnPanelLeft	= 0.5
;resP@gsnPanelRight	= 1.
;do nl=0, dimsizes (iLag)-1 
;  plot(nl)	= gsn_csm_contour_map_ce (wks,colr(6,{iLag(nl)},:,:),res)
;end do
;gsn_panel (wks, plot, (/7,1/), resP)
;------------------------------------------------------------------------------
res@gsnLeftString	= ""
res@gsnRightString	= ""
do np=0, nPhase-1
  plot(np)	= gsn_csm_contour_map_ce (wks,colr(np,{0},:,:),res)
end do
;resP@lbAutoManage = True
;resP@lbOrientation = "Vertical"
;resP@lbTitleString = "[W m:S:-2:N:]"
gsn_panel (wks, plot, (/nPhase,1/), resP)
frame (wks)
;==============================================================================
end
