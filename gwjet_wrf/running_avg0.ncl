setvalues NhlGetWorkspaceObjectId()
 "wsMaximumSize" : 256000000
end setvalues

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************
begin
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("eps","fig/running_avg0")         ; open a ncgm file
;  gsn_define_colormap(wks,"StepSeq25")
  gsn_define_colormap(wks,"cmp_b2r")

  res                   = True
  res@mpProjection      = "Orthographic"       ; choose projection
;  res@mpProjection      = "Satellite"       ; choose projection
;  res@mpSatelliteDistF  = 1.5
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
  res@mpCenterLonF      = 0.
  res@mpCenterLatF      = 30.
  res@mpLimitMode       = "Angles"
  res@mpTopAngleF       = 60
  res@mpBottomAngleF    = 0
  res@mpLeftAngleF      = 0
  res@mpRightAngleF     = 60

  res@cnFillOn          = True              ; color plot desired
  res@cnLineLabelsOn    = False
  res@cnLinesOn         = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels     = fspan(-10,10,21)*3.
;  res@cnFillColors = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
  res@cnFillColors = ispan(2,65,3)
  res@lbLabelBarOn = False

;  res@gsnSpreadColors     = True            ; use full range of color map
;  res@gsnSpreadColorStart = 2 
;  res@gsnSpreadColorEnd   = -3 
;  res@lbLabelFontHeightF  = 0.015           ; label bar font height
 
  res@mpGridAndLimbOn = True
;  res@mpGridLatSpacingF = 0.1*10.
;  res@mpGridLonSpacingF = 0.1*10.
  res@mpGridLineThicknessF = 0.1

  res@mpOutlineBoundarySets = "NoBoundaries"

  fdat = "/data18/GW_jet_wrf/Hist/auxhist_files/auxhist_d01_2000-01-10_12:00:00.nc"
  f1 = addfile(fdat,"r")
  v = f1->W_300(144,972:,:)*1.e3
  v@_FillValue = 1.e20
  v(:,1458:) = v@_FillValue
;  v = (/ (/-1.,-1./), (/-1.,-1./), (/-1.,-1./) /)
;  v@_FillValue = -1.
  v!0 = "lat"
  v!1 = "lon"
  v&lat = fspan(0,90,972)
  v&lon = fspan(0,360.*(3887./3888.),3888)
;  v&lat = fspan(-90,90,3)
;  v&lon = fspan(0,180,2)
  v&lat@units = "degrees_north"
  v&lon@units = "degrees_east"

  plot = gsn_csm_contour_map(wks,v,res)  ; create the plot
;  plot = gsn_csm_contour_map(wks,v(::3,::3),res)  ; create the plot

end

