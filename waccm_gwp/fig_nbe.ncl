load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; Load the NCL file that
                     ; contains the gsn_* functions used below.
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
                     ; contains the gsn_* functions used below.  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;----------------------------------------------------------------------------------
begin

  imon = 0

  smon = (/"DJF","MAM","JJA","SON"/)
  smons= (/"djf","mam","jja","son"/)

  f = addfile("dat/dnbe_xyz."+smons(imon)+".0001-0005.exp2.nc","r")   ; file addition
  ftg_climatology = f->D_NBE(0,3:,:,:)
  ftg_climatology = (/ftg_climatology/)*1.e9

; ----------------------------------------------------------------------------------
  wks = gsn_open_wks("eps","nbe_climatology_"+smon(imon)+"_revision_370")  ; Open an X11 workstation. 
  plot_1 = new(2,graphic); icreate graphic array
;  plot_2 = new(2,graphic); icreate graphic array
  gsn_define_colormap(wks,"prcp_3")  ; select  colormap 
  res                       = True             ; plot mods desired
  res@gsnDraw               = False            ; don't draw yet
  res@gsnFrame              = False            ; don't advance frame yet
  
;  res@tiMainString          = "MERRA Jan. 1980 ~ 2011 Average"  ; title
  res@gsnLeftString         = ""
  res@gsnCenterString       = ""
  res@gsnRightString        = smon(imon)
;  res@gsnCenterStringFontHeightF  = 0.02
  res@gsnRightStringFontHeightF  = 0.04
  res@tmXBLabelFontHeightF  = 0.025
  res@tmYLLabelFontHeightF  = 0.025
;  res@cnLevelSpacingF       = 15.              ; contour interval
;  res@cnLineLabelsOn        = True             ; turn on line label  
   res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels              = ispan(0,11,1)*0.3
  res@cnFillOn              = True
  res@cnFillColors =  (/0,5,6,7,8,10,12,14,15,16,17,19,21/)
;  res@gsnSpreadColors       = True             ; use full colormap
;  res@cnLevelSpacingF       = 0.005
  res@lbLabelFontHeightF =.009                 ; make labels larger
   res@lbLabelBarOn          = True
   res@cnLineLabelsOn        = False             ; turn on line label  
   res@cnLinesOn             = False
  res@lbLabelFontThicknessF  = 2
  res@lbLabelFontHeightF  = 0.016
  res@lbTitleString    = "(10:S:-9:N: s:S:-2:N:)"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF  = 0.020
  res@lbTopMarginF  = 0.4
  res@vpOn             = True
  res@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
  res@vpYF             = 0.8                 ; default is 0.8
  res@vpWidthF         = 1.         ;Default: 0.6
  res@vpHeightF        = 1.             ;Default: 0.6 
  res@tmYLBorderOn   =  True
   res@mpMinLatF  =  -75.0
   res@mpMaxLatF  =  75.0
   res@mpFillOn              = False 
   res@tmXBMode              = "Explicit"
  res@tmXBValues            = (/-180,-90,0,90,180/) ; y axis value i want
  res@tmXBLabels            = (/"180W","90W","0","90E","180E"/)
;************************************************
; use mask function to mask out a range of the data 
;************************************************
  b  = ftg_climatology(0,:,:)                               ; trick to keep cv's and atts
;  b  = mask(ftg_climatology(0,:,:),(ftg_climatology(0,:,:).ge. -2000.),True)
   ; shown region : caution
; --------------------------------------------------------------------------------
 b&lon@units = "degrees_east"
 b&lat@units = "degrees_north" 
 plot_1(0) = gsn_csm_contour_map_ce(wks,b,res)     ; draw global map 
  
   res3                       = True             ; plot mods desired
  res3@gsnDraw               = False            ; don't draw yet
  res3@gsnFrame              = False            ; don't advance frame yet
  res3@tmXBLabelFontHeightF  = 0.025
  res3@tmYLLabelFontHeightF  = 0.025
;  res@tiMainString          = "MERRA Jan. 1980 ~ 2011 Average"  ; title
  res3@gsnLeftString         = ""
  res3@gsnCenterString       = ""
  res3@gsnRightString        = smon(imon)
;   res3@gsnCenterStringFontHeightF  = 0.02
  res3@gsnRightStringFontHeightF  = 0.04
;  res@cnLevelSpacingF       = 15.              ; contour interval
;  res3@cnLineLabelsOn        = True             ; turn on line label  
;   res3@cnLevelSelectionMode  = "ExplicitLevels"
  res3@cnLevels              = ispan(0,11,1)*0.3
  res3@cnFillOn              = True             ; turn on color
  res@cnFillColors =  (/0,5,6,7,8,10,12,14,15,16,17,19,21/)
  res3@lbLabelFontHeightF =.009                 ; make labels larger
   res3@cnLineLabelsOn        = False             ; turn on line label  
   res3@cnLinesOn             = True
 res3@lbLabelOn    = False
;  res3@vpOn             = True
  res3@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
  res3@vpYF             = 0.8                 ; default is 0.8
  res3@vpWidthF         = 1.         ;Default: 0.6
  res3@vpHeightF        = 1.             ;Default: 0.6 
   res3@tmYLBorderOn   =  True
   res3@mpMinLatF  =  -75.0
   res3@mpMaxLatF  =  75.0
   res3@mpFillOn              = False 
   res3@lbLabelBarOn       = True
   res3@lbLabelsOn         = True
   res3@lbTitleOn          = True
   res3@lbTitlePosition    = "Bottom"
   res3@lbTitleString      =" "
   res3@gsnContourZeroLineThicknessF = 3.
  res3@lbTopMarginF  =0.4

   d = ftg_climatology(0,:,:)
   d&lon@units = "degrees_east"
   d&lat@units = "degrees_north"
 res3@tmXBMode              = "Explicit"
  res3@tmXBValues            = (/-180,-90,0,90,180/) ; y axis value i want
 res3@tmXBLabels            = (/"180W","90W","0","90E","180E"/) 

; plot_2(0)  = gsn_csm_contour_map_ce(wks,d,res3); jet in 200hPa
;   plot_2(0)  = ShadeLtContour(plot_2(0),-2000,17)      ; shade contours > 30.
   res@lbLabelBarOn       = True
  plot_1(0) = gsn_csm_contour_map_ce(wks,b,res)     ; draw global map       
  
 
res2                       = True             ; plot mods desired
  res2@gsnDraw               = False            ; don't draw yet
  res2@gsnFrame              = False            ; don't advance frame yet

;  res2@tiMainString          = "MERRA Jan. 1980 ~ 2011 Average"  ; title
  res2@gsnLeftString         = "600 hPa"
  res2@gsnCenterString       = ""
  res2@gsnRightString        = smon(imon)
  res2@gsnRightStringFontHeightF  = 0.04
  res2@gsnLeftStringFontHeightF  = 0.04
;  res2@cnLevelSpacingF       = 15.              ; contour interval
;  res2@cnLineLabelsOn        = True             ; turn on line label  
   res2@cnLevelSelectionMode  = "ExplicitLevels"
   res2@cnLevels              = (/0.045,0.05,0.055,0.06,0.07,\
                               0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2/)  
res2@cnFillOn              = True             ; turn on color
  res2@tmXBLabelFontHeightF  = 0.025
  res2@tmYLLabelFontHeightF  = 0.025 
 res2@gsnSpreadColors       = True             ; use full colormap
;  res2@cnLevelSpacingF       = 0.005
  res2@gsnSpreadColorStart  = 72 ;15                    ; start at color 10
  res2@gsnSpreadColorEnd    = 200 ;90  
  res2@lbLabelFontHeightF =.006                 ; make labels larger
   res2@lbLabelFontHeightF =.009                  ; make labels larger
  res2@lbLabelFontThicknessF  = 2
  res2@lbTopMarginF  =0.4
   res2@lbLabelBarOn          = True
   res2@cnLineLabelsOn        = False             ; turn on line label  
   res2@cnLinesOn             = False
  res2@tmXBMode              = "Explicit"
  res2@tmXBValues            = (/-180,-90,0,90,180/) ; y axis value i want
  res2@tmXBLabels            = (/"180W","90W","0","90E","180E"/)
  res2@vpOn             = True
  res2@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
  res2@vpYF             = 0.8                 ; default is 0.8
  res2@vpWidthF         = 1.         ;Default: 0.6
  res2@vpHeightF        = 1.             ;Default: 0.6 
  res2@tmYLBorderOn   =  True
   res2@mpMinLatF  =  -75.0
   res2@mpMaxLatF  =  75.0
   res2@mpFillOn              = False
  
   NhlSetColor(wks,72,1,1,1)

;  gsn_define_colormap(wks,"wind_17lev")  ; select  colormap  
  res4                     = True             ; plot mods desired
  res4@gsnDraw              = False           ; Do not draw plot
  res4@gsnFrame             = False  
  res4@cnLevelSelectionMode = "ExplicitLevels"  ; set manual contour levels
  res4@cnLevels              = (/25,35,45,55,65,75/)
  res4@cnMinLevelValF       =  25.             ; set min contour level
  res4@cnMaxLevelValF       =  80.             ; set max contour level
  res4@cnLevelSpacingF      =  10.            ; set contour spacing
  res4@cnLineThicknessF     =  2.3 
  
  res4@cnFillOn              = False
  res4@cnLinesOn             = True
  res4@gsnSpreadColors       = True             ; use full colormap
;  res4@cnFillOpacityF    = 0.8     ; Half transparent
  res4@tmXBBorderOn      = True
  res4@tmXBOn            = True
  res4@tmXTBorderOn      = True
  res4@tmXTOn            = True
  res4@tmYLBorderOn      = True
  res4@tmYLOn            = True
  res4@tmYRBorderOn      = True
  res4@tmYROn            = True
  res4@vpOn             = True
  res4@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
  res4@vpYF             = 0.8                 ; default is 0.8
  res4@vpWidthF         = 1.          ;Default: 0.6
  res4@vpHeightF        = 1.             ;Default: 0.6 
  res4@lbLabelBarOn          = True
  res4@mpFillOn              = False       ; turn off gray fill
  res4@gsnLeftString         = ""
  res4@gsnCenterString       = ""
  res4@gsnRightString        = ""
;  res4@gsnRightStringFontHeightF  = 10
  res4@tmXBLabelFontHeightF  = 0.025
  res4@tmYLLabelFontHeightF  = 0.025
  res4@lbLabelFontHeightF =.02   ;Default: 0.33 
  res4@cnConstFLabelFontThicknessF = 10
  res4@mpMinLatF  =  -75.0
  res4@mpMaxLatF  =  75.0
     res4@lbLabelFontHeightF =.009                  ; make labels larger
  res4@lbLabelFontThicknessF  = 2
  res4@lbTopMarginF  = 0.4
  res4@mpOutlineOn = False
  res4@tmXBMode              = "Explicit"
  res4@tmXBValues            = (/-180,-90,0,90,180/) ; y axis value i want
 res4@tmXBLabels            = (/"180W","90W","0","90E","180E"/)  
   ; shown region : caution 
;   a&lon@units = "degrees_east"
;   a&lat@units = "degrees_north"
;   plot_1(1) = gsn_csm_contour_map_ce(wks,a, res4)
;plot_1(1) = gsn_csm_contour_map_ce(wks,a, res4)
 ;************************************************
; use mask function to mask out a range of the data 
;************************************************
  c  = ftg_climatology(0,:,:)                               ; trick to keep cv's and atts
  c&lon@units = "degrees_east"
  c&lat@units = "degrees_north"
;  c  = mask(ftg_climatology(0,:,:),(ftg_climatology(0,:,:).ge. 0.044),True)
   ; shown region : caution
; --------------------------------------------------------------------------------  

; plot_2(1) = gsn_csm_contour_map_ce(wks,c, res2)
  
;************************************************
; create panel
;************************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelLabelBar = False                   ; add common colorbar
;   resP@vpOn             = True
;  resP@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
;  resP@vpYF             = 0.8                 ; default is 0.8
;  resP@vpWidthF         = 1.0         ;Default: 0.6
;  resP@vpHeightF        = 5.0             ;Default: 0.6 
 
;  resP@txFont            = 21
  resP@gsnPanelBottom   = 0.05                   ; add space at bottom
;  resP@gsnPanelFigureStrings= (/"1979","1980","1981","1982","1983","1984","1985","1986","1987"/) ; add strings to panel
  resP@amJust   = "TopCenter"
  resP@gsnPanelYWhiteSpacePercent = 0.5
  resP@gsnPanelXWhiteSpacePercent = 1
  ;-----------------
  resP1                  = True                   ; modify the panel plot
  resP1@gsnFrame         = False                  ; don't advance panel plot
  resP1@gsnPanelLabelBar = False                   ; add common colorbar
;   resP1@vpOn             = True
;  resP1@vpXF             = 0.2             ; default is 0.2 (aspect ratio)
;  resP1@vpYF             = 0.8                 ; default is 0.8
;  resP1@vpWidthF         = 1.0         ;Default: 0.6
;  resP1@vpHeightF        = 5.0             ;Default: 0.6 

;  resP1@txFont            = 21
  resP1@gsnPanelBottom   = 0.05                   ; add space at bottom
;  resP1@gsnPanelFigureStrings= (/"1979","1980","1981","1982","1983","1984","1985","1986","1987"/) ; add strings to panel
  resP1@amJust   = "TopCenter"
  resP1@gsnPanelYWhiteSpacePercent = 0.5
  resP1@gsnPanelXWhiteSpacePercent = 1
;  gsn_panel(wks,plot_2,(/2,1/),resP)               ; now draw as one plot
  gsn_panel(wks,plot_1,(/2,1/),resP)
 frame(wks)
end
