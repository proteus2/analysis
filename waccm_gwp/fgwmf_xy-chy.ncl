load "$UTIL/ncl/header.ncl"

xspace = 83;170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 7, 6 ;1, 4

begin

;--------------------------------------- parameters
 hdname  = "test0";_500"
 expname = "0";"500"
 varname = "TAU"+(/"E","W"/)
 titname = (/"East","West"/)+"ward-MF"

 y1 = 1
 y2 = 10
 p = 300  ;5 ;70 ;300
; monselect = 1  ; 0,1,2,3,4 - Annual, DJF, MAM, JJA, SON
 monselect = loopctl

;--------------------------------------- workstation
 if (monselect .eq. 0) then
   mon = ispan(1,12,1)
   tstr = "Annual"
 else
;   monlist = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
;   tstrlist = (/"DJF","MAM","JJA","SON"/)
;   mon = monlist(monselect-1,:)
;   tstr = tstrlist(monselect-1)
   tstrlist = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                "Oct","Nov","Dec"/)
   mon = monselect
   tstr = tstrlist(monselect-1)
 end if

 wks = cr_wks("fig/fgwmf_xy-p"+p+"-"+tstr+"-"+expname,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)

 fname = new(nt,string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(n) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
              sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0),"r")
 x = f0->lon
 y = f0->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 print( f0->ilev({p}) )

 ; get variables
 var = new((/nv,ny,nx/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n),"r")
   tmp = abs( dim_avg_n( f1->$varname(iv)$(:,{p},:,:), 0 ) )
   var(iv,:,:) = var(iv,:,:) + tmp(:,:)
 end do
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)*1.e3

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
; var = where( var .eq. 0., var@_FillValue, var )
 xa = array_append_record(x,(/360.d0/),0)
 vara = new((/nv,ny,nx+1/),float,1.e20)
 vara(:,:,:nx-1) = var(:,:,:)
 vara(:,:,nx   ) = var(:,:,0)

;--------------------------------------- cn level
; cnlev = fspan(-4,4,21)
; cnlev = (/0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56/)
 cnlev = ispan(1,12,1)*0.2
 cnlev(dimsizes(cnlev)-2:) = (/3,4/)
 cnfil = (/5,6,7,8,10,12,14,15,16,17,19,21,0/)

;--------------------------------------- contour
 wi = 60
 he = 60
 xf = 15 + (/0,1,0,1/)*(wi+10)
 yf = 10 + (/0,0,1,1/)*(he+12)

 lbopt = new(nv,string)
 lbopt(:) = ""
 lbopt(0) = "(mPa)"

do iv=0, nv-1

 map = create "mp" mapPlotClass wks
  "mpProjection"          : "Stereographic" ;"AzimuthalEquidistant"
  "mpGridAndLimbOn"       : True
  "mpGridLineColor"       : colind
  "mpGridLineThicknessF"  : 1.2
  "mpGridLineDashPattern" : 2
  "mpEllipticalBoundary"  : True
;  "mpPerimOn" : True

  "mpLimitMode"    : "LatLon"
  "mpCenterLonF"   : 0
  "mpCenterLatF"   : -90
  "mpMinLonF"      : -180
  "mpMaxLonF"      : 180
  "mpMinLatF"      : -90
  "mpMaxLatF"      : -45
  "mpGridLatSpacingF" : 15
  "mpGridLonSpacingF" : 30
  "mpGridPolarLonSpacingF" : 90
 end create
 set_size_mp(map,xf(iv),yf(iv),wi,he)

 contour1 = cnshadeplot(wks,xa,y,vara(iv,:,:),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.02
  "pmLabelBarParallelPosF" : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.64
  "pmLabelBarHeightF"        : 0.08
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,-180,180,-90,-45,False,False)
 axis_str(contour1,titname(iv),"","")

 cnfinalize(contour1,wks)

 overlay(map,contour1)

 draw(map)

end do


 frame(wks)

end

end do  ; nnn

