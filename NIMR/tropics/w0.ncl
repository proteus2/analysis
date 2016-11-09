load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"u","rhov","rhodz_dt","theta"/)

 yy = "1983"
 mm = "05"
 z0 = 18.59 ;16.95 ;20.95 ;17.41 ;20.95
 lat0 = 90;25
 n0 = 4  ;23

 wks = cr_wks("fig/w0",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

;--------------------------------------- read data
 nv = dimsizes(varname)
 nw = 4

; fname1 = "/data15/HG2AMIP/L60CGW-t/pa/uanuja.pa_"+yy+mm+"2800.nc"
 fname1 = "/data15/HG2AMIP/L60CGW-t/pa/uanuja.pa_"+yy+"060100.nc"
 f1 = addfile(fname1,"r")
 w0 = 0.5*( f1->dz_dt(n0,{z0*1.e3-600.},{-lat0:lat0},:) + \
            f1->dz_dt(n0,{z0*1.e3+600.},{-lat0:lat0},:) )
 t0 = 0.5*( f1->theta(n0,{z0*1.e3-600.},{-lat0:lat0},:) + \
            f1->theta(n0,{z0*1.e3+600.},{-lat0:lat0},:) )
 t0m = avg(t0)
 t0 = t0m - t0

 lon = f1->longitude
 lat = f1->latitude({-lat0:lat0})
 nx0 = dimsizes(lon)
 ny = dimsizes(lat)

; fw = cfftf(w0,0.,0)/tofloat(nx0)
; fw20 = fw(:,:,:nk*2-1)
; fw20(:,:,3:) = 0.
; fw20(:,:,0:1) = 0.

;--------------------------------------- process 0
 xx = fspan(0.,12.5663706,nx0+1)
; do j=0, ny-1
;   w0(j,:) = 5.e-3*sin(xx(:nx0-1))
; end do
; do j=0, ny-1
; do jj=0, j
;   w0(j,j*2) = 1.e-1
; end do
; end do

 wiy = w0
 wiy(1:ny-2,:) = 0.5*w0(1:ny-2,:) + 0.25*(w0(0:ny-3,:)+w0(2:ny-1,:))
 wix = w0
 wix(:,1:nx0-2) = 0.5*w0(:,1:nx0-2) + 0.25*(w0(:,0:nx0-3)+w0(:,2:nx0-1))
 wix(:,0) = 0.5*w0(:,0) + 0.25*(w0(:,nx0-1)+w0(:,1))
 wix(:,nx0-1) = 0.5*w0(:,nx0-1) + 0.25*(w0(:,nx0-2)+w0(:,0))

 wixy = 0.5*(wix+wiy)

 nw = 4
 vvw = (/ w0, wix, wiy, wixy /)
 vvt = (/ t0, t0, t0, t0 /)

; fw = cfftf(vvw,0.,0)
; fw(:,:,:,:2) = 0.
; fw(:,:,:,4:) = 0.
; fwr = cfftb(fw,1)*2.

; ft = cfftf(vvt,0.,0)
; ft(:,:,:,:2) = 0.
; ft(:,:,:,4:) = 0.
; ftr = cfftb(ft,1)*2.

; vvw = fwr
; vvt = ftr

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)*2.*5.*4.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlevh = fspan(-1.9,1.9,20)*4.*5.

;--------------------------------------- contour
 wi = 32*3
 he = 20*3
 xf = 8 + (/0,1,0,1/)*40
 yf = 10 + (/0,0,1,1/)*30

 mstr  = new(nw,string)
 xstr  = new(nw,string)
 ystr  = new(nw,string)
 lbopt = new(nw,string)

 mstr = (/"w0","w1_x","w1_y","w1_xy"/)
 xstr(:) = ""
 ystr(:) = ""
 lbtxt = (/"","","","","(mm s:S:-1:N:)"/)

do il=0, 0;nw-1

 contour1 = cnshadeplot(wks,lon,lat,vvw(il,:,:)*1.e3,cnlev,cnfil,lbtxt(il))
; cnaddline(contour1,wks,lon,lat,vvt(il,:,:),cnlevh,0,"D",0,"T",0.5,"C",(/colindr,-1,colindb/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.15
;  "pmLabelBarParallelPosF"   : 0.50-0.40
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.0
;  "pmLabelBarHeightF"        : 0.15
;  "lbLabelFontHeightF" : 0.025
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
;  "lbTitleFontHeightF" : 0.030
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,360,-90,90,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))

 cnfinalize(contour1,wks)

; vf = create "epf_vec" vectorFieldClass defaultapp
;  "vfDataArray" : vv2(il,0:1,2::3,:)
;  "vfXArray"    : lon
;  "vfYArray"    : lat(2::3)
; end create
;
; vec1 = create "vector" vectorPlotClass wks
;  "vcVectorFieldData" : vf
;;  "vcRefMagnitudeF" : ref
;  "vcRefLengthF"    : 0.05
;  "vcRefAnnoOn"     : False
;;  "vcLineArrowThicknessF" :  thc(nn)
;;  "vcLineArrowColor"      :  col(nn)
;;  "vcRefAnnoString1" : ref*tmpc(nn)
;  "vcRefAnnoString2On" : False
;;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
; end create
; set_size(vec1,xf(il),yf(il),wi,he)
; axis_range(vec1,-180/kc(0),180/kc(0),-20,20,False,False)
; axis_str(vec1,"","","")
;
; overlay(contour1,vec1)
 draw(contour1)
; NhlRemoveOverlay(contour1,vec1,False)

end do

 frame(wks)

end



