load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 3.6

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 var2 = "U"
 lev = 300
 n1s = ispan(719,1151,6)
 var1name = "prt_"+var1+"_"+lev
 var2name = "prt_"+var2+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
; wks = cr_wks("fig/prt16_xy-"+var1name+"-"+var2name,600,0)
 wks = cr_wks("fig/prt-all_xy-"+var1name+"-box.ps",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname1b= fdir+"/x1/prt_d300km/"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname2 = fdir+"/x1/prt_d300km/"+var2name+"__x1__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 f1b= addfile(fname1b,"r")
 f2 = addfile(fname2,"r")
 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

nt = dimsizes(n1s)

do il=0, nt-1

 dt = 5
 nmin = 60/dt
 tday  = (n1s(il)+1)/(nmin*24)
 tmin  = (n1s(il)+1) - tday*(nmin*24)
 thour = tmin/nmin
 tmin  = (tmin - thour*nmin)*dt
 tday  = 4 + tday

 v1 = f1->$var1name$(n1s(il),:,:)
 v1b= f1b->$var1name$(n1s(il),:,:nx/2-1)
; v2 = f2->$var2name$(n1s(il),:,:)

; fname0 = "/data18/GW_jet_wrf/dat/x1/"+var1+"/"+var1+"_"+lev+"__x1__08_00_05__12_00_00.nc"
; f0 = addfile(fname0,"r")
; v0 = f0->$(var1+"_"+lev)$(n1s(il),53:52+ny,:)
; v0 = v0 - v1

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

 xx = ispan(0,nx+nx/2-1,1)*dx
 nxx = dimsizes(xx)

 vv1 = new((/ny,nxx/),float)
 vv1(:,:nx-1) = v1 (:,:)
 vv1(:,nx:  ) = v1b(:,:)

 y!0 = "lat"
 y&lat = y
 xx!0 = "lon"
 xx&lon = xx

 vv1!0 = "lat"
 vv1!1 = "lon"
 vv1&lat = y
 vv1&lon = xx

 pi = 3.141592
;vv1 = 1./30.;yh

; W2a
 dspeed = 0.54
 x10 = -3 + il*dspeed
 x20 = 10 + il*dspeed
 y10 = 52
 y20 = 67
 x1b = 3.
 x2b = 3.
 y1b = 3.
 y2b = 3.
 xb = x1b + x2b
 yb = y1b + y2b
 x1 = x10 + x1b
 x2 = x20 - x2b
 y1 = y10 + y1b
 y2 = y20 - y2b
 if ( x10
 x10 = max((/xx(0), x10/))
 x20 = min((/xx(nxx-1), x20/))
 x1 = max((/xx(0), x1/))
 x2 = min((/xx(nxx-1), x2/))
;print("AAA   "+x10+"  "+x1+"  "+x2+"  "+x20)
 
 ws2a = vv1
 ws2a({y20:},:) = 0.
 ws2a({:y10},:) = 0.
 ws2a(:,{x20:}) = 0.
 ws2a(:,{:x10}) = 0.
 ws2a({y2:y20},{x10:x20}) = ws2a({y2:y20},{x10:x20})* \
     conform(ws2a({y2:y20},{x10:x20}),abs(cos((y({y2:y20})-y({y2}))/yb*pi)),0)
 ws2a({y10:y1},{x10:x20}) = ws2a({y10:y1},{x10:x20})* \
     conform(ws2a({y10:y1},{x10:x20}),abs(cos((y({y10:y1})-y({y1}))/yb*pi)),0)
 ws2a({y10:y20},{x2:x20}) = ws2a({y10:y20},{x2:x20})* \
     conform(ws2a({y10:y20},{x2:x20}),abs(cos((xx({x2:x20})-xx({x2}))/xb*pi)),1)
 ws2a({y10:y20},{x10:x1}) = ws2a({y10:y20},{x10:x1})* \
     conform(ws2a({y10:y20},{x10:x1}),abs(cos((xx({x10:x1})-xx({x1}))/xb*pi)),1)

; W2b
 x10 = -3 + il*dspeed + 40
 x20 = 10 + il*dspeed + 40
 y10 = 52
 y20 = 67
 x1b = 3.
 x2b = 3.
 y1b = 3.
 y2b = 3.
 xb = x1b + x2b
 yb = y1b + y2b
 x1 = x10 + x1b
 x2 = x20 - x2b
 y1 = y10 + y1b
 y2 = y20 - y2b
 x10 = min( (/xx(nxx-1), max((/xx(0), x10/))/) )
 x20 = min((/xx(nxx-1), x20/))
 x1 = min( (/xx(nxx-1), max((/xx(0), x1/))/) )
 x2 = min((/xx(nxx-1), x2/))
 
 ws2b = vv1
 ws2b({y20:},:) = 0.
 ws2b({:y10},:) = 0.
 ws2b(:,{x20:}) = 0.
 ws2b(:,{:x10}) = 0.
 ws2b({y2:y20},{x10:x20}) = ws2b({y2:y20},{x10:x20})* \
     conform(ws2b({y2:y20},{x10:x20}),abs(cos((y({y2:y20})-y({y2}))/yb*pi)),0)
 ws2b({y10:y1},{x10:x20}) = ws2b({y10:y1},{x10:x20})* \
     conform(ws2b({y10:y1},{x10:x20}),abs(cos((y({y10:y1})-y({y1}))/yb*pi)),0)
 ws2b({y10:y20},{x2:x20}) = ws2b({y10:y20},{x2:x20})* \
     conform(ws2b({y10:y20},{x2:x20}),abs(cos((xx({x2:x20})-xx({x2}))/xb*pi)),1)
 ws2b({y10:y20},{x10:x1}) = ws2b({y10:y20},{x10:x1})* \
     conform(ws2b({y10:y20},{x10:x1}),abs(cos((xx({x10:x1})-xx({x1}))/xb*pi)),1)

; W2c
 x10 = -3 + il*dspeed - 40
 x20 = 10 + il*dspeed - 40
 y10 = 52
 y20 = 67
 x1b = 3.
 x2b = 3.
 y1b = 3.
 y2b = 3.
 xb = x1b + x2b
 yb = y1b + y2b
 x1 = x10 + x1b
 x2 = x20 - x2b
 y1 = y10 + y1b
 y2 = y20 - y2b
 x10 = max((/xx(0), x10/))
 x20 = max( (/0, min((/xx(nxx-1), x20/))/) )
 x1 = max((/xx(0), x1/))
 x2 = max( (/0, min((/xx(nxx-1), x2/))/) )
print("BBB   "+x10+"  "+x1+"  "+x2+"  "+x20)
 
 ws2c = vv1
 ws2c({y20:},:) = 0.
 ws2c({:y10},:) = 0.
 ws2c(:,{x20:}) = 0.
 ws2c(:,{:x10}) = 0.
 ws2c({y2:y20},{x10:x20}) = ws2c({y2:y20},{x10:x20})* \
     conform(ws2c({y2:y20},{x10:x20}),abs(cos((y({y2:y20})-y({y2}))/yb*pi)),0)
 ws2c({y10:y1},{x10:x20}) = ws2c({y10:y1},{x10:x20})* \
     conform(ws2c({y10:y1},{x10:x20}),abs(cos((y({y10:y1})-y({y1}))/yb*pi)),0)
 ws2c({y10:y20},{x2:x20}) = ws2c({y10:y20},{x2:x20})* \
     conform(ws2c({y10:y20},{x2:x20}),abs(cos((xx({x2:x20})-xx({x2}))/xb*pi)),1)
 ws2c({y10:y20},{x10:x1}) = ws2c({y10:y20},{x10:x1})* \
     conform(ws2c({y10:y20},{x10:x1}),abs(cos((xx({x10:x1})-xx({x1}))/xb*pi)),1)

 ws2 = ws2a + ws2b + ws2c

; W4
 x10 = 15
 x20 = 40
 y10 = 32
 y20 = 47
 x1b = 3.
 x2b = 3.
 y1b = 3.
 y2b = 3.
 xb = x1b + x2b
 yb = y1b + y2b
 x1 = x10 + x1b
 x2 = x20 - x2b
 y1 = y10 + y1b
 y2 = y20 - y2b
 
 ws4 = vv1
 ws4({y20:},:) = 0.
 ws4({:y10},:) = 0.
 ws4(:,{x20:}) = 0.
 ws4(:,{:x10}) = 0.
 ws4({y2:y20},{x10:x20}) = ws4({y2:y20},{x10:x20})* \
     conform(ws4({y2:y20},{x10:x20}),abs(cos((y({y2:y20})-y({y2}))/yb*pi)),0)
 ws4({y10:y1},{x10:x20}) = ws4({y10:y1},{x10:x20})* \
     conform(ws4({y10:y1},{x10:x20}),abs(cos((y({y10:y1})-y({y1}))/yb*pi)),0)
 ws4({y10:y20},{x2:x20}) = ws4({y10:y20},{x2:x20})* \
     conform(ws4({y10:y20},{x2:x20}),abs(cos((xx({x2:x20})-xx({x2}))/xb*pi)),1)
 ws4({y10:y20},{x10:x1}) = ws4({y10:y20},{x10:x1})* \
     conform(ws4({y10:y20},{x10:x1}),abs(cos((xx({x10:x1})-xx({x1}))/xb*pi)),1)

; Waves
 vv1 = 0.
 vv1 = ws2 + ws4

;--------------------------------------- cn level
 cnlev1 = fspan(-10,10,21)*3.
 cnlev2 = fspan(-10,10,21)*4.
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = ispan(2,65,3)
 cnfil(10:11) = 0

;--------------------------------------- contour
 wi = 128;40
 he = 154;48
 xf = 28 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+3)
 yf = 18 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+6)

 mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 lbopt = "w [mm s:S:-1:N:]"
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

;do il=0, nt-1

 contour1 = cnshadeplot(wks,xx,y,vv1*1.e3,cnlev1,cnfil,lbopt)
; cnaddline(contour1,wks,x,y,v0*1.e3,cnlev1,0,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.60
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.018
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.022
 end setvalues

 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,0,60,26,74,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; dr_txt(wks,mstr,xf(0)+wi/2.,yf(0)-2.5,22,5.5,0.)

 frame(wks)

end do

end



