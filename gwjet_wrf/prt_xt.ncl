load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 var0 = "U"
 lev = 300
 lat = 60
; n1s = ispan(0,1151,6)
 var1name = var1+"_"+lev
 var0name = var0+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
; wks = cr_wks("fig/prt_xt-y"+lat+"-"+var1name,600,0)
 wks = cr_wks("fig/prt_xt-y"+lat+"W2-"+var1name,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/prt_"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname1b= fdir+"/x2/prt_d300km/prt_"+var1name+"__x2__08_00_05__12_00_00.nc"
 fname0 = fdir+"/x1/prt_d300km/prt_"+var0name+"__x1__08_00_05__12_00_00.nc"
 fname0b= fdir+"/x2/prt_d300km/prt_"+var0name+"__x2__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 f1b= addfile(fname1,"r")
 f0 = addfile(fname0,"r")
 f0b= addfile(fname0,"r")

 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

 nt = 1152 ;dimsizes(n1s)

 t = 4. + ispan(1,nt,1)/(24.*12.)
; dt = 5
; nmin = 60/dt
; tday  = (n1s(il)+1)/(nmin*24)
; tmin  = (n1s(il)+1) - tday*(nmin*24)
; thour = tmin/nmin
; tmin  = (tmin - thour*nmin)*dt
; tday  = 4 + tday

 j1 = toint((lat - ys)/dy)

 v1 = f1->$("prt_"+var1name)$(:,j1,:)
 v1b= f1b->$("prt_"+var1name)$(:,j1,:nx/2-1)
 v0 = f0->$("mean_"+var0name)$(:,j1,:)
 v0b= f0b->$("mean_"+var0name)$(:,j1,:nx/2-1)

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx

 xx = ispan(0,nx+nx/2-1,1)*dx

 vv1 = new((/nt,nx+(nx/2)/),float)
 vv1(:,:nx-1) = v1 (:,:)
 vv1(:,nx:  ) = v1b(:,:)

 vv0 = new((/nt,nx+(nx/2)/),float)
 vv0(:,:nx-1) = v0 (:,:)
 vv0(:,nx:  ) = v0b(:,:)

;--------------------------------------- cn level
 cnlev1 = ispan(-10,10,1)*4.
 cnlev0 = ispan(-20,50,1)*1.
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = ispan(2,65,3)

;--------------------------------------- contour
 wi = 40
 he = 48
 xf = 28 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+3)
 yf = 18 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+6)

 mstr = "" ;"Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lon [deg]"
 ystr = "Day"
 lbopt = "w [mm s:S:-1:N:]"
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)


 contour1 = cnshadeplot(wks,xx,t,vv1*1.e3,cnlev1,cnfil,lbopt)
 cnaddline(contour1,wks,xx,t,vv0,cnlev0,1,"D",0,"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "tmYLMinorPerMajor" : 3
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
 axis_range(contour1,0,60,4,999,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; dr_txt(wks,mstr,xf(0)+wi/2.,yf(0)-2.5,22,5.5,0.)

 frame(wks)

end



