load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "V"
 lev = 300
 ds = (/300/)
 n1s = toint(288*(( 7.5 )-4))-1
 var1name = "prt_"+var1+"_"+lev
 var2name = "mean_"+var1+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/prt_r"+"-exm1-"+var1+"_"+lev,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
nd = dimsizes(ds)

do il=0, nd-1

 fname1 = fdir+"/x1/prt_d"+ds(il)+"km/"+var1name+"__x1__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

 dt = 5
 nmin = 60/dt
 tday  = (n1s+1)/(nmin*24)
 tmin  = (n1s+1) - tday*(nmin*24)
 thour = tmin/nmin
 tmin  = (tmin - thour*nmin)*dt
 tday  = 4 + tday

 v1 = f1->$var1name$(n1s,:,:)
 v2 = f1->$var2name$(n1s,:,:)

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

;--------------------------------------- cn level
 if (var1 .eq. "W") then
   scl = 1.e-3
   cnlev1 = ispan(-10,10,1)*3
   cnlev2 = ispan(-10,10,1)*4
   cnlab = 0
   lbopt = "[mm s:S:-1:N:]"
 end if
 if (var1 .eq. "U") then
   scl = 1.
   cnlev1 = ispan(-10,10,1)*0.25
   cnlev2 = ispan(1,20,1)*3
   cnlab = 1
   lbopt = "[m s:S:-1:N:]"
 end if
 if (var1 .eq. "V") then
   scl = 1.
   cnlev1 = ispan(-10,10,1)*0.25
   cnlev2 = ispan(-10,10,1)*2
   cnlab = 1
   lbopt = "[m s:S:-1:N:]"
 end if
 cnfil = ispan(2,65,3)
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 128/1.5 / 4. *1.5
 he = 154 / 4. *0.8
 xf = 20 + (/0,1,2,3,4,0,1,2,3,4,0,1,2,3,4/)*(wi+6)
 yf = 20 + (/0,0,0,0,0,1,1,1,1,1,2,2,2,2,2/)*(he+12)

 mstr = "";"Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 if (yf(il) .ne. max(yf(:nd-1))) then
   xstr = ""
 end if
 if (xf(il) .ne. min(xf(:nd-1))) then
   ystr = ""
 end if
 if ( xf(il) .ne. min(xf(:nd-1)) .or. yf(il) .ne. max(yf(:nd-1)) ) then
   lbopt = ""
 end if

 contour1 = cnshadeplot(wks,x,y,v1/scl,cnlev1,cnfil,lbopt)
 cnaddline(contour1,wks,x,y,v2/scl,cnlev2,cnlab,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.020
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.025
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,40,34,74,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; dr_txt(wks,mstr,xf(0)+wi/2.,yf(0)-2.5,22,5.5,0.)

end do

 frame(wks)

end

