load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 xya = "a"
 var1 = "W"
 var2 = "U"
 var3 = "V"
 lev = 300
 ds = 300
 n1s = toint(288*(fspan(6.5,8-2.25/24.,16)-4))-1
 var1name = "prt_"+var1+"_"+lev
 var2name = "prt_"+var2+"_"+lev
 var3name = "prt_"+var3+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/prt"+ds+"_xy-mf"+xya+"_"+lev,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d"+ds+"km/"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname2 = fdir+"/x1/prt_d"+ds+"km/"+var2name+"__x1__08_00_05__12_00_00.nc"
 fname3 = fdir+"/x1/prt_d"+ds+"km/"+var3name+"__x1__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f3 = addfile(fname3,"r")
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
 if (xya .ne. "y") then
   v2 = f2->$var2name$(n1s(il),:,:)
   mfx = v1*v2
   mf = mfx
 end if
 if (xya .ne. "x") then
   v3 = f3->$var3name$(n1s(il),:,:)
   mfy = v1*v3
   mf = mfy
 end if
 if (xya .eq. "a") then
   mf = sqrt(mfx^2 + mfy^2)
 end if

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

if (il .eq. 0) then
 ev = (lev/1.e3 - 0.252)*3.141592/2.
 sinlat = sin(y*3.141592/180.)
 coslat = cos(y*3.141592/180.)
 tz = 288.*(lev/1.e3)^(287.*5.e-3/9.80616)

 ty = tz + 0.75*(lev/1.e3)*3.141592*35./287.*sin(ev)*sqrt(cos(ev))* \
      ( (-2*sinlat^6*(coslat^2 + 1./3.) + 10./63.)*2.*35.*cos(ev)^1.5 + \
        (1.6*coslat^3*(sinlat^2 + 2./3.) - 3.141592/4.)*6371.229e3*7.29212e-5 )

 rho = conform( mf, lev*1.e2/287./ty, 0 )
end if

 mf = mf*rho*1.e3  ; [mPa]

;--------------------------------------- cn level
 if (xya .eq. "a") then
   cnlev1 = ispan(1,10,1)*3
   cnfil = ispan(35,65,3)
   cnfil(0) = 0
 else
   cnlev1 = ispan(-10,10,1)*2
   cnfil = ispan(2,65,3)
   nfil = dimsizes(cnfil)
   cnfil(nfil/2-1:nfil/2) = 0
;   cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 end if

;--------------------------------------- contour
 wi = 128/1.5 / 4.
 he = 154 / 4.
 xf = 20 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 20 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+9)

 mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 lbopt = "[mPa]"
 if (yf(il) .ne. max(yf(:nt-1))) then
   xstr = ""
 end if
 if (xf(il) .ne. min(xf(:nt-1))) then
   ystr = ""
 end if
; if ( xf(il) .ne. min(xf(:nt-1)) .or. yf(il) .ne. max(yf(:nt-1)) ) then
 if ( xf(il) .ne. max(xf(:nt-1)) .or. yf(il) .ne. min(yf(:nt-1)) ) then
   lbopt = ""
 end if

 contour1 = cnshadeplot(wks,x,y,mf,cnlev1,cnfil,lbopt)
; cnaddline(contour1,wks,x,y,v2/scl,cnlev2,cnlab,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.040
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.045
  "tmXBLabelsOn" : True
  "tmYLLabelsOn" : True
 end setvalues
 setvalues contour1
  "lbOrientation" : "Vertical"
  "pmLabelBarOrthogonalPosF" : -1.10
  "pmLabelBarParallelPosF"   : 1.3
  "pmLabelBarWidthF"         : 0.2
  "pmLabelBarHeightF"        : 2.1
  "lbLabelFontHeightF" : 0.040
  "lbTitleOffsetF"     : 0.006
  "lbTitleFontHeightF" : 0.045
 end setvalues
; if (mod(il,4) .ne. 0) then
;  setvalues contour1
;   "tmYLLabelsOn" : False
;  end setvalues
; end if
 if (il .lt. 12) then
  setvalues contour1
   "tmXBLabelsOn" : False
  end setvalues
 end if

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,40,26,74,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)
 
 draw(contour1)

 delete(contour1)

; dr_txt(wks,mstr,xf(0)+wi/2.,yf(0)-2.5,22,5.5,0.)

end do

 frame(wks)

end

