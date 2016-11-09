load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 ds = 300
 n1s = toint(4*((7.25)-4))
 ipt0 = (/116,252/) + (/0,0/)
 iitv = (/1,-1/)
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/prt"+ds+"_hz5-mfx",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d"+ds+"km/prt_W__x1__08_00_00__12_00_00__6hrly.nc"
 fname2 = fdir+"/x1/prt_d"+ds+"km/prt_U__x1__08_00_00__12_00_00__6hrly.nc"
 fname3 = fdir+"/x1/prt_d"+ds+"km/prt_V__x1__08_00_00__12_00_00__6hrly.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f3 = addfile(fname3,"r")
 nx = 432
 ny = 648+2
 nz = 44
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

nt = dimsizes(n1s)

do il=0, nt-1

 dt = 6
 nhr = 24/dt
 tday  = n1s(il)/nhr
 thour = (n1s(il) - tday*nhr)*dt
 tday  = 4 + tday

 v0_w = f1->prt_W(n1s(il),:nz-1,:,:)
 v0_u = 0.5*( f2->prt_U(n1s(il),0:nz-2,:,:) + f2->prt_U(n1s(il),1:nz-1,:,:) )

 v1 = v0_w
 v1(1:,:,:) = v0_w(1:,:,:)*v0_u(:,:,:)
 v1(0,:,:) = 0
 h1 = f1->h(:nz-1)*990. + 10.
 v2 = f2->mean_U(n1s(il),:nz-1,:,:)
 v3 = f3->mean_V(n1s(il),:nz-1,:,:)
 h23 = f2->h(:nz-1)*990. + 10.

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

 nxy = 73
 xx = x(ipt0(0)-nxy/2*iitv(0):ipt0(0)+nxy/2*iitv(0):abs(iitv(0)))
 yy = y(ipt0(1)-nxy/2*iitv(1):ipt0(1)+nxy/2*iitv(1):abs(iitv(1)))

 v1i = new((/nz,nxy/),float)
 v2i = new((/nz,nxy/),float)
 v3i = new((/nz,nxy/),float)
 do ij=0, nxy-1
   v1i(:,ij) = v1(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
   v2i(:,ij) = v2(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
   v3i(:,ij) = v3(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
 end do
 v23i = (v2i*iitv(0) + v3i*iitv(1))/sqrt(iitv(0)^2+iitv(1)^2)
 
 ev = (h1/1.e3 - 0.252)*3.141592/2.
 sinlat = sin(yy*3.141592/180.)
 coslat = cos(yy*3.141592/180.)
 tz = 288.*(h1/1.e3)^(287.*5.e-3/9.80616)
 rho = new((/nz,nxy/),float)
 do k=0, nz-1
   ty = tz(k) + 0.75*(h1(k)/1.e3)*3.141592*35./287.*sin(ev(k))*sqrt(cos(ev(k)))* \
      ( (-2*sinlat(:)^6*(coslat(:)^2 + 1./3.) + 10./63.)*2.*35.*cos(ev(k))^1.5 + \
        (1.6*coslat(:)^3*(sinlat(:)^2 + 2./3.) - 3.141592/4.)*6371.229e3*7.29212e-5 )
   rho(k,:) = h1(k)*1.e2/287./ty(:)
 end do

 v1i = v1i*rho
 v1i(0,:) = 0.

;--------------------------------------- cn level
 scl = 1.e-3
; cnlev1 = ispan(-10,10,1)*1.0
 cnlev1 = (/-100,-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50,100/)
 lbopt = "[mPa]"
 cnlev2 = ispan(-10,10,1)*4
 cnlab = 1
; cnfil = ispan(2,65,3)
 cnfil = ispan(3,63,4)
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 20 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 20 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+9)

 mstr = "";"Day "+tday+"  "+sprinti("%2.2i",thour)+"::00"
 xstr = "lon [deg]"
 xstr2 = "lat [deg]"
 ystr = "p* [hPa]"
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

 contour1 = cnshadeplot(wks,xx,h1,v1i/scl,cnlev1,cnfil,lbopt)
; cnaddline(contour1,wks,xx,h23,v23i,cnlev2,cnlab,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
;  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
  "tmXBLabelsOn" : True
  "tmYLLabelsOn" : True
 end setvalues
 labt = ispan(toint(yy(0))/2*2,toint(yy(nxy-1)),2)
 minort = fspan(min(labt)-2,max(labt)+2,(max(labt)-min(labt)+4)*2+1)
 setvalues contour1
  "tmXUseBottom" : False
  "tmXTLabelsOn" : True
  "tmXTMode"   : "Explicit"
  "tmXTValues" : linint1(yy(::-1),xx(::-1),False,labt(::-1),0)
  "tmXTLabels" : labt(::-1)
  "tmXTMinorValues" : linint1(yy(::-1),xx(::-1),False,minort,0)
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-999,999,999,-999,False,True)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)
 
 draw(contour1)

 delete(contour1)

 dr_txt(wks,xstr2,xf(0)+wi/2.,yf(0)-6.5,22,6.0,0.)

end do

 frame(wks)

end

