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
 ds = 300
 n1s = (/ toint(288*((7.25)-4)), toint(288*((7.25)-4)), toint(288*((7.25)-4)), \
          toint(288*((7.25)-4)), toint(288*((7.25)-4)) /)
 ipt00 = (/ (/356,324/) + (/0,0/), \
            (/216,448/) + (/0,0/), \
            (/62 ,438/) + (/0,0/), \
            (/162,244/) + (/0,0/), \
            (/116,252/) + (/0,0/) /)
 iitv0 = (/ (/1,-1/), \
            (/1,0 /), \
            (/1,1 /), \
            (/1,1 /), \
            (/1,-1/) /)
 nxy0 = (/81,101,101,69,73/)
 scl = (/1,2,1,6,2/)
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/"+var1name+"_hz-5",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 nw = dimsizes(n1s)

do iw=0, nw-1
 ipt0 = ipt00(iw,:)
 iitv = iitv0(iw,:)
 nxy = nxy0(iw)
;--------------------------------------- read data
 nx = 432
 ny = 648+2
 nz = 101
 dz = 250
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

 dt = 6
 nhr = 24/dt
 tday  = n1s(iw)/nhr
 thour = (n1s(iw) - tday*nhr)*dt
 tday  = 4 + tday

 v1 = new((/nz,ny,nx/),float)
 v2 = new((/nz,ny,nx/),float)
 v3 = new((/nz,ny,nx/),float)

 z = ispan(0,(nz-1)*dz,dz)
 do k=0, nz-1
   fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d"+ds+"km/"+var1name+"__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fname2 = fdir+"/x6_5min/U_z/prt_d"+ds+"km/mean_U__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fname3 = fdir+"/x6_5min/V_z/prt_d"+ds+"km/mean_V__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   f1 = addfile(fname1,"r")
   f2 = addfile(fname2,"r")
   f3 = addfile(fname3,"r")

   v1(k,:,:) = f1->$var1name$(n1s(iw),:,:)
   v2(k,:,:) = f2->mean_U(n1s(iw),:,:)
   v3(k,:,:) = f3->mean_V(n1s(iw),:,:)
 end do

;--------------------------------------- process 0
 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

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
 
;--------------------------------------- cn level
 if (var1 .eq. "W") then
   cnlev1 = ispan(-10,10,1)*0.005
   cnlev1 = cnlev1*scl(iw)
   lbopt = "[m s:S:-1:N:]"
 end if
 if (var1 .eq. "U") then
   cnlev1 = ispan(-10,10,1)*0.25
   lbopt = "[m s:S:-1:N:]"
 end if
 if (var1 .eq. "V") then
   cnlev1 = ispan(-10,10,1)*0.25
   lbopt = "[m s:S:-1:N:]"
 end if
 cnlev2 = ispan(-10,10,1)*4
 cnlab = 1
 cnfil = ispan(2,65,3)
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 20 + (/0,1,2,0,1,2/)*(wi+8)
 yf = 20 + (/0,0,0,1,1,1/)*(he+27)

 mstr = "";"Day "+tday+"  "+sprinti("%2.2i",thour)+"::00"
 xstr = "lon [deg]"
 xstr2 = "lat [deg]"
 ystr = "z [km]"
 if (iitv(1) .eq. 0) then
   xstr2 = ""
 end if
 if (xf(iw) .ne. xf(0)) then
   ystr = ""
 end if

 contour1 = cnshadeplot(wks,xx,z/1.e3,v1i,cnlev1,cnfil,lbopt)
; cnaddline(contour1,wks,xx,z/1.e3,v23i,cnlev2,cnlab,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues
 if (iitv(1) .ne. 0) then
 setvalues contour1
  "tmXBLabelsOn" : True
  "tmYLLabelsOn" : True
 end setvalues
 labt = ispan(toint(yy(0))/2*2,toint(yy(nxy-1)),2)
 minort = fspan(min(labt)-2,max(labt)+2,(max(labt)-min(labt)+4)*2+1)
 setvalues contour1
  "tmXUseBottom" : False
  "tmXTLabelsOn" : True
  "tmXTMode"   : "Explicit"
  "tmXTValues" : linint1(yy(::iitv(1)),xx(::iitv(1)),False,labt(::iitv(1)),0)
  "tmXTLabels" : labt(::iitv(1))
  "tmXTMinorValues" : linint1(yy(::iitv(1)),xx(::iitv(1)),False,minort,0)
 end setvalues
 delete(labt)
 delete(minort)
 end if

 set_size(contour1,xf(iw),yf(iw),wi,he)
 axis_range(contour1,-999,999,-999,999,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)
 
 draw(contour1)

 delete(contour1)

 if (iitv(1) .ne. 0) then
   dr_txt(wks,xstr2,xf(iw)+wi/2.,yf(iw)-6.5,22,6.0,0.)
 end if

 delete([/xx,yy,v1i,v2i,v3i,v23i/])

end do  ; iw

 frame(wks)

end

