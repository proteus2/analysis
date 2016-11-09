load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 ds = 300
 ftail = (/"10_00_00__11_12_00","10_00_00__11_12_00","10_00_00__11_12_00", \
           "11_00_00__12_00_00","11_00_00__12_00_00"/)
 ttt = (/7,7,7,7.5,7.5/)
 n1s = toint(24*(ttt-(/6,6,6,7,7/)))
 n1ms = toint(288*(ttt-4))
 phi = (/-30,10,45,60,-20/)

 xy0p = (/ (/42.1,50/), \
           (/29,59/), \
           (/54,59/), \
           (/30,41.0/), \
           (/28.9,45.1/) /)
 iitv0 = (/ (/1,0/), \
            (/1,0 /), \
            (/1,1 /), \
            (/1,0 /), \
            (/1,-1/) /)
 xyl0 = (/6,6,6,24,4.2/)

 ys0 = (/48.00925,48.00925,48.00925,21.06481,21.06481/)
 ipt00 = (/303,303,303,12,12/)
 scl = (/10,7,6,8,10/)
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/"+var1name+"_hz-5x",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 nw = dimsizes(n1s)

do iw=0, nw-1
 ys = ys0(iw)
 iitv = iitv0(iw,:)
;--------------------------------------- read data
 nx = 432
 ny = 270
 nym = 648+2
 nz = 85 ;101
 dz = 250
 dx = 1./10.8
 dy = 1./10.8
 xs = 20.0463

 v1 = new((/nz,nym,nx/),float)
; v1 = new((/nz,ny ,nx/),float)
 v2 = new((/nz,nym,nx/),float)
 v3 = new((/nz,nym,nx/),float)

 z = ispan(0,(nz-1)*dz,dz)
 do k=0, nz-1
   fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d"+ds+"km/"+var1name+"__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
;   fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d"+ds+"km/recon_xy/"+var1name+"__z"+sprinti("%5.5i",z(k))+"__x6__"+ftail(iw)+"__1hrly-w"+(iw+1)+".nc"
   fname2 = fdir+"/x6_5min/U_z/prt_d"+ds+"km/mean_U__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fname3 = fdir+"/x6_5min/V_z/prt_d"+ds+"km/mean_V__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   f1 = addfile(fname1,"r")
   f2 = addfile(fname2,"r")
   f3 = addfile(fname3,"r")

   v1(k,:,:) = f1->$var1name$(n1ms(iw),:,:)
;   v1(k,:,:) = f1->$var1name$(n1s(iw),:,:)
   v2(k,:,:) = f2->mean_U(n1ms(iw),:,:)
   v3(k,:,:) = f3->mean_V(n1ms(iw),:,:)
 end do

;--------------------------------------- process 0
 x = xs + ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

 nxy = toint(xyl0(iw)/dx)

 ipt0 = (/ closest_val(xy0p(iw,0),x), closest_val(xy0p(iw,1),y) /)
 ipt0m = ipt0 + (/0,ipt00(iw)/)

 xx = x(ipt0(0):ipt0(0)+(nxy-1)*iitv(0):abs(iitv(0)))
 yy = y(ipt0(1):ipt0(1)+(nxy-1)*iitv(1):abs(iitv(1)))

 v1i = new((/nz,nxy/),float)
 v2i = new((/nz,nxy/),float)
 v3i = new((/nz,nxy/),float)
 ec = exp(-z/(2.*7.e3))
 do ij=0, nxy-1
   v1i(:,ij) = v1(:,ipt0m(1)+ij*iitv(1),ipt0m(0)+ij*iitv(0))*ec(:)
;   v1i(:,ij) = v1(:,ipt0 (1)+ij*iitv(1),ipt0 (0)+ij*iitv(0))*ec(:)
   v2i(:,ij) = v2(:,ipt0m(1)+ij*iitv(1),ipt0m(0)+ij*iitv(0))
   v3i(:,ij) = v3(:,ipt0m(1)+ij*iitv(1),ipt0m(0)+ij*iitv(0))
 end do
 v23i = v2i*cos(phi(iw)/180.*3.141592) + v3i*sin(phi(iw)/180.*3.141592)
 
;--------------------------------------- cn level
; cnlev1 = ispan(-10,10,1)
 cnlev1 = ispan(-5,5,1)*1.
 cnlev1 = cnlev1*scl(iw)
 lbopt = "[mm s:S:-1:N:]"
 cnlev2 = ispan(-10,10,1)*4
 cnlab = 1
; cnfil = ispan(2,65,3)
 cnfil = array_append_record( ispan(2,32,6), ispan(35,65,6), 0 )
 nfil = dimsizes(cnfil)
; cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 40
 he = 30
 xf = 20 + (/0,1,2,0,1,2/)*(wi+8)
 yf = 20 + (/0,0,0,1,1,1/)*(he+25)

 tday = toint(ttt(iw))
 tmp = (ttt(iw) - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)
 mstr = "Day "+tday+",  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 if (iw .eq. 0) then
   mstr1 = mstr
 end if
 if (iw .eq. 3) then
   mstr2 = mstr
 end if
 xstr = "longitude [deg]"
 xstr2 = "latitude [deg]"
 ystr = "z [km]"
 if (iitv(1) .eq. 0) then
   xstr2 = ""
 end if
 if (xf(iw) .ne. xf(0)) then
   ystr = ""
 end if

 contour1 = cnshadeplot(wks,xx,z/1.e3,v1i*1.e3,cnlev1,cnfil,lbopt)
 cnaddline(contour1,wks,xx,z/1.e3,v23i,cnlev2,cnlab,"D",(/0,0,14/),"T",0.5,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
;  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
;  "lbTitleFontHeightF" : 0.015
 end setvalues
; if (iitv(1) .ne. 0) then
; setvalues contour1
;  "tmXBLabelsOn" : True
;  "tmYLLabelsOn" : True
; end setvalues
; labt = ispan(toint(yy(0))/2*2,toint(yy(nxy-1)),2)
; minort = fspan(min(labt)-2,max(labt)+2,(max(labt)-min(labt)+4)*2+1)
; setvalues contour1
;  "tmXUseBottom" : False
;  "tmXTLabelsOn" : True
;  "tmXTMode"   : "Explicit"
;  "tmXTValues" : linint1(yy(::iitv(1)),xx(::iitv(1)),False,labt(::iitv(1)),0)
;  "tmXTLabels" : labt(::iitv(1))
;  "tmXTMinorValues" : linint1(yy(::iitv(1)),xx(::iitv(1)),False,minort,0)
; end setvalues
; delete(labt)
; delete(minort)
; end if

 set_size(contour1,xf(iw),yf(iw),wi,he)
 axis_range(contour1,-999,999,-999,999,False,False)
 axis_str(contour1,"",xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; if (iitv(1) .ne. 0) then
;   dr_txt(wks,xstr2,xf(iw)+wi/2.,yf(iw)-6.5,22,6.0,0.)
; end if

 dr_txt(wks,(iw+1),xf(iw)+1.5*1.8,yf(iw)+1.5*1.8,22,4.5*1.8,0.)

 delete([/xx,yy,v1i,v2i,v3i,v23i/])

end do  ; iw

 dr_txt(wks,"(a)",xf(0)-6*1.8,yf(0)-3.1*1.8,25,6*1.8,0.)
 dr_txt(wks,"(b)",xf(0)-6*1.8,yf(3)-3.1*1.8,25,6*1.8,0.)
 dr_txt(wks,mstr1,xf(0)+6*1.8,yf(0)-2.8*1.8,22,5*1.8,0.)
 dr_txt(wks,mstr2,xf(0)+6*1.8,yf(3)-2.8*1.8,22,5*1.8,0.)

 frame(wks)

end

