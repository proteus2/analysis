load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = "merra"
 varname = (/"u","v","omega","t"/)  ; do not change order

 ew = "e"
 nw = 3

 yy = "1980"
 mm = "01"
 ndays = 90
 p0 = 70
 n0 = (/1,1,1,1/)*90 ;(/115,115,60,60/)   ; 60 - 120
 lat0 = 27
 h_scale = 6.6e3

 nk = 10
 no = ndays/1  ; 1-day  ; cut-off period

 wks = cr_wks("fig/rce-ra1.ps",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

 pi = 3.141592
 deg2rad = pi/180.

;--------------------------------------- read data
 nv = dimsizes(varname)

 fdir  = "/data4/kyh/dat/MERRA/fcoef_ko"
 fname = fdir+"/"+yy+"/"+expname+".fft_"+varname+"_ko."+yy+"."+mm+".nc"
print(fname(2)+"")

; axis
 f0 = addfile(fname(0),"r")
 s_v0 = "fcr_"+varname(0)
 lat = f0->$(f0->$s_v0$!1)$({-lat0:lat0})
 lath = f0->$(f0->$s_v0$!1)$({0:lat0})
 ny2 = dimsizes(lat)
 ny = dimsizes(lath)
 nx = nk*2
 noa = dimsizes(f0->$(f0->$s_v0$!2)$)

 f1 = addfile("$DATD/clim/MERRA/monthly/"+yy+"/"+expname+".u_yp."+yy+mm+".nc","r")
 f2 = addfile("$DATD/clim/MERRA/monthly/"+yy+"/"+expname+".t_yp."+yy+mm+".nc","r")
 um = f1->u(0,{p0},{-lat0-1.25:lat0+1.25})
 pa = tofloat( f2->levels )
 kp0 = minind(abs(log(pa/p0)))
 dptmdz = ( f2->t(0,kp0+1,{-lat0:lat0})*(1.e3/pa(kp0+1))^0.286 - \
            f2->t(0,kp0-1,{-lat0:lat0})*(1.e3/pa(kp0-1))^0.286 )/ \
          log(pa(kp0-1)/pa(kp0+1))/h_scale

 dumdy = new(ny2,float)
 coslat = tofloat(cos(f1->latitude({-lat0-1.25:lat0+1.25})*deg2rad))
 do j=0, ny2-1
   dumdy(j) = (um(j+2)*coslat(j+2) - um(j)*coslat(j)) \
              /(2.5*deg2rad)/6371.e3/coslat(j+1)
 end do

 if (ew .eq. "w") then
   io0 = 1
   io9 = no
 else 
   io0 = noa-1
   io9 = noa-no
 end if

 cf = new((/2,nv,2,ny,no,nk/),float)
 cf2 = new((/2,nv,ny2,no,nk/),float)
 si = (/"r","i"/)
 do iv=0, nv-1
   f0 = addfile(fname(iv),"r")
   do i=0, 1
     s_v0 = "fc"+si(i)+"_"+varname(iv)
     cf2(i,iv,:,:,:) = f0->$s_v0$({p0},{-lat0:lat0},io0:io9,1:nk)
   end do
 end do
 cf2(:,2,:,:,:) = cf2(:,2,:,:,:)*(-h_scale/(p0*1.e2))
 cf2(:,3,:,:,:) = cf2(:,3,:,:,:)*(1.e3/p0)^0.286

; Kelvin, E-MRG
; do n=no/2+1, no-1
; do i=1, nk
;   do j=ny, ny2-1
;     if (abs(mf(j,n,i)) .lt. abs(hf(j,n,i))) then
;       cf2(:,:,j:,n,i) = 0.
;     end if
;     if (abs(mf(ny2-1-j,n,i)) .lt. abs(hf(ny2-1-j,n,i))) then
;       cf2(:,:,:ny2-1-j,n,i) = 0.
;     end if
;;     if ( abs(mf(j,n,i)) .lt. abs(hf(j,n,i)) .or. \
;;          abs(mf(ny2-1-j,n,i)) .lt. abs(hf(ny2-1-j,n,i))) then
;;       cf2(:,:,j:,n,i) = 0.
;;       cf2(:,:,:ny2-1-j,n,i) = 0.
;;     end if
;   end do
; end do
; end do

;; W-MRG
; do n=oc_wmi(0), oc_wmi(1) ; 0, no/2-1
; do i=1, nk
;   do j=ny, ny2-1
;     if (mf(j,n,i)*fhf(j,n,i) .gt. 0.) then
;       cf2(:,:,j:,n,i) = 0.
;     end if
;     if (mf(ny2-1-j,n,i)*fhf(ny2-1-j,n,i) .gt. 0.) then
;       cf2(:,:,:ny2-1-j,n,i) = 0.
;     end if
;   end do
;;   if ( sum(abs(cf2(:,:,ny-2,n,i))) .eq. 0. .and. \
;;        sum(abs(cf2(:,:,ny  ,n,i))) .eq. 0. ) then
;;     cf2(:,:,:,n,i) = 0.
;;   end if
; end do
; end do

 do iv=0, nv-1
 do i=0, 1
   cf(i,iv,:,:,:,:) = 0.5*(/ cf2(i,iv,ny:ny2-1,:,:)+cf2(i,iv,ny-1:0,:,:), \
                             cf2(i,iv,ny:ny2-1,:,:)-cf2(i,iv,ny-1:0,:,:) /)
   tmpcf = cf(i,1,0,:,:,:)
   cf(i,1,0,:,:,:) = cf(i,1,1,:,:,:)
   cf(i,1,1,:,:,:) = tmpcf
 end do
 end do 

; mf = -uw  ;  hf = vth/th_z*(f-u_y)  ;  hf = vth/th_z*f
 mf0 = -(cf(0,0,:,:,:,:)*cf(0,2,:,:,:,:) + cf(1,0,:,:,:,:)*cf(1,2,:,:,:,:))
 hf0 = cf(0,1,:,:,:,:)*cf(0,3,:,:,:,:) + cf(1,1,:,:,:,:)*cf(1,3,:,:,:,:)

 mf = new((/2,ny2,no,nk/),float)
 mf(:,ny:,:,:) = mf0(:,:,:,:)
 mf(:,0:ny-1,:,:) = mf0(:,ny-1:0,:,:)
 hf = new((/2,ny2,no,nk/),float)
 hf(:,ny:,:,:) = hf0(:,:,:,:)
 hf(:,0:ny-1,:,:) = hf0(:,ny-1:0,:,:)*(-1.)
 fhf = hf

 fcor = 2.*7.292116e-5*sin(lat*deg2rad)
 do j=0, ny2-1
   fhf(:,j,:,:) = fhf(:,j,:,:)/dptmdz(j)*fcor(j)
   hf (:,j,:,:) = hf (:,j,:,:)/dptmdz(j)*(fcor(j)-dumdy(j))
 end do

 tmp1 = mf
 tmp2 = hf
 do j=2, ny2-3
   mf(:,j,:,:) = dim_avg_n( tmp1(:,j-2:j+2,:,:), 1 )
   hf(:,j,:,:) = dim_avg_n( tmp2(:,j-2:j+2,:,:), 1 )
 end do

 fz = mf + hf

 rt = fz

 rt(0,:,:,:) = log10(abs(mf(0,:,:,:)/hf(0,:,:,:)))

 rt(1,:,:,:) = -mf(1,:,:,:)/fhf(1,:,:,:)

;--------------------------------------- process 0
 vv = new((/nv,no,nk,2,ny,nx/),float)
 xi = ispan(0,nx-1,1)*1./tofloat(nx)*2.*pi
 do iv=0, nv-1
 do io=0, no-1
 do i=0, nk-1
 do s=0, 1
 do j=0, ny-1
   vv(iv,io,i,s,j,:) = \
      cf(0,iv,s,j,io,i)*cos(xi) - cf(1,iv,s,j,io,i)*sin(xi)
 end do
 end do
 end do
 end do
 end do
 vv = vv*2.

 padd = new((/2,no,nk,ny2,nx/),float)
 do s=0, 1
 do io=0, no-1
 do i=0, nk-1
 do j=0, ny2-1
   padd(s,io,i,j,:) = rt(s,j,io,i)
 end do
 end do
 end do
 end do

 vv2 = new((/nv,no,nk,2,ny2,nx+1/),float)
 vv2(:,:,:,:,ny:,:nx-1) = vv(:,:,:,:,:,:)
; symm
 vv2(:,:,:,0,0:ny-1,:nx-1) = vv(:,:,:,0,ny-1:0,:)
; anti-symm
 vv2(:,:,:,1,0:ny-1,:nx-1) = vv(:,:,:,1,ny-1:0,:)*(-1.)
; v'
 vv2(1,:,:,:,0:ny-1,:nx-1) = vv2(1,:,:,:,0:ny-1,:nx-1)*(-1.)

 vv2(:,:,:,:,:,nx) = vv2(:,:,:,:,:,0)

 lon = new((/nk,nx+1/),float)
 do i=0, nk-1
   lon(i,:) = ( fspan(-180.,180.,nx+1) - 0.5 )/tofloat(i+1)
 end do

; tmp = vv2(:,2,:,:)
; do j=1, ny2-2
;   vv2(:,2,j,:) = (2.*tmp(:,j,:) + tmp(:,j+1,:) + tmp(:,j-1,:))/4.
; end do

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlevh = fspan(-1.9,1.9,20)*2.

 cnlevr = log10( (/0.1,0.2,0.5,1.,2.,5.,10/) )
 cnfilr = (/24,25,26,0,0,11,10,9/)

;--------------------------------------- contour
 wi = 40
 he = 20
 xf = 8 + (/0,2,1,3/)*46
 yf = 2 + ispan(0,9,1)*24

 mstr = ""
 xstr = ""
 ystr = ""
 lbtxt = (/"","","","(mm s:S:-1:N:)"/)

do io=0, 14 ;no-1

  ioo = (io+1)*3 - 1

do i=0, nk-1

  x0 = min(lon(i,:))
  x9 = max(lon(i,:))

do s=0, 1

 contour1 = cnshadeplot(wks,lon(i,:),lat,vv2(2,ioo,i,s,:,:)*1.e3,cnlev/(i+1),cnfil,lbtxt(s))
 cnaddline(contour1,wks,lon(i,:),lat,vv2(3,ioo,i,s,:,:),cnlevh/(i+1),0,"D",0,"T",0.5,"C",(/colindr,1,colindb/))
; setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.15
;  "pmLabelBarParallelPosF"   : 0.50-0.40
;  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.0
;  "pmLabelBarHeightF"        : 0.15
;  "lbLabelFontHeightF" : 0.025
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.10
;  "lbTitleFontHeightF" : 0.030
; end setvalues

 set_size(contour1,xf(s),yf(i),wi,he)
 axis_range(contour1,x0,x9,-22,22,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vv2(0:1,ioo,i,s,2::3,::2)
  "vfXArray"    : lon(i,::2)
  "vfYArray"    : lat(2::3)
 end create
 ref = 1./(i+1)
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(s),yf(s),wi,he)
 axis_range(vec1,x0,x9,-22,22,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

 contour2 = cnshadeplot(wks,lon(i,:),lat,padd(s,ioo,i,:,:),cnlevr,cnfilr,"")
 cnaddline(contour2,wks,lon(i,:),lat,padd(s,ioo,i,:,:),cnlevr,1,"D",0,"T",(/0.5,1.5,0.5/),"C",(/colindr,1,colindb/))

 set_size(contour2,xf(2+s),yf(i),wi,he)
 axis_range(contour2,x0,x9,-22,22,False,False)
 axis_str(contour2,"","","")

 cnfinalize(contour2,wks)

 draw(contour2)

end do  ; s

end do  ; i

 dr_txt2(wks,ioo+"-"+(90/ioo)+"d",xf(3),yf(0),22,6,0.,"CenterLeft",1)

 frame(wks)

end do  ; io

end

