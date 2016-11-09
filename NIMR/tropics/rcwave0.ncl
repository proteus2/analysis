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
 ndays = 90
 z0 = 18.59;20.95 ;20.95 ;17.41 ;20.95
 n0 = (/115,115,60,60/)   ; 60 - 120
 lat0 = 25

 kc = (/3,3/)
 nk = 20
 no = 1*(2*ndays)  ; 1-day  ; cut-off period

 oc_pn = (/0.5,0.5/) ; 0.75  ; cut-off freq. for waves [cyc/day]

 oc_k  = (/ 1./90. , 0.5 /)
 oc_r1 = (/ 1./90. , 0.2 /)
 oc_em = (/ 0.25 , 0.5 /)
 oc_wm = (/ 0.15 , 0.4 /)

 wks = cr_wks("fig/rcwave0",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

;--------------------------------------- read data
 nv = dimsizes(varname)
 nw = 4

 fdir  = "$DATD/L60CGW-t/fcoef_ko"
 fname = fdir+"/"+yy+"/uan"+expname+".fft_"+varname+"_ko."+yy+"."+mm+".nc"
 fname(3) = fdir+"/"+yy+"/uan"+expname+".fft_"+varname(3)+"_ko_zr."+yy+"."+mm+".nc"
print(fname+"")

; axis
 f0 = addfile(fname(0),"r")
 s_v0 = "fcr_"+varname(0)
 lat = f0->$(f0->$s_v0$!1)$({-lat0:lat0})
 lath = f0->$(f0->$s_v0$!1)$({0:lat0})
 ny2 = dimsizes(lat)
 ny = dimsizes(lath)
 nx = nk*2
 noa = dimsizes(f0->$(f0->$s_v0$!2)$)
 oc_pn_i = (/ toint(oc_pn(1)*ndays), no - toint(oc_pn(0)*ndays) /)

 f1 = addfile("$DATD/L60CGW-t/xm/xm.uan"+expname+"a.pj_"+yy+mm+".nc","r")
 rhob = f1->rho(0,{z0*1.e3},{-lat0:lat0})

 cf = new((/2,nv,2,ny,no,nk+1/),float)
 cf2 = new((/2,nv,ny2,no,nk+1/),float)
 si = (/"r","i"/)
 do iv=0, nv-1
   f0 = addfile(fname(iv),"r")
   do i=0, 1
     s_v0 = "fc"+si(i)+"_"+varname(iv)
if (varname(iv) .ne. "rhodz_dt") then
     cf2(i,iv,:,:no/2,:) = f0->$s_v0$({z0*1.e3},{-lat0:lat0},:no/2,:nk)
     cf2(i,iv,:,no/2+1:,:) = f0->$s_v0$({z0*1.e3},{-lat0:lat0},noa-no/2+1:,:nk)
else
     cf2(i,iv,:,:no/2,:) = 0.5*( f0->$s_v0$({z0*1.e3-600.},{-lat0:lat0},:no/2,:nk) + f0->$s_v0$({z0*1.e3+600.},{-lat0:lat0},:no/2,:nk) )
     cf2(i,iv,:,no/2+1:,:) = 0.5*( f0->$s_v0$({z0*1.e3-600.},{-lat0:lat0},noa-no/2+1:,:nk) + f0->$s_v0$({z0*1.e3+600.},{-lat0:lat0},noa-no/2+1:,:nk) )
end if
     cf(i,iv,:,:,:,:) = 0.5*(/ cf2(i,iv,ny-1:ny2-1,:,:)+cf2(i,iv,ny-1:0,:,:), \
                               cf2(i,iv,ny-1:ny2-1,:,:)-cf2(i,iv,ny-1:0,:,:) /)
     if (varname(iv) .eq. "rhov") then
       tmp = cf(i,iv,0,:,:,:)
       cf(i,iv,0,:,:,:) = cf(i,iv,1,:,:,:)
       cf(i,iv,1,:,:,:) = tmp
     end if
;     cf(i,iv,:,:,:,:) = (/ cf2(i,iv,0:ny-1,:,:), cf2(i,iv,0:ny-1,:,:) /)
;     cf(i,iv,:,:,:,:) = (/ cf2(i,iv,ny-1:ny2-1,:,:), cf2(i,iv,ny-1:ny2-1,:,:) /)
   end do
 end do
 cf(:,:,:,:,no/2,:) = 0.
 cf2(:,:,:,no/2,:) = 0.

;--------------------------------------- process 0
 rho0 = avg(rhob({-10:10}))
 oc_ki = no - toint( (/ oc_k(1), oc_k(0) /)*ndays )
 oc_emi = no - toint( (/ oc_em(1), oc_em(0) /)*ndays )
 oc_wmi = toint( oc_wm*ndays )
 oc_r1i = toint( oc_r1*ndays )

 vv = new((/nw,nv,ny,nx/),float)
 cfx = new((/2,nw,ny,no,nx/),float)
 do iv=0, nv-1
   cf(:,iv,:,:,oc_pn_i(0):oc_pn_i(1),:) = 0.  ; high-freq
   cfe = cf(:,iv,:,:,:,:)
   cfw = cf(:,iv,:,:,:,:)
   cfe(:,:,:,:no/2,:) = 0.
   cfw(:,:,:,no/2:,:) = 0.

   cfk  = cfe(:,0,:,:,:)
   cfr1 = cfw(:,0,:,:,:)
   cfrge= cfe(:,1,:,:,:)
   cfrgw= cfw(:,1,:,:,:)

;cfk  (:,:,:oc_ki(0)-1,:) = 0.
;cfk  (:,:,oc_ki(1)+1:,:) = 0.
cfrge(:,:,:oc_emi(0)-1,:) = 0.
cfrge(:,:,oc_emi(1)+1:,:) = 0.
cfrgw(:,:,:oc_wmi(0)-1,:) = 0.
cfrgw(:,:,oc_wmi(1)+1:,:) = 0.
cfr1 (:,:,:oc_r1i(0)-1,:) = 0.
cfr1 (:,:,oc_r1i(1)+1:,:) = 0.
;cfrge(:,:,153:,:) = 0.
;cfrgw(:,:,:26,:) = 0.
;cfrgw(:,:,30:,:) = 0.
;   cfr2 = cfw(:,1,:,:,:)

   cfx = 0.
   do i=kc(0), kc(1)
     cfx(:,0,:,:,i) = cfftb(cfk  (:,:,:,i),0)
     cfx(:,1,:,:,i) = cfftb(cfr1 (:,:,:,i),0)
     cfx(:,2,:,:,i) = cfftb(cfrge(:,:,:,i),0)
     cfx(:,3,:,:,i) = cfftb(cfrgw(:,:,:,i),0)
   end do

   do j=0, 3
     vv(j,iv,:,:) = cfftb(cfx(:,j,:,n0(j),:),1)*2. * (nx*no)
   end do
 end do

 vv2 = new((/nw,nv,ny2,nx/),float)
 vv2(:,:,ny-1:,:) = vv(:,:,:,:)
 vv2(0:1,:,0:ny-2,:) = vv(0:1,:,ny-1:1,:)
 vv2(2:3,:,0:ny-2,:) = vv(2:3,:,ny-1:1,:)*(-1.)
 do iv=0, nv-1
   if (varname(iv) .eq. "rhov") then
     vv2(:,iv,0:ny-2,:) = vv2(:,iv,0:ny-2,:)*(-1.)
   end if
   if (varname(iv) .eq. "rhov" .or. varname(iv) .eq. "rhodz_dt") then
     do j=0, ny2-1
       vv2(:,iv,j,:) = vv2(:,iv,j,:)/rho0
     end do
   end if
 end do

 lon = fspan(0,360.-(360./nx),nx) - 180.

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlevh = fspan(-1.9,1.9,20)*2.
 cnfac = (/2.,0.5,0.5,0.2/)

;--------------------------------------- contour
 wi = 32
 he = 20
 xf = 8 + (/0,1,0,1/)*40
 yf = 10 + (/0,0,1,1/)*30

 mstr  = new(nw,string)
 xstr  = new(nw,string)
 ystr  = new(nw,string)
 lbopt = new(nw,string)

 mstr = (/"Kelvin","n=1 Rossby","E-MRG","W-MRG"/)
 xstr(:) = ""
 xstr(2:3) = "longitude (deg)"
 ystr(:) = ""
 ystr((/0,2/)) = "latitude (deg)"
 lbtxt = (/"","","","(mm s:S:-1:N:)"/)

do il=0, nw-1

 do iv=0, nv-1
   if (varname(iv) .eq. "rhodz_dt") then
     contour1 = cnshadeplot(wks,lon,lat,vv2(il,iv,:,:)*1.e3,cnlev*cnfac(il),cnfil,lbtxt(il))
   end if
 end do
 do iv=0, nv-1
   if (varname(iv) .eq. "theta") then
     cnaddline(contour1,wks,lon,lat,vv2(il,iv,:,:),cnlevh*cnfac(il),0,"D",0,"T",0.5,"C",(/colindr,1,colindb/))
   end if
 end do
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.15
  "pmLabelBarParallelPosF"   : 0.50-0.40
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.0
  "pmLabelBarHeightF"        : 0.15
  "lbLabelFontHeightF" : 0.025
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
  "lbTitleFontHeightF" : 0.030
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-180/kc(0),180/kc(0),-22,22,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))

 cnfinalize(contour1,wks)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vv2(il,0:1,2::3,:)
  "vfXArray"    : lon
  "vfYArray"    : lat(2::3)
 end create
 ref = cnfac*1.
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref(il)
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(il),yf(il),wi,he)
 axis_range(vec1,-180/kc(0),180/kc(0),-22,22,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

end do

 frame(wks)

end

