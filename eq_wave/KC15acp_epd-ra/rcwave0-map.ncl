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

 yy = "1958";"1983"
 mm = "02";"05"
 ndays = 90
 z0 = 21. ;18.59 ;20.95 ;17.41 ;20.95
; n0 = (/115,115,60,60/)   ; 60 - 120
 lat0 = 26
 wv = (/ (/"k",""/), (/"","mrg"/), (/"r_s","r_a"/), (/"ig_s","ig_a"/) /)

; kc = (/3,3/)
 kc = (/1,4/)
 nk = 20
 no = 1*(2*ndays)  ; 1-day  ; cut-off period

 wks = cr_wks("fig/rcwave0-map",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)
 colindg = NhlNewColor(wks,0.2,0.9,0.2)
 colind0 = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = dimsizes(varname)
 nw = 4

 fdir  = "$DATD/L60CGW-t/fcoef_ko-wav4"
 fname = fdir+"/"+yy+"/uan"+expname+".fft_var4_ko."+yy+"."+mm+".nc"
print(fname+"")

; axis
 f0 = addfile(fname(0),"r")
 s_v0 = "fcr_"+varname(0)+"_k"
 lath = f0->$(f0->$s_v0$!1)$({0:lat0})
 lat = array_append_record(-lath({lat0:0.1}), lath, 0)
 ny2 = dimsizes(lat)
 ny = dimsizes(lath)
 nx = nk*2
 noa = dimsizes(f0->$(f0->$s_v0$!2)$)

 f1 = addfile("$DATD/L60CGW-t/xm/xm.uan"+expname+"a.pj_"+yy+mm+".nc","r")
 rhob = f1->rho(0,{z0*1.e3},{-lat0:lat0})
 ub   = f1->u  (0,{z0*1.e3},{-lat0:lat0})
print(ub)

 cfs = new((/2,nw,nv,ny,no,nk+1/),float)
 cfa = new((/2,nw,nv,ny,no,nk+1/),float)
 cfs = 0.
 cfa = 0.
 si = (/"r","i"/)
 do iw=0, nw-1
 do iv=0, nv-1
   f0 = addfile(fname,"r")
   do i=0, 1
     if (wv(iw,0) .ne. "") then
       s_v0 = "fc"+si(i)+"_"+varname(iv)+"_"+wv(iw,0)
       cfs(i,iw,iv,:,:no/2,:) = f0->$s_v0$({z0*1.e3},{0:lat0},:no/2,:nk)
       cfs(i,iw,iv,:,no/2+1:,:) = f0->$s_v0$({z0*1.e3},{0:lat0},noa-no/2+1:,:nk)
     end if
     if (wv(iw,1) .ne. "") then
       s_v0 = "fc"+si(i)+"_"+varname(iv)+"_"+wv(iw,1)
       cfa(i,iw,iv,:,:no/2,:) = f0->$s_v0$({z0*1.e3},{0:lat0},:no/2,:nk)
       cfa(i,iw,iv,:,no/2+1:,:) = f0->$s_v0$({z0*1.e3},{0:lat0},noa-no/2+1:,:nk)
     end if
   end do
 end do
 end do
 cfs(:,:,:,:,no/2,:) = 0.
 cfa(:,:,:,:,no/2,:) = 0.

;--------------------------------------- process 0
 rho0 = avg(rhob({-10:10}))

 nt = no/3
 vvs0 = new((/nw,nv,nt,ny,nx/),float)
 vva0 = new((/nw,nv,nt,ny,nx/),float)
 vvs0 = 0.
 vva0 = 0.
 cfsx = new((/2,nw,ny,no,nx/),float)
 cfax = new((/2,nw,ny,no,nx/),float)
 cfsx = 0.
 cfax = 0.
 do iv=0, nv-1

;cfs(:,iv,:,no/2:,:) = 0.
;cfa(:,iv,:,no/2:,:) = 0.

   cfsx = 0.
   cfax = 0.
   do i=kc(0), kc(1)
     cfsx(:,:,:,:,i) = cfftb(cfs(:,:,iv,:,:,i),0)
     cfax(:,:,:,:,i) = cfftb(cfa(:,:,iv,:,:,i),0)
   end do

   do n=0, nt-1
   do j=0, nw-1
     vvs0(j,iv,n,:,:) = cfftb(cfsx(:,j,:,n+no/3,:),1)*2. * (nx*no)
     vva0(j,iv,n,:,:) = cfftb(cfax(:,j,:,n+no/3,:),1)*2. * (nx*no)
   end do
   end do

 end do

 vv_xy = new((/nw,nv,ny2,nx*2/),float)
 vv2 = new((/nw,nv,nt,ny2,nx/),float)
 vv2(:,:,:,ny-1:,:) = vvs0(:,:,:,:,:) + vva0(:,:,:,:,:)
 vv2(:,:,:,0:ny-2,:) = vvs0(:,:,:,ny-1:1,:) - vva0(:,:,:,ny-1:1,:)
 do iv=0, nv-1
   if (varname(iv) .eq. "rhov") then
     vv2(:,iv,:,0:ny-2,:) = vv2(:,iv,:,0:ny-2,:)*(-1.)
   end if
   if (varname(iv) .eq. "rhov" .or. varname(iv) .eq. "rhodz_dt") then
     do j=0, ny2-1
       vv2(:,iv,:,j,:) = vv2(:,iv,:,j,:)/rho0
     end do
   end if
 end do

 detr = abs(vv2(0,3,:,ind(lat .eq. -15):ind(lat .eq. 15),:))

 ndetr = dimsizes(detr)
 detr1d = ndtooned(detr)
 inds = ind_resolve(maxind(detr1d), ndetr)
 nn = inds(0,0)
 jj = inds(0,1)
 ii = inds(0,2)    

;nn = 58
print(nn)

 lon = fspan(0,720.-(360./nx),nx*2) - 360.
 lon!0 = "longitude"
 lon&longitude = lon

 vv_xy(:,:,:,:nx-1) = vv2(:,:,nn,:,:)
 vv_xy(:,:,:,nx:  ) = vv2(:,:,nn,:,:)
 vv_xy!3 = "longitude"
 vv_xy&longitude = lon

; tmp = vv2(:,2,:,:)
; do j=1, ny2-2
;   vv2(:,2,j,:) = (2.*tmp(:,j,:) + tmp(:,j+1,:) + tmp(:,j-1,:))/4.
; end do

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)
 cnlevu = ispan(-10,10,1)*1.5
 cnlevv = ispan(-10,10,1)*0.5
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
; cnfil = (/22,23,24,25,26,0,0,11,10,9,8,7/)
 cnlevh = ispan(-10,10,1)*0.5
 cnfac = (/1.2,0.25,1.2,0.25/)

;--------------------------------------- contour
 wi = 90
 he = 25 * 0.75
 xf = 8 + (/0,0,0,0/)*52
 yf = 10 + (/0,2,1,3/)*(he+3) + (/0,1,0,1/)*8

 mstr  = new(nw,string)
 xstr  = new(nw,string)
 ystr  = new(nw,string)
 lbopt = new(nw,string)

 mstr = "";"Rossby"
 xstr = (/"","","longitude (deg)","longitude (deg)"/)
 ystr = (/"latitude (deg)","","latitude (deg)",""/)
 xstr = ""
 ystr = ""
 lbtxt = "(K)"

do il=0, nw-1

 map = create "mp" mapPlotClass wks
;  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 0.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"

  "mpShapeMode"     : "FreeAspect"
  "mpGeophysicalLineColor" : colindg ;"Orange4"
  "mpGeophysicalLineThicknessF" : 2.5
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : -180
  "mpMaxLonF"      : 180
  "mpMinLatF"      : -25
  "mpMaxLatF"      : 25
 end create
 set_size_mp(map,xf(il),yf(il),wi,he)

; do iv=0, nv-1
iv=3
   if (varname(iv) .eq. "theta") then
     contour1 = cnshadeplot(wks,lon({-180:180}),lat,vv_xy(il,iv,:,{-180:180}), \
                cnlevh*cnfac(il),cnfil,lbtxt)
;     cnaddline(contour1,wks,lon({-180:180}),lat,vv_xy(il,iv,:,{-180:180}), \
;               cnlevh*cnfac(il),0,"D",0,"T",0.5,"C",(/colindr,1,colindb/))
   end if
; end do

; do iv=0, nv-1
;   if (varname(iv) .eq. "rhodz_dt") then
;;     contour1 = cnshadeplot(wks,lon({-180:180}),lat,vv_xy(il,iv,:,{-180:180})*1.e3, \
;;                cnlevh*cnfac(il),cnfil,lbtxt)
;     cnaddline(contour1,wks,lon({-180:180}),lat,vv_xy(il,iv,:,{-180:180})*1.e3, \
;               cnlevh*cnfac(il),0,"D",0,"T",0.5,"C",(/colindr,-1,colindb/))
;   end if
; end do

 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.15-0.12
  "pmLabelBarParallelPosF"   : 0.50;-0.40
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.08
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10-0.10
  "lbTitleFontHeightF" : 0.014
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-180/kc(0),180/kc(0),-25,25,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))

; cnfinalize(contour1,wks)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vv_xy(il,0:1,::4,{-180:180})
  "vfXArray"    : lon({-180:180})
  "vfYArray"    : lat(::4)
 end create
 ref = cnfac*2.
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref(il)
  "vcRefLengthF"    : 0.010
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(il),yf(il),wi,he)
 axis_range(vec1,-180/kc(0),180/kc(0),-25,25,False,False)
 axis_str(vec1,"","","")

; overlay(contour1,vec1)
; draw(contour1)
; NhlRemoveOverlay(contour1,vec1,False)

 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)
 draw(map)
 NhlRemoveOverlay(map,contour1,False)
 NhlRemoveOverlay(map,vec1    ,False)

end do

 frame(wks)

end

