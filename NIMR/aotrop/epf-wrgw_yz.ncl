load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "tc"
 varname = (/"epd","f_y","f_z"/)

 y1 = 1951
 y2 = 1970
 lat1 = 25.
 z1 = 16.
 z2 = 44.

 k0 = 0

 z0 = (/21.0,23.3,26.9/)  ; 50, 30, 20 hPa
 z = z0(k0)

 if (k0 .eq. 0) then
   udist = fspan(-19.5,10.5,11)
 end if
 if (k0 .eq. 1) then
   udist = fspan(-26,14,11)
 end if
 if (k0 .eq. 2) then
   udist = fspan(-30.5,19.5,11)
 end if
; udist = fspan(-33,27,11)
; udist = fspan(-40,32,10)
 n_int = 10 ; 3-day interpol.

 wks = cr_wks("fig/epf-wrgw_yz"+numeric2int(z,1)+"-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2 - y1 + 1)*12 - 1
 nv = dimsizes(varname) + 1

 fname = "$DATD/AOL60CGW/xm.uan"+expname+"a.pj.nc"
 f1 = addfile(fname,"r")
 lat = f1->latitude({-lat1:lat1})
 ht = f1->hybrid_ht({z1*1.e3:z2*1.e3})
 ny = dimsizes(lat)
 nz = dimsizes(ht )
 kk = closest_val(z*1.e3,ht)
 if (ht(kk) .ne. ht({z*1.e3})) then
   print("wrong index")
 end if
 u0 = dim_avg( f1->u(:nt-1,{z*1.e3},{-5:5}) )
 u  = f1->u(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("$DATD/AOL60CGW/tem/uan"+expname+".tem_yzt.nc","r")
 ut0 = dim_avg( f2->u_tend(:nt-1,{z*1.e3},{-5:5}) )
 ut  = f2->u_tend(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 f2a = addfile("$DATD/AOL60CGW/epf/uan"+expname+".epfsa0_yz_recon0.nc","r")
 f2b = addfile("$DATD/AOL60CGW/epf/uan"+expname+".epf_yz_recon0.nc","r")
 wv = new((/nv-1,nt,nz,ny/),float)
 do iv=0, nv-2
   wv(iv,:,:,:) = f2a->$varname(iv)$(:nt-1,5,{z1*1.e3:z2*1.e3},{-lat1:lat1})
 end do
; wv(1:2,:,:,:) = wv(1:2,:,:,:)*(-1.)  ; for eastward waves

; kw   = f2a->epd(:nt-1,2,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; rw   = f2a->epd(:nt-1,3,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; ergw = f2a->epd(:nt-1,4,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; wrgw = f2a->epd(:nt-1,5,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; kas  = f2a->epd(:nt-1,6,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; eigw = f2b->epd(:nt-1,8,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; wigw = f2b->epd(:nt-1,9,{z1*1.e3:z2*1.e3},{-lat1:lat1})
; wtw  = f2b->epd(:nt-1,2,{z1*1.e3:z2*1.e3},{-lat1:lat1}) + \
;        f2b->epd(:nt-1,3,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 f3 = addfile("$DATD/AOL60CGW/xm.uan"+expname+"a.pe.nc","r")
 gw = f3->du_gw(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})*72.
 if (cexp .ne. "CTL") then
   gwc = f3->dudt_gwconv(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})*86400.
   gw  = gw + gwc
 end if

print(max(u0)+"  "+min(u0))

;--------------------------------------- process 0
 ndist = dimsizes(udist)
 itv = udist(1) - udist(0)
 uitv = (/udist-0.5*itv,udist+0.5*itv/)

 nti = (nt-1)*n_int+1
 ti  = fspan(t(0),t(nt-1),nti)

 ui0  = linint1(t,u0 ,False,ti,0)
 uti0 = linint1(t,ut0,False,ti,0)
 ui   = linint1_n(t,u ,False,ti,0,0)
 uti  = linint1_n(t,ut,False,ti,0,0)
 wvi  = linint1_n(t,wv,False,ti,0,1)
 gwi  = linint1_n(t,gw,False,ti,0,0)

 vari = new((/nv,nti,nz,ny/),float)
 vari(0,:,:,:) = ui
 vari(1:,:,:,:) = wvi

; determine criterion vari.
 cri = new((/2,nti/),float)

 zeropt = new(nti,integer)
 zeropt = -999

 tag = 0
 nlast = -999
 nlast2 = -999
 cnt = 0
 do n=1, nti-1
   if ( ui0(n)*ui0(n-1) .le. 0. ) then 
     zeropt(cnt) = n
     if (n-nlast2 .lt. 6*n_int) then
       zeropt(cnt-1:cnt) = -999
       cnt = cnt - 1
     else
       cnt = cnt + 1
     end if
     nlast2 = nlast
     nlast = n
   end if
 end do
 nzero = cnt

 maxmin = new((/2,nzero-1/),integer)
 do cnt=0, nzero-2
   if ( sum(ui0(zeropt(cnt):zeropt(cnt+1))) .gt. 0. ) then
     maxmin(0,cnt) = maxind(ui0(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = 1 
   else
     maxmin(0,cnt) = minind(ui0(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = -1
   end if
 end do

 cri = 0.
 do cnt=0, nzero-3
   cri(0,maxmin(0,cnt):maxmin(0,cnt+1)) = tofloat(maxmin(1,cnt))*(-1.0)
   ; test
   if (maxmin(1,cnt)*maxmin(1,cnt+1) .ge. 0) then
     print("Check here !!")
   end if
 end do
 ncyc = new(2,float)
 ncyc(:) = tofloat((nzero-2)/2)
 if (mod(nzero-2,2) .ne. 0) then
   if (maxmin(1,0) .lt. 0) then
     ncyc(0) = ncyc(0) + 1.0
   else
     ncyc(1) = ncyc(1) + 1.0
   end if
 end if
print(ncyc)

 vard = new((/2,ndist,nv,nz,ny/),float)
 vard = 0.
 ncom = new((/2,ndist/),float)
 ncom = 0.
 do n=0, nti-1
   if (cri(0,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ui0(n) .ge. uitv(0,ii) .and. ui0(n) .lt. uitv(1,ii) ) then
         ncom(0,ii) = ncom(0,ii) + 1.
         vard(0,ii,:,:,:) = vard(0,ii,:,:,:) + vari(:,n,:,:)
       end if
     end do
   end if
   if (cri(0,n) .lt. 0.0) then
     do ii=0, ndist-1
       if ( ui0(n) .ge. uitv(0,ii) .and. ui0(n) .lt. uitv(1,ii) ) then
         ncom(1,ii) = ncom(1,ii) + 1.
         vard(1,ii,:,:,:) = vard(1,ii,:,:,:) + vari(:,n,:,:)
       end if
     end do
   end if
 end do
 ncom@_FillValue = 0.

 do ip=0, 1
 do ii=0, ndist-1
 do iv=0, nv-1
   vard(ip,ii,iv,:,:) = vard(ip,ii,iv,:,:)/ncom(ip,ii)
 end do
 end do
 end do

print(ncom)

;--------------------------------------- cn level
 cnlev = fspan(-20,20,21)*0.01
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 nfig = 2
 nu   = 5
 nstrd= 2
 ustt = 1

 lat1 = 25
 z1p = 17
 z2p = 42

 wi = 24
 he = 36
 xf = 15 + (/0,1,2,3,4/)*32
 yf = 10 + (/0,1/)*45

 r_yz = 3.141592/180.*6370.
 vard(:,:,3,:,:) = vard(:,:,3,:,:)*r_yz
 do k=0, nz-1
   vard(:,:,2:3,k,:) = vard(:,:,2:3,k,:)*exp(ht(k)/(2.*7.e3))
 end do

 mstr  = (/"Easterly to westerly","Westerly to easterly"/)
 xstr  = (/"","lat (deg)"/)
 ystr  = (/"z (km)","","","",""/)

do ifig=0, nfig-1
do iu=0, nu-1

 iu1 = iu*nstrd + ustt
 iu2 = (/iu1,ndist-1-iu1/)
 contour1 = cnshadeplot(wks,lat,ht/1.e3,vard(ifig,iu2(ifig),1,:,:),cnlev,cnfil,0)
 contour0 = cnlineplot(wks,lat,ht/1.e3,vard(ifig,iu2(ifig),0,:,:),ispan(-50,50,5),0)
 cnline_att(contour0,"D",(/0,0,16/),"T",(/0.8,2.0,1.0/),"C",1)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vard(ifig,iu2(ifig),2:3,::2,::4)
  "vfXArray"    : lat(::4)
  "vfYArray"    : ht(::2)/1.e3
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : 1.e6
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
  "vcGlyphStyle" : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
 end create

 set_size(contour1,xf(iu),yf(ifig),wi,he)
 axis_range(contour1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour1,"",xstr(ifig),ystr(iu))

 set_size(contour0,xf(iu),yf(ifig),wi,he)
 axis_range(contour0,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour0,"","","")

 set_size(vec1,xf(iu),yf(ifig),wi,he)
 axis_range(vec1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)
 cnfinalize(contour0,wks)

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

; draw(contour1)
; draw(contour0)

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

 draw(contour0)

end do
end do

 frame(wks)

end

