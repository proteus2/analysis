load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "tc"
 varname = (/"epd","f_y","f_z"/)

 y1 = 1953
 y2 = 2006
 lat1 = 25.
 z1p = 17
 z2p = 39
 z1 = z1p - 1.5

 k0 = 1

 z0 = (/21.0,23.3,26.9/)  ; 50, 30, 20 hPa
 z = z0(k0)

; iwvs = (/7,8/)
; swvs = "cgw"
 iwvs = (/0,5/)
 swvs = "lw"
 ; 0 - 6   : K, R, E-MRG, W-MRG, E-IG, W-IG, mig.-tide
 ; 7 - 10  : E-CGW, W-CGW, E-BGW, W-BGW
 ; 11 - 12 : E-GW, W-GW
 swv = (/"k","r","emrg","wmrg","eig","wig","tid", \
         "ecgw","wcgw","ebgw","wbgw","egw","wgw"/)
 scl = (/1.,2.,0.1,0.2,1.,1.,1., \
         2.,2.,2.,2.,3.,3./)

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

 wks = cr_wks("fig/epf-"+swvs+"_yz"+numeric2int(z,1)+"-m3-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

do iw=0, 1
  iwv = iwvs(iw)
;--------------------------------------- read data
 nt = (y2 - y1 + 1)*12 - 1
 nv = dimsizes(varname) + 1

 z2 = z2p + 4

 fname = "$DATD/AOL60CGW/xm.uan"+expname+"a.pj_53-06.nc"
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

 f2 = addfile("$DATD/AOL60CGW/tem/uan"+expname+".tem-wc_yzt_53-06.nc","r")
 ut0 = dim_avg( f2->u_tend(:nt-1,{z*1.e3},{-5:5}) )
 ut  = f2->u_tend(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 f2a = addfile("$DATD/AOL60CGW/epf-wc/uan"+expname+".epfsa_yz_recon0_53-06.nc","r")
 f2b = addfile("$DATD/AOL60CGW/epf-wc/uan"+expname+".epf_yz_recon0_53-06.nc","r")

 if (iw .eq. 0) then
   wv  = new((/nv-1,nt,nz,ny/),float)
   wv2 = new((/nv-1,nt,nz,ny/),float)
 end if
 ; 0 - 6   : K, R, E-MRG, W-MRG, E-IG, W-IG, mig.-tide, MRG
 ; 7 - 10  : E-CGW, W-CGW, E-BGW, W-BGW
 ; 11 - 12 : E-GW, W-GW

 do iv=0, nv-2
   if (iwv .le. 5) then
     wv2(iv,:,:,:) = f2a->$varname(iv)$(:nt-1,iwv+2, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1})
     if (iwv .eq. 1) then
       wv2(iv,:,:,:) = wv2(iv,:,:,:) + f2a->$varname(iv)$(:nt-1,iwv+12, \
          {z1*1.e3:z2*1.e3},{-lat1:lat1})
     end if
     if (iwv .ge. 4) then
       wv2(iv,:,:,:) = wv2(iv,:,:,:) + f2b->$varname(iv)$(:nt-1,iwv+4, \
          {z1*1.e3:z2*1.e3},{-lat1:lat1}) \
                                   + f2a->$varname(iv)$(:nt-1,iwv+12, \
          {z1*1.e3:z2*1.e3},{-lat1:lat1})
     end if
   end if
   if (iwv .eq. 6) then
     wv2(iv,:,:,:) = f2b->$varname(iv)$(:nt-1,2, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1}) + \
                    f2b->$varname(iv)$(:nt-1,3, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1})
   end if
 end do
 if ( iwv .le. 5 .and. mod(iwv,2) .eq. 0 ) then
   wv2(1:2,:,:,:) = wv2(1:2,:,:,:)*(-1.)  ; for eastward waves
 end if

 if (iwv .ge. 7) then
   f0 = addfile("$DATD/AOL60CGW/xm.uan"+expname+"a.pj_53-06.nc","r")
   rho0 = f0->rho(:nt-1,{z1*1.e3:z2*1.e3},{-lat1:lat1})
   if (iwv .eq. 7 .or. iwv .eq. 11) then
     f3 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pd_53-06.nc","r")
     varn = "F_e_gwconv"
   end if
   if (iwv .eq. 8 .or. iwv .eq. 12) then
     f3 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pd_53-06.nc","r")
     varn = "F_w_gwconv"
   end if
   if (iwv .eq. 9) then
     f3 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pg_53-06.nc","r")
     varn = "Fp_e_gwussp"
   end if
   if (iwv .eq. 10) then
     f3 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pg_53-06.nc","r")
     varn = "Fp_w_gwussp"
   end if
   mf = f3->$varn$(:nt-1,{z1*1.e3-600.:z2*1.e3+1000.},{-lat1:lat1})
   ht1= mf&$mf!1$
   if (iwv .eq. 10) then
     mf = mf*(-1.)
   end if
   if (iwv .eq. 11) then
     f4 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pg_53-06.nc","r")
     mf = mf + f4->Fp_e_gwussp(:nt-1,{z1*1.e3-600.:z2*1.e3+1000.},{-lat1:lat1})
   end if
   if (iwv .eq. 12) then
     f4 = addfile("$DATD/AOL60CGW/mmxm.uan"+expname+"a.pg_53-06.nc","r")
     mf = mf - f4->Fp_w_gwussp(:nt-1,{z1*1.e3-600.:z2*1.e3+1000.},{-lat1:lat1})
   end if
   do k=0, nz-1
     wv2(2,:,k,:) = 0.5*(mf(:,k,:) + mf(:,k+1,:))
     wv2(0,:,k,:) = -(mf(:,k+1,:) - mf(:,k,:))/(ht1(k+1)-ht1(k))/rho0(:,k,:)
   end do
   wv2(0,:,:,:) = wv2(0,:,:,:)*86400.
   wv2(1,:,:,:) = 0.
   wv2(2,:,:,:) = wv2(2,:,:,:)*6371.e3
   if ( mod(iwv,2) .eq. 0 ) then
     wv2(2,:,:,:) = wv2(2,:,:,:)*(-1.)
   end if
   do j=0, ny-1
     wv2(1:2,:,:,j) = wv2(1:2,:,:,j)*cos(lat(j)*3.141592/180.)
   end do
 end if

 if (iw .eq. 0) then
   wv = wv2
   wv2 = 0.
 end if

end do  ; iw
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
 wv2i  = linint1_n(t,wv2,False,ti,0,1)

 vari = new((/nv,nti,nz,ny/),float)
 vari(0,:,:,:) = ui
 vari(1:,:,:,:) = wvi

 vari2 = new((/nv,nti,nz,ny/),float)
 vari2(0,:,:,:) = ui
 vari2(1:,:,:,:) = wv2i

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

 tim = ti - tofloat(toint(ti))

 vard = new((/2,ndist,nv,nz,ny/),float)
 vard = 0.
 vard2 = vard
 ncom = new((/2,ndist/),float)
 ncom = 0.
 do n=0, nti-1
if ( tim(n) .gt. 11./12. .or. tim(n) .lt. 2./12. ) then  ; DJF
;if ( tim(n) .gt. 5./12. .and. tim(n) .lt. 8./12. ) then  ; JJA
   if (cri(0,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ui0(n) .ge. uitv(0,ii) .and. ui0(n) .lt. uitv(1,ii) ) then
         ncom(0,ii) = ncom(0,ii) + 1.
         vard (0,ii,:,:,:) = vard (0,ii,:,:,:) + vari (:,n,:,:)
         vard2(0,ii,:,:,:) = vard2(0,ii,:,:,:) + vari2(:,n,:,:)
       end if
     end do
   end if
   if (cri(0,n) .lt. 0.0) then
     do ii=0, ndist-1
       if ( ui0(n) .ge. uitv(0,ii) .and. ui0(n) .lt. uitv(1,ii) ) then
         ncom(1,ii) = ncom(1,ii) + 1.
         vard (1,ii,:,:,:) = vard (1,ii,:,:,:) + vari (:,n,:,:)
         vard2(1,ii,:,:,:) = vard2(1,ii,:,:,:) + vari2(:,n,:,:)
       end if
     end do
   end if
end if
 end do
 ncom@_FillValue = 0.

 do ip=0, 1
 do ii=0, ndist-1
 do iv=0, nv-1
   vard (ip,ii,iv,:,:) = vard (ip,ii,iv,:,:)/ncom(ip,ii)
   vard2(ip,ii,iv,:,:) = vard2(ip,ii,iv,:,:)/ncom(ip,ii)
 end do
 end do
 end do

print(ncom)

 vard(:,:,1,:,:) = vard(:,:,1,:,:)*30.  ; m/s/month
 vard(:,:,2:3,:,:) = vard(:,:,2:3,:,:) / 6371.e3 * 1.e3  ; mPa
 vard2(:,:,1,:,:) = vard2(:,:,1,:,:)*30.  ; m/s/month
 vard2(:,:,2:3,:,:) = vard2(:,:,2:3,:,:) / 6371.e3 * 1.e3  ; mPa

;--------------------------------------- cn level
 cnlev = fspan(-5,5,21)*scl(0)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 nfig = 2
 nu   = 5
 nstrd= 2
 ustt = 1

 wi = 24
 he = 36
 xf = 15 + (/0,1,2,3,4/)*30
 yf = 10 + (/0,1,2,3/)*42 + (/0,0,1,1/)*6

 r_yz = 3.141592/180.*6371.
 vard (:,:,2,:,:) = vard (:,:,2,:,:)/r_yz
 vard2(:,:,2,:,:) = vard2(:,:,2,:,:)/r_yz
 do k=0, nz-1
   vard (:,:,2:3,k,:) = vard (:,:,2:3,k,:)*exp(ht(k)/(2.*6.6e3))
   vard2(:,:,2:3,k,:) = vard2(:,:,2:3,k,:)*exp(ht(k)/(2.*6.6e3))
 end do

do iw=0, 1
  iwv = iwvs(iw)

 mstr  = (/"Easterly to westerly","Westerly to easterly"/)
 xstr  = (/"","lat (deg)"/)
 ystr  = (/"z (km)","","","",""/)
 refann = new((/nfig,nu/),logical,False)
 if (iw .eq. 1) then
   refann(nfig-1,nu-1) = True
 end if

  vvv = vard
  if (iw .eq. 1) then
    vvv = vard2
  end if

do ifig=0, nfig-1
do iu=0, nu-1

 iu1 = iu*nstrd + ustt
 iu2 = (/iu1,ndist-1-iu1/)
 contour1 = cnshadeplot(wks,lat,ht/1.e3,vvv(ifig,iu2(ifig),1,:,:),cnlev,cnfil,0)
 contour0 = cnlineplot(wks,lat,ht/1.e3,vvv(ifig,iu2(ifig),0,:,:),ispan(-50,50,5),0)
 cnline_att(contour0,"D",(/0,0,16/),"T",(/0.6,1.8,0.7/),"C",(/1,colind,1/))

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vvv(ifig,iu2(ifig),2:3,::2,::4)
  "vfXArray"    : lat(::4)
  "vfYArray"    : ht(::2)/1.e3
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
;  "vcLineArrowThicknessF" :  vthc
;  "vcLineArrowColor"      :  vcol
  "vcRefMagnitudeF" : scl(iwv)
  "vcRefLengthF"    : 0.05
  "vcGlyphStyle"              : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.12
  "vcFillArrowWidthF"         : 0.1
  "vcPositionMode"            : "ArrowTail"
  "vcRefAnnoOn"          : refann(ifig,iu)
  "vcRefAnnoArrowAngleF" : 90.
  "vcRefAnnoFontHeightF" : 0.024
  "vcRefAnnoString1On"   : False
;  "vcRefAnnoString1"     : "1 mPa"
;  "vcRefAnnoString2On"   : False
  "vcRefAnnoString2"     : scl(iwv)+" mPa"
  "vcRefAnnoPerimSpaceF" : 0.6
  "vcRefAnnoOrthogonalPosF" : 0.15
;  "vcRefAnnoParallelPosF"   :  1. - nn*0.1
 end create

 set_size(contour1,xf(iu),yf(ifig+iw*2),wi,he)
 axis_range(contour1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour1,"",xstr(ifig),ystr(iu))

 set_size(contour0,xf(iu),yf(ifig+iw*2),wi,he)
 axis_range(contour0,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour0,"","","")

 set_size(vec1,xf(iu),yf(ifig+iw*2),wi,he)
 axis_range(vec1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)
 cnfinalize(contour0,wks)

 draw(contour1)

 overlay(contour0,vec1)
 draw(contour0)
 NhlRemoveOverlay(contour0,vec1,False)

end do
end do

end do  ; iw

 label = create "lab" labelBarClass wks
  "lbOrientation" : "horizontal"
  "lbBoxCount" : dimsizes(cnfil)
  "lbFillColors" : cnfil
  "lbMonoFillPattern" : True
  "lbLabelsOn" : True
  "lbLabelStrings" : cnlev
  "lbAutoManage" : False
  "lbLabelAutoStride" : True
;  "lbLabelPosition" : "Bottom" ;"Right"
  "lbLabelAlignment" : "InteriorEdges" ;"BoxCenters" ;"ExternalEdges"
;  "lbLabelAngleF" : 0.
  "lbPerimOn" : False
  "lbTitleString" : "(m s:S:-1:N: month:S:-1:N:)"
  "lbTitlePosition" : "Right" ;"Top"
  "lbTitleDirection" : "Across"
 end create
 set_size_lb(label,xf(0),yf(3)+he+10,90,4)

 draw(label)

 dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
 dr_txt(wks,"(b)",xf(0)-8,yf(2)-2,25,8,0.)

 frame(wks)

end

