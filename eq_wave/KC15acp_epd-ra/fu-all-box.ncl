load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 0.85

begin

;--------------------------------------- parameters
 expname = (/"era-int","merra","jra55","cfsr","era-int_f"/)
 datname = (/"ERA-I"  ,"MERRA","JRA55","CFSR","ERA-I-nr"/)
 uname = (/"U_GDS0_ISBL","u","U_GRD_GDS0_ISBL","UGRD","u"/)

 y0 = 1979
 feb_1st = True
 y1 = 1979
 y2 = 2010
 lat0 = 5
 p = 30
 p_em = 19  ; (50,30,20,10) : (22,19,17,14)

 wks = cr_wks("fig/fu-all-box-p"+p,600,0)
; gsn_define_colormap(wks,"amwg")
; gsn_define_colormap(wks,"cosam")
; gsn_define_colormap(wks,"radar")
 colind = NhlNewColor(wks,0.70,0.70,0.70)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nw = 8
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12
 ntr = nt
 if ( feb_1st .and. ( n0 .eq. 0 ) ) then
   ntr = nt - 1
 end if

 fnam0 = "$DATD/clim/"+datname+"/"+expname+".u_ypt.nc"
 fnam1 = "$DATD/tem/"+datname+"/epf/"+expname+".epfsa_yz_recon0.nc"
 fnam2 = "$DATD/tem/"+datname+"/epf/"+expname+".epf_yz_recon0.nc"
 fnam3 = "$DATD/tem/"+datname+"/"+expname+".tem_ypt.nc"

 latrng = (/-lat0,lat0/)

 do ie=0, ne-1
   f0 = addfile(fnam0(ie),"r")
   f1 = addfile(fnam1(ie),"r")
   f2 = addfile(fnam2(ie),"r")
   f3 = addfile(fnam3(ie),"r")

   lat = f1->lat({-lat0:lat0})
   dy05 = 0.5*(lat(1) - lat(0))
   wgt = where( lat .ge. latrng(0)+dy05 .and. lat .le. latrng(1)-dy05, 1., 0.)
   wgt = where( abs(lat-latrng(0)) .lt. dy05, 0.5*(1.+(lat-latrng(0))/dy05), wgt)
   wgt = where( abs(lat-latrng(1)) .lt. dy05, 0.5*(1.-(lat-latrng(1))/dy05), wgt)

   print(expname(ie)+"  "+f1->p({p})+"  "+abs(lat(0))+"  "+wgt(0))

   rw_i = new((/ntr,nw/),float)
   rw_i(:,:) = dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,:nw-1,{p},{-lat0:lat0}), wgt,1)
   rw_i(:,nw-2:) = rw_i(:,nw-2:) + \
                   dim_avg_wgt_Wrap( f2->epd(n0:n0+ntr-1,8:9,{p},{-lat0:lat0}), wgt,1)
; migrating tide
   rw_i(:,nw-1) = rw_i(:,nw-1) + \
                  ( dim_avg_wgt_Wrap( f2->epd(n0:n0+ntr-1,2,{p},{-lat0:lat0}) + \
                            f2->epd(n0:n0+ntr-1,3,{p},{-lat0:lat0}) , wgt,1) )

   rw_i(:,3  ) = rw_i(:,3  ) + \
                 dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,13   ,{p},{-lat0:lat0}), wgt,1)
   rw_i(:,6:7) = rw_i(:,6:7) + \
                 dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,16:17,{p},{-lat0:lat0}), wgt,1)

   rw0 = new((/nw,nt/),float,1.e20)
   rw0(:,nt-ntr:nt-1) = transpose(rw_i)

   if (expname(ie) .eq. "era-int_f") then
     u0 = dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{p_em},{-lat0:lat0}) ), wgt,1)
   else
     u0 = dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{p},{-lat0:lat0}) ), wgt,1)
   end if
   if (expname(ie) .eq. "era-int") then
     uc = u0
   end if

   ut   = dim_avg_wgt_Wrap( f3->u_tend (n0:n0+nt-1,{p},{-lat0:lat0}), wgt,1)
   advz = dim_avg_wgt_Wrap( f3->uadv_z (n0:n0+nt-1,{p},{-lat0:lat0}), wgt,1)
   uf   = dim_avg_wgt_Wrap( f3->u_force(n0:n0+nt-1,{p},{-lat0:lat0}), wgt,1)
   ue   = dim_avg_wgt_Wrap( f3->epd    (n0:n0+nt-1,{p},{-lat0:lat0}), wgt,1)
   advy = dim_avg_wgt_Wrap( f3->uadv_y (n0:n0+nt-1,{p},{-lat0:lat0}) + \
                            f3->cor    (n0:n0+nt-1,{p},{-lat0:lat0}) , wgt,1)

;   ureq = (/ ut - advz - advy /)

   if (ie .eq. 0) then
     un = new((/9,ne,nt/),float)
     fu = new((/9,ne,nt/),float)
     t = y1 + (fspan(1,nt,nt)-0.5)/12.
;     tt = y1 + (fspan(1,nt-1,nt-1))/12.
   end if
   un(:,ie,:) = conform_dims((/9,nt/),ut,1)
   fu(:,ie,:) = (/ rw0(2,:), rw0(3,:), rw0(4,:)+rw0(5,:), \
                   rw0(6,:)+rw0(7,:), ue, advz, uf, ut, ut - uf /)
   iadvz = 5
   irsby = 1
   imrgw = 2

   delete(lat)
   delete(wgt)
 end do  ; ie

 if ( (nt - ntr) .ne. 0 ) then
   fu(0:3,:,0) = fu@_FillValue
;   fu(0:3,:,0) = fu(0:3,:,1)
 end if

;--------------------------------------- process 0
 nv = 9

 fu = fu*30.
 ut = ut*30.
 advz = advz*30.
; ureq = ureq*30.

; u3 = (/ut, advz, ureq/)

; uc = 

; determine criterion vari.
; cr = sh
 cr = new((/2,nt/),float)

 zeropt = new(nt,integer)
 zeropt = -999

 tag = 0
 nlast = -999
 nlast2 = -999
 cnt = 0
 do n=1, nt-1
   if ( uc(n)*uc(n-1) .le. 0. ) then
     zeropt(cnt) = n
     if (n-nlast2 .lt. 6) then
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

 ;yh +1 for 30 hPa
 nzero = cnt + 1
 zeropt(nzero-1) = nt-1

 maxmin = new((/2,nzero-1/),integer)
 do cnt=0, nzero-2
   if ( sum(uc(zeropt(cnt):zeropt(cnt+1))) .gt. 0. ) then
     maxmin(0,cnt) = maxind(uc(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = 1
   else
     maxmin(0,cnt) = minind(uc(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = -1
   end if
 end do

 cr = 0.
 do cnt=0, nzero-3
   cr(0,maxmin(0,cnt):maxmin(0,cnt+1)) = tofloat(maxmin(1,cnt))*(-1.0)
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

 npt = toint(max(ncyc))
 maxfu = new((/2,nv,ne,npt/),float,-999.)
 maxfu = -999.
 cnt1 = 0
 cnt2 = 0
 do cnt=0, nzero-3
   if ( cr(0,maxmin(0,cnt)) .gt. 0. ) then
     maxfu(0,:,:,cnt1) = dim_max(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     maxfu(0,iadvz,:,cnt1) = (-1.)*dim_min(fu(iadvz,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     maxfu(0,irsby,:,cnt1) = (-1.)*dim_min(fu(irsby,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     maxfu(0,imrgw,:,cnt1) = (-1.)*dim_min(fu(imrgw,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,:,:,cnt1) = dim_max(abs(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1))))
     cnt1 = cnt1 + 1
   else
     maxfu(1,:,:,cnt2) = dim_min(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     maxfu(1,iadvz,:,cnt2) = (-1.)*dim_max(fu(iadvz,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(1,:,:,cnt2) = dim_max(abs(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1))))
     cnt2 = cnt2 + 1
   end if
 end do

;--------------------------------------- cn level
 lcol = (/2,3,13,4,1/)  ; default
; lcol = (/15,9,13,3,1,17/)  ; amwg
; lcol = (/11,7,10,4,1,17/)  ; cosam
; lcol = (/11,7,9,4,1,17/)  ; radar
 ltck = (/1,1,1,1,1/)*0.8
; ltck = 1.2 + 0.2*fspan(0,-6,7)
 imrk = (/1,1,1,1,1/)*2

;--------------------------------------- contour
 wi = 60
 he = 40
 xf = 15 + (/0,1/)*(wi+12)
 yf = 10 + (/0,0/)*(he+10)

 ystr = "(m s:S:-1:N: month:S:-1:N:)"
; txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz","REQ","Kelvin","Rossby","MRG","IG","CGW+BGW","D:B:u:H-10V15:_:N:"/)
 txtstr = (/"Kelvin","Rossby","MRG","IG","EPD_total","ADVz","F:B:U:N:_total",":F34:6:F21:U/:F34:6:F21:t",":F34:6:F21:U/:F34:6:F21:t - F:B:U:N:_total"/)

do is=0, 1

print(max(maxfu(is,:,:,:))+"  "+min(maxfu(is,:,:,:)))

do iv=0, nv-1

 do ie=0, ne-1
   if (ie+iv .eq. 0) then
     plot = xylineplot(wks,iv+(ie-ne/2)*0.15,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
   else
     xyaddline(plot,iv+(ie-ne/2)*0.15,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
   end if
 end do

end do

  lcol2 = new((/nv*ne,npt/),integer)
  do iv=0, nv-1
  do i=0, npt-1
    lcol2(iv*ne:(iv+1)*ne-1,i) = lcol(:)
  end do
  end do

  xyaddmarker(plot,imrk,"S",0.8,"T",5.0,"C",lcol2)

  do iv=0, nv-2
    xyaddline(plot,iv+0.5,(/-999,999/),"D",0,"T",0.5,"C",colind)
  end do
  xyaddline(plot,ispan(-1,nv,1),0,"D",0,"T",0.5,"C",1)

 set_size(plot,xf(is),yf(is),wi,he)
 if (is .eq. 0) then
   axis_range(plot,-0.5,nv-0.5,-3,25,False,False)
 else
   axis_range(plot,-0.5,nv-0.5,3,-25,False,False)
 end if
 axis_str(plot,"","",ystr)
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(0,nv-1,1)
  "tmXBLabels" : txtstr
  "tmXBLabelAngleF" : 30
  "tmXBMinorPerMajor" : 0
  "tmXTOn" : False
;  "tmXBMajorLengthF" : -0.002
 end setvalues

 draw(plot)

; dr_txt2(wks,txtstr(iv),xf+1,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; is

 frame(wks)

end

