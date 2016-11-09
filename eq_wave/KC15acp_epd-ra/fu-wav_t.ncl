load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 expname = (/"era-int","merra","cfsr","jra55"/);,"era-int_f"/)
 datname0= (/"ERA-I"  ,"MERRA","CFSR","JRA55"/);,"ERA-I-nr"/)
 datname = (/"ERA-I"  ,"MERRA","CFSR","JRA-55"/);,"ERA-I_ml"/)

 y0 = 1979
 feb_1st = True
 y1 = 1979
 y2 = 2010
 lat0 = 5 ;10
 lat0u = 5
 lata = lat0 + 2.5
 p = 30 ;30
; p_em = 19  ; (50,30,20,10) : (22,19,17,14)

; ymax = (/7,0.5,6,2,10,30,30,20/)  ; for lat0 = 10
 ymax = (/12,2,6,2,12,30,30,20/)   ; for lat0 = 5
 ymin = (/-2,-2.4,-4,-4,-6,30,30,20/)   ; for lat0 = 5

 wks = cr_wks("fig/fu-wav_t-p"+p+"ya"+lat0,600,0)
; gsn_define_colormap(wks,"amwg")
; gsn_define_colormap(wks,"cosam")
; gsn_define_colormap(wks,"radar")
 icol_e = NhlNewColor(wks,32/256.,120/256.,180/256.)
 icol_m = NhlNewColor(wks,226/256.,31/256.,38/256.)
 icol_c = NhlNewColor(wks,52/256.,160/256.,72/256.)
 icol_j = NhlNewColor(wks,245/256.,126/256.,32/256.)
 colind = NhlNewColor(wks,0.80,0.80,0.80)
 colind2 = NhlNewColor(wks,0.65,0.65,0.65)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nw = 8
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12
 ntr = nt
 if ( feb_1st .and. ( n0 .eq. 0 ) ) then
   ntr = nt - 1
 end if

 u_expname = (/"era-int","merra","cfsr","jra55"/);,"era-int_f"/)
 u_datname0= (/"ERA-I"  ,"MERRA","CFSR","JRA55"/);,"ERA-I-nr"/)
 u_vname = (/"U_GDS0_ISBL","u","UGRD","U_GRD_GDS0_ISBL"/);,"u"/)
 u_fnam0 = "$DATD/clim/"+u_datname0+"/"+u_expname+".u_ypt.nc"
 uc = new(nt,float)
 uc = 0.
 do ie=0, dimsizes(u_expname)-1
   f0 = addfile(u_fnam0(ie),"r")
   s_lat = f0->$u_vname(ie)$!2
   lat = tofloat( f0->$s_lat$({-lata:lata}) )
   u0 = new(nt,float)
;   if (u_expname(ie) .eq. "era-int_f") then
;     do n=0, nt-1
;       u0(n) = ftcurvi(-lat0u,lat0u,lat, \
;               tofloat( f0->$u_vname(ie)$(n0+n,{p_em},{-lata:lata}) ))
;     end do
;   else
     do n=0, nt-1
       u0(n) = ftcurvi(-lat0u,lat0u,lat, \
               tofloat( f0->$u_vname(ie)$(n0+n,{p},{-lata:lata}) ))
     end do
;   end if
   u0 = u0 / (lat0u+lat0u)
   uc = uc + u0/tofloat(dimsizes(u_expname))
   delete(lat)
 end do

 fnam0 = "$DATD/clim/"+datname0+"/"+expname+".u_ypt.nc"
 fnam1 = "$DATD/tem/"+datname0+"/epf/"+expname+".epfsa_yz_recon0.nc"
 fnam2 = "$DATD/tem/"+datname0+"/epf/"+expname+".epf_yz_recon0.nc"
 fnam3 = "$DATD/tem/"+datname0+"/"+expname+".tem_ypt.nc"

 do ie=0, ne-1
   f0 = addfile(fnam0(ie),"r")
   f1 = addfile(fnam1(ie),"r")
   f2 = addfile(fnam2(ie),"r")
   f3 = addfile(fnam3(ie),"r")

   lat = f1->lat({-lata:lata})
   z = 18.6 + 6.6*log( 70./ f1->p({p}) )
   zp = 18.6 + 6.6*log( 70./ p )
if ( abs(z-zp) .lt. 0.5 ) then
   print(z+"")
else   ; ERA-I-nr can be
   z0 = z
   k0 = ind( f1->p .eq. f1->p({p}) )
   if ( z0 .gt. zp ) then
     k1 = -1
   else
     k1 = 1
   end if
   if ( f1->p(1) - f1->p(0) .gt. 0. ) then
     k1 = k1*(-1)
   end if
   k1 = k1 + k0
   z1 = 18.6 + 6.6*log( 70./ f1->p(k1) )
   print(z0+"  "+zp+"  "+z1)
   c0 = abs((z1-zp)/(z1-z0))
   c1 = 1. - c0
end if

   rw_i = new((/ntr,nw/),float)
if ( abs(z-zp) .lt. 0.5 ) then
   do n=0, ntr-1
     do iw=0, nw-1
       rw_i(n,iw) = ftcurvi(-lat0,lat0,lat, f1->epd(n0+n,iw,{p},{-lata:lata}))
     end do
     rw_i(n,nw-2) = rw_i(n,nw-2) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0+n,8,{p},{-lata:lata}))
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0+n,9,{p},{-lata:lata}))

; migrating tide
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0+n,2,{p},{-lata:lata}) + \
                                            f2->epd(n0+n,3,{p},{-lata:lata}))

     rw_i(n,3) = rw_i(n,3) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0+n,13,{p},{-lata:lata}))
     rw_i(n,6) = rw_i(n,6) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0+n,16,{p},{-lata:lata}))
     rw_i(n,7) = rw_i(n,7) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0+n,17,{p},{-lata:lata}))
   end do
else   ; ERA-I-nr can be
   do n=0, ntr-1
     do iw=0, nw-1
       rw_i(n,iw) = ftcurvi(-lat0,lat0,lat, \
                            c0*f1->epd(n0+n,iw,k0,{-lata:lata}) + \
                            c1*f1->epd(n0+n,iw,k1,{-lata:lata}))
     end do
     rw_i(n,nw-2) = rw_i(n,nw-2) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0+n,8,k0,{-lata:lata}) + \
                            c1*f2->epd(n0+n,8,k1,{-lata:lata}))
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0+n,9,k0,{-lata:lata}) + \
                            c1*f2->epd(n0+n,9,k1,{-lata:lata}))

; migrating tide
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0+n,2,k0,{-lata:lata}) + \
                            c1*f2->epd(n0+n,2,k1,{-lata:lata}) + \
                            c0*f2->epd(n0+n,3,k0,{-lata:lata}) + \
                            c1*f2->epd(n0+n,3,k1,{-lata:lata}))

     rw_i(n,3) = rw_i(n,3) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0+n,13,k0,{-lata:lata}) + \
                         c1*f1->epd(n0+n,13,k1,{-lata:lata}))
     rw_i(n,6) = rw_i(n,6) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0+n,16,k0,{-lata:lata}) + \
                         c1*f1->epd(n0+n,16,k1,{-lata:lata}))
     rw_i(n,7) = rw_i(n,7) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0+n,17,k0,{-lata:lata}) + \
                         c1*f1->epd(n0+n,17,k1,{-lata:lata}))
   end do
end if

   rw_i = rw_i / (lat0+lat0)
   rw0 = new((/nw,nt/),float,1.e20)
   rw0(:,nt-ntr:nt-1) = transpose(rw_i)

   ut   = new(nt,float)
   advz = new(nt,float)
   uf   = new(nt,float)
   ue   = new(nt,float)
   advy = new(nt,float)
if ( abs(z-zp) .lt. 0.5 ) then
   do n=0, nt-1
     ut  (n) = ftcurvi(-lat0,lat0,lat, f3->u_tend (n0+n,{p},{-lata:lata}))
     advz(n) = ftcurvi(-lat0,lat0,lat, f3->uadv_z (n0+n,{p},{-lata:lata}))
     uf  (n) = ftcurvi(-lat0,lat0,lat, f3->u_force(n0+n,{p},{-lata:lata}))
     ue  (n) = ftcurvi(-lat0,lat0,lat, f3->epd    (n0+n,{p},{-lata:lata}))
     advy(n) = ftcurvi(-lat0,lat0,lat, f3->uadv_y (n0+n,{p},{-lata:lata}) + \
                                       f3->cor    (n0+n,{p},{-lata:lata}))
   end do
else   ; ERA-I-nr can be
   do n=0, nt-1
     ut  (n) = ftcurvi(-lat0,lat0,lat, \
                       c0*f3->u_tend (n0+n,k0,{-lata:lata}) + \
                       c1*f3->u_tend (n0+n,k1,{-lata:lata}))
     advz(n) = ftcurvi(-lat0,lat0,lat, \
                       c0*f3->uadv_z (n0+n,k0,{-lata:lata}) + \
                       c1*f3->uadv_z (n0+n,k1,{-lata:lata}))
     uf  (n) = ftcurvi(-lat0,lat0,lat, \
                       c0*f3->u_force(n0+n,k0,{-lata:lata}) + \
                       c1*f3->u_force(n0+n,k1,{-lata:lata}))
     ue  (n) = ftcurvi(-lat0,lat0,lat, \
                       c0*f3->epd    (n0+n,k0,{-lata:lata}) + \
                       c1*f3->epd    (n0+n,k1,{-lata:lata}))
     advy(n) = ftcurvi(-lat0,lat0,lat, \
                       c0*f3->uadv_y (n0+n,k0,{-lata:lata}) + \
                       c1*f3->uadv_y (n0+n,k1,{-lata:lata}) + \
                       c0*f3->cor    (n0+n,k0,{-lata:lata}) + \
                       c1*f3->cor    (n0+n,k1,{-lata:lata}))
   end do
end if
   ut   = ut   / (lat0+lat0)
   advz = advz / (lat0+lat0)
   uf   = uf   / (lat0+lat0)
   ue   = ue   / (lat0+lat0)
   advy = advy / (lat0+lat0)

;   ureq = (/ ut - advz - advy /)

   if (ie .eq. 0) then
     fu = new((/8,ne,nt/),float)
     t = y1 + (fspan(1,nt,nt)-0.5)/12.
;     tt = y1 + (fspan(1,nt-1,nt-1))/12.
   end if
   fu(:,ie,:) = (/ rw0(2,:), rw0(4,:)+rw0(5,:), \
                   rw0(6,:)+rw0(7,:), rw0(3,:), ue, advz, uf, ut - uf /)

   delete(lat)
 end do  ; ie

 if ( (nt - ntr) .ne. 0 ) then
   fu(0:3,:,0) = fu@_FillValue
;   fu(0:3,:,0) = fu(0:3,:,1)
 end if

;--------------------------------------- process 0
 nv = 8

 fu = fu*30.
 ut = ut*30.
 advz = advz*30.
; ureq = ureq*30.

 polylx = new(nt*2,float)
 polyly = new((/nv,ne,nt*2/),float)
 do i=0, nt-1
   polylx(i*2  ) = t(i) - 0.5*(t(1) - t(0))
   polylx(i*2+1) = t(i) + 0.5*(t(1) - t(0))
 end do
 do i=0, nt-1
   polyly(:,:,i*2  ) = fu(:,:,i)
   polyly(:,:,i*2+1) = fu(:,:,i)
 end do

 polygx = new(nt*4+1,float)
 polygy = new((/nv,nt*4+1/),float)
 polygx(:nt*2-1) = polylx(:)
 polygx(nt*2:nt*4-1) = polylx(::-1)
 polygx(nt*4) = polygx(0)
 polygy(:,:nt*2-1) = dim_max_n(polyly(:,:,:),1)
 polygy(:,nt*2:nt*4-1) = dim_min_n(polyly(:,:,::-1),1)
 polygy(:,nt*4) = polygy(:,0)

 polygx2 = polygx
 polygy2 = conform_dims((/nv,nt*4+1/),ymax(:)*1.,0)
 polygy2(:,nt*2:nt*4-1) = polygy2(:,nt*2:nt*4-1) + \
      polygy(:,nt*2:nt*4-1) - polygy(:,0:nt*2-1:-1)

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

 iv = 0  ; 0-4
 cnt0 = 0  ; E-W
; cnt0 = 1  ; W-E
 do cnt=cnt0, nzero-3, 2
   print( maxmin(1,cnt)+"  "+t(maxmin(0,cnt))+"  "+t(maxmin(0,cnt+1))+" :  "+ \
          min( dim_max(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) )+"  "+ \ 
          max( dim_max(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) ) )
 end do
;  p-lev, +
;          min( dim_max(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) )+"  "+ \ 
;          max( dim_max(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) ) )
;  p-lev, -
;          min( dim_min(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) )+"  "+ \ 
;          max( dim_min(fu(iv,:3,maxmin(0,cnt):maxmin(0,cnt+1))) ) )
;X  m-lev, +
;X          max(fu(iv,4,maxmin(0,cnt):maxmin(0,cnt+1))) )
;X  m-lev, -
;X          min(fu(iv,4,maxmin(0,cnt):maxmin(0,cnt+1))) )

;--------------------------------------- cn level
 lcol = (/4,2,3,13,1,7/)  ; default
; lcol = (/icol_e,icol_m,icol_c,icol_j,1,7/)
; lcol = (/13,15,9,3,1,17/)  ; amwg
; lcol = (/10,11,7,4,1,17/)  ; cosam
; lcol = (/9,11,7,4,1,17/)  ; radar
 ltck = (/1,1,1,1,1/)*0.75
; ltck = 1.2 + 0.2*fspan(0,-6,7)
 imrk = (/1,1,1,1,1/)

;--------------------------------------- contour
 wi = 140
 he = 15
 xf = 15
 yf = 10 + ispan(0,nv-1,1)*(he+5)

 xstr = new(nv,string)
 xstr(:) = ""
; xstr(nv-1) = "year"
 xstr(4) = "year"
 ystr = "(m s:S:-1:N: month:S:-1:N:)"
; txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz","REQ","Kelvin","Rossby","MRG","IG","CGW+BGW","D:B:u:H-10V15:_:N:"/)
 txtstr = (/"Kelvin","MRG","IG","Rossby","Net_resolved","ADVz","F:B:U:N:_total","X"/)

do iv=0, 4 ;nv-1

 plot = xylineplot(wks,t,fu(iv,:,:),"D",0,"T",0.5,"C",-1)
; plot = xylineplot(wks,t,fu(iv,:,:),"D",0,"T",0.4,"C",lcol)
; xyaddmarker(plot,imrk,"S",0.8,"T",10.0,"C",lcol)
 respg = True
 respg@gsFillColor = colind
; dumpg = gsn_add_polygon(wks,plot,polygx,polygy(iv,:),respg)
 respg2 = True
 respg2@gsFillColor = colind2
 dumpg2 = gsn_add_polygon(wks,plot,polygx2,polygy2(iv,:),respg2)
 dumpl = new(ne,graphic)
 ie4 = (/2,1,3,0/)
 do ie0=0, ne-1
   respl = True
   respl@gsLineColor      = lcol(ie4(ie0))
   respl@gsLineThicknessF = ltck(ie4(ie0))
;   respl@tfPolyDrawOrder = "Draw"
   dumpl(ie0) = gsn_add_polyline(wks,plot,polylx,polyly(iv,ie4(ie0),:),respl)
 end do

 xyaddline(plot,t,0.,"D",0,"T",0.6,"C",1)
 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,y1,y2+1,ymin(iv),ymax(iv),False,False)
 axis_str(plot,"",xstr(iv),ystr)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 do cnt=0, nzero-3
   if (maxmin(1,cnt) .eq. 1) then
     vldash = 0
   else
     vldash = 2
   end if
   xyaddline(plot,t(maxmin(0,cnt)),(/-999,999/),"D",vldash,"T",0.5,"C",1)
 end do

 draw(plot)

 dr_txt2(wks,txtstr(iv),xf+1,yf(iv)+he-1.5,22,6,0.,"CenterLeft",1)
; dr_txt2(wks,txtstr(iv),xf+1,yf(iv)+1.8,22,6,0.,"CenterLeft",1)

; delete(dumpg)
 delete(dumpg2)
 delete(dumpl)

end do

 axisfontratio = 1.2

 do ie=0, ne-1
 lg = lgplot(wks,datname(ie),"D",0,"T",1.5,"C",lcol(ie),2)
 set_size_lg(lg,xf+ie*35,yf(5)+2,wi*0.15,3)
 draw(lg)
 end do

 frame(wks)

end

