load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute0.ncl"

;axisfontratio = 0.85

begin

;--------------------------------------- parameters
 expname = (/"era-int","merra","merra2","cfsr","jra55","era-int_f","merra2_f" ,"jra55_f"/)
 datname0= (/"ERA-I"  ,"MERRA","MERRA2","CFSR","JRA55","ERA-I-nr" ,"MERRA2-nr","JRA55-nr"/)
 datname = (/"ERA-I"  ,"MERRA","MERRA2","CFSR","JRA55","ERA-I_ml" ,"MERRA2_ml","JRA55_ml"/)
 uname = (/"U_GDS0_ISBL","u","U","UGRD","U_GRD_GDS0_ISBL","u","U","U_GRD_GDS4_HYBL"/)
 l_ml = (/False,False,False,False,False,True,True,True/)
 y0 = (/1979,1979,1980,1979,1979,1979,1980,1979/)

 feb_1st = True
 y1 = 1981
 y2 = 2010
 lat0 = 5.
 lata = lat0 + 2.5
 p = 20 ;30

 wks = cr_wks("fig/fu-all2-p"+p,600,0)
; gsn_define_colormap(wks,"amwg")
; gsn_define_colormap(wks,"cosam")
; gsn_define_colormap(wks,"radar")
 load "./ftn_get_color_S-RIP.ncl"
 icol = get_col_srip(wks)
 colind = NhlNewColor(wks,0.80,0.80,0.80)
 colind2 = NhlNewColor(wks,0.65,0.65,0.65)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nw = 8
 nt = (y2-y1+1)*12
 n0s = (y1-y0)*12

 fnam0 = "$DATD/clim/"+datname0+"/"+expname+".u_ypt.nc"
 fnam1 = "$DATD/tem/"+datname0+"/epf/"+expname+".epfsa_yz_recon0.nc"
 fnam2 = "$DATD/tem/"+datname0+"/epf/"+expname+".epf_yz_recon0.nc"
 fnam3 = "$DATD/tem/"+datname0+"/"+expname+".tem_ypt.nc"

 do ie=0, ne-1
   n0 = n0s(ie)

   n0r = n0
   if ( feb_1st ) then
     n0r = n0 - 1
   end if

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

   rw_i = new((/nt,nw/),float)
if ( abs(z-zp) .lt. 0.5 ) then
   do n=0, nt-1
     if (n0r .lt. 0) then
       rw_i(n,:) = rw_i@_FillValue
       continue
     end if
     do iw=0, nw-1
       rw_i(n,iw) = ftcurvi(-lat0,lat0,lat, f1->epd(n0r+n,iw,{p},{-lata:lata}))
     end do
     rw_i(n,nw-2) = rw_i(n,nw-2) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0r+n,8,{p},{-lata:lata}))
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0r+n,9,{p},{-lata:lata}))

; migrating tide
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, f2->epd(n0r+n,2,{p},{-lata:lata}) + \
                                            f2->epd(n0r+n,3,{p},{-lata:lata}))

     rw_i(n,3) = rw_i(n,3) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0r+n,13,{p},{-lata:lata}))
     rw_i(n,6) = rw_i(n,6) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0r+n,16,{p},{-lata:lata}))
     rw_i(n,7) = rw_i(n,7) + \
                 ftcurvi(-lat0,lat0,lat, f1->epd(n0r+n,17,{p},{-lata:lata}))
   end do
else   ; ERA-I-nr can be
   do n=0, nt-1
     if (n0r .lt. 0) then
       continue
     end if
     do iw=0, nw-1
       rw_i(n,iw) = ftcurvi(-lat0,lat0,lat, \
                            c0*f1->epd(n0r+n,iw,k0,{-lata:lata}) + \
                            c1*f1->epd(n0r+n,iw,k1,{-lata:lata}))
     end do
     rw_i(n,nw-2) = rw_i(n,nw-2) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0r+n,8,k0,{-lata:lata}) + \
                            c1*f2->epd(n0r+n,8,k1,{-lata:lata}))
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0r+n,9,k0,{-lata:lata}) + \
                            c1*f2->epd(n0r+n,9,k1,{-lata:lata}))

; migrating tide
     rw_i(n,nw-1) = rw_i(n,nw-1) + \
                    ftcurvi(-lat0,lat0,lat, \
                            c0*f2->epd(n0r+n,2,k0,{-lata:lata}) + \
                            c1*f2->epd(n0r+n,2,k1,{-lata:lata}) + \
                            c0*f2->epd(n0r+n,3,k0,{-lata:lata}) + \
                            c1*f2->epd(n0r+n,3,k1,{-lata:lata}))

     rw_i(n,3) = rw_i(n,3) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0r+n,13,k0,{-lata:lata}) + \
                         c1*f1->epd(n0r+n,13,k1,{-lata:lata}))
     rw_i(n,6) = rw_i(n,6) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0r+n,16,k0,{-lata:lata}) + \
                         c1*f1->epd(n0r+n,16,k1,{-lata:lata}))
     rw_i(n,7) = rw_i(n,7) + \
                 ftcurvi(-lat0,lat0,lat, \
                         c0*f1->epd(n0r+n,17,k0,{-lata:lata}) + \
                         c1*f1->epd(n0r+n,17,k1,{-lata:lata}))
   end do
end if

   rw_i = rw_i / (lat0+lat0)
   rw0 = transpose(rw_i)

   if ( .not. l_ml(ie) ) then
     u0 = new(nt,float)
     do n=0, nt-1
       u0(n) = ftcurvi(-lat0,lat0,lat, tofloat( f0->$uname(ie)$(n0+n,{p},{-lata:lata}) ))
     end do
     u0 = u0 / (lat0+lat0)
     if ( .not. isdefined("uc") ) then
       uc = new(nt,float)
       uc = 0.
     end if
     uc = uc + u0/tofloat(ne-1)
   end if

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
     fu = new((/7,ne,nt/),float)
     fua = new((/3,ne,nt/),float)
     t = y1 + (fspan(1,nt,nt)-0.5)/12.
;     tt = y1 + (fspan(1,nt-1,nt-1))/12.
   end if
   fu(:,ie,:) = (/ rw0(2,:), rw0(3,:), rw0(4,:)+rw0(5,:), \
                   rw0(6,:)+rw0(7,:), ue, advz+advy, ut - uf /)
   fua(:,ie,:) = (/advz+advy, ut, uf/)

   delete(lat)
 end do  ; ie

; uf modification
 do ie=0, ne-1
   if (datname0(ie) .eq. "ERA-I-nr") then
     i_eraf = ie
   end if
 end do
 do ie=0, ne-1
   if (datname0(ie) .ne. "ERA-I-nr") then
;     fu(6,ie,:) = fu(6,ie,:) + fua(0,i_eraf,:) - fua(0,ie,:)
     fu(5,ie,:) = fu(5,ie,:) - fua(0,i_eraf,:) + fua(0,ie,:)
     fu(5,ie,:) = fu(5,ie,:) + fua(1,i_eraf,:) - fua(1,ie,:)
   end if
 end do

;--------------------------------------- process 0
 nv = 7

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
     maxfu(0,:,:,cnt1) = dim_avg(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,iadvz,:,cnt1) = (-1.)*dim_avg(fu(iadvz,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,irsby,:,cnt1) = (-1.)*dim_avg(fu(irsby,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,imrgw,:,cnt1) = (-1.)*dim_avg(fu(imrgw,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,:,:,cnt1) = dim_max(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,1,:,cnt1) = dim_min(fu(1,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,2,:,cnt1) = dim_min(fu(2,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(0,5,:,cnt1) = dim_min(fu(5,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     cnt1 = cnt1 + 1
   else
     maxfu(1,:,:,cnt2) = dim_avg(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(1,iadvz,:,cnt2) = (-1.)*dim_avg(fu(iadvz,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(1,:,:,cnt2) = dim_min(fu(:,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(1,0,:,cnt1) = dim_max(fu(0,:,maxmin(0,cnt):maxmin(0,cnt+1)))
;     maxfu(1,5,:,cnt1) = dim_max(fu(5,:,maxmin(0,cnt):maxmin(0,cnt+1)))
     cnt2 = cnt2 + 1
   end if
 end do

;--------------------------------------- cn level
; lcol = (/4,2,3,13,1/)  ; default
 lcol = (/icol@ei,icol@m,icol@m2,icol@c,icol@j5,icol@ei,icol@m2,icol@j5/)
; lcol = (/15,9,13,3,1,17/)  ; amwg
; lcol = (/11,7,10,4,1,17/)  ; cosam
; lcol = (/11,7,9,4,1,17/)  ; radar
 ltck = (/1,1,1,1,1,1,1,1/)*0.8
; ltck = 1.2 + 0.2*fspan(0,-6,7)
 imrk = (/2,2,2,2,2,2,2,2/)

;--------------------------------------- contour
 wi = 70
 he = (/32,32/)
 xf = 15 + (/0,0/)*(wi+12)
 yf = 10 + (/0,1/)*(he(0)+10)

 ystr = "force (m s:S:-1:N: month:S:-1:N:)"
; txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz","REQ","Kelvin","Rossby","MRG","IG","CGW+BGW","D:B:u:H-10V15:_:N:"/)
 txtstr = (/"Kelvin","Rossby","MRG","IG","Net_resolved","ADV","X"/)

do is=0, 1

print(max(maxfu(is,:,:,:))+"  "+min(maxfu(is,:,:,:)))

do iv=0, nv-1

 do ie=0, ne-1
   if (ie+iv .eq. 0) then
;     plot = xylineplot(wks,iv+(ie-ne/2)*0.10,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
     plot = xylineplot(wks,iv+(ie-ne/2+0.5)*0.10,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
   else
;     xyaddline(plot,iv+(ie-ne/2)*0.10,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
     xyaddline(plot,iv+(ie-ne/2+0.5)*0.10,maxfu(is,iv,ie,:),"D",0,"T",1,"C",-1)
   end if
 end do

end do

  lcol2 = new((/nv*ne,npt/),integer)
  do iv=0, nv-1
  do i=0, npt-1
    lcol2(iv*ne:(iv+1)*ne-1,i) = lcol(:)
  end do
  end do

  xyaddmarker(plot,imrk,"S",0.8,"T",8.0,"C",lcol2)

  do iv=0, nv-2
    xyaddline(plot,iv+0.5,(/-999,999/),"D",0,"T",0.5,"C",1)
  end do
  xyaddline(plot,ispan(-1,nv,1),0,"D",0,"T",0.5,"C",1)

 pgres = True
 pgres@gsFillColor = colind
 dumm = new(nv,graphic)
 do iv=0, nv-1
;   tmp = iv+(-ne/2)*0.10
   tmp = iv-0.5
   tmp2 = iv+(dimsizes(ind(.not. l_ml))-ne/2)*0.10
   dumm(iv) = gsn_add_polygon(wks,plot,(/tmp,tmp,tmp2,tmp2/),(/-100,100,100,-100/),pgres)
 end do

 set_size(plot,xf(is),yf(is),wi,he(is))
 if (is .eq. 0) then
   if (p .ge. 30) then
     axis_range(plot,-0.5,nv-0.5,-8,10,False,False) ; 30
   else if (p .eq. 20) then
     axis_range(plot,-0.5,nv-0.5,-11,10,False,False) ; 20
   else
     axis_range(plot,-0.5,nv-0.5,-20,15,False,False) ; 10
   end if
   end if
 else
   if (p .ge. 30) then
     axis_range(plot,-0.5,nv-0.5,8,-10,False,False) ; 30
   else if (p .eq. 20) then
     axis_range(plot,-0.5,nv-0.5,11,-10,False,False) ; 20
   else
     axis_range(plot,-0.5,nv-0.5,15,-20,False,False) ; 10
   end if
   end if
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

