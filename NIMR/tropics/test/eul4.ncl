load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 dzsml = 10.

 wks = cr_wks("./eul4",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 n0 = 23
 n1 = 0
 n2 = 240
 files0 = systemfunc("ls ../dat/arch/uj/pa/uanuja.pa_195801??00.nc ../dat/arch/uj/pa/uanuja.pa_1958020100.nc")
 filesc = systemfunc("ls ../dat/arch/uj/pc/uanuja.pc_195801??00.nc ../dat/arch/uj/pc/uanuja.pc_1958020100.nc")
; files0 = systemfunc("ls /data7/kyh/arch/up/pa/uanupa.pa_194401??00.nc /data7/kyh/arch/up/pa/uanupa.pa_1944020100.nc")
; filesc = systemfunc("ls /data7/kyh/arch/up/pc/uanupa.pc_194401??00.nc /data7/kyh/arch/up/pc/uanupa.pc_1944020100.nc")

 f00 = addfiles(files0,"r")
 f0c = addfiles(filesc,"r")

do n=n1, n2
print("n="+n)

print("read")
 u3 = f00[:]->u    (n+n0,:,:,:)  ; 23 - 263
 v3 = f00[:]->v    (n+n0,:,:,:)
 w3 = f00[:]->dz_dt(n+n0,:,:,:)
 r3 = f00[:]->rho  (n+n0,:,:,:)
 p3 = f0c[:]->p    (n+n0,:,:,:)
print(".")

 if (n .eq. 0) then
   lon = f00[0]->longitude
   lat = f00[0]->latitude
   z   = f00[0]->hybrid_ht
   zt  = f00[0]->hybrid_ht_1
   nx = dimsizes(lon)
   ny = dimsizes(lat)
   nz = dimsizes(z  )
   nd = (/nz,ny/)
   r = z + 6371229.
   rcos = conform_dims(nd,r,0)*conform_dims(nd,cos(lat/180.*3.141592),1)
   dx = (lon(2)-lon(0))/180.*3.141592*rcos
   dy = (lat(2)-lat(0))/180.*3.141592*conform_dims(nd,r,0)
   dz0 = new(nz,float)
   do k=1, nz-2
     dz0(k) = z(k+1) - z(k-1)
   end do
   dz0(0   ) = dz0(1   )
   dz0(nz-1) = dz0(nz-2)
   dz = conform_dims(nd,dz0,0)
   f_r = 2.*(2.*3.141592/86400.)*sin(lat/180.*3.141592)
   f_y = 2.*(2.*3.141592/86400.)*cos(lat/180.*3.141592)
   f_y(0   ) = 0.
   f_y(ny-1) = 0.
   tanphi = tan(lat/180.*3.141592)

   cor  = new(nd,float)
   cor = 0.
   yadv  = cor
   zadv  = cor
   yadvp = cor
   zadvp = cor
   yadvm = cor
   zadvm = cor
   sml01 = cor
   sml02 = cor
   sml03 = cor
   sml04 = cor
 end if

; u4 = runave_n(u3,3,0,0)
; u3 = u4*1.5 - u3*0.5

; v4 = runave_n(v3,3,0,0)
; v3 = v4*1.5 - v3*0.5

 cor1 = dim_avg(v3)*conform_dims(nd,f_r,1)

 u3m = dim_avg(u3)
 u3p = u3 - conform_dims(dimsizes(u3),u3m,(/0,1/))
 v3m = dim_avg(v3)
 v3p = v3 - conform_dims(dimsizes(v3),v3m,(/0,1/))
 w3m = dim_avg(w3)
 w3p = w3 - conform_dims(dimsizes(w3),w3m,(/0,1/))

 tmp = (/u3/)
 tmpm = (/u3m/)
 do j=1, ny-2
   tmp(:,j,:) = u3p(:,j+1,:) - u3p(:,j-1,:)
   tmpm(:,j) = u3m(:,j+1) - u3m(:,j-1)
 end do
 tmp(:,0   ,:) = 0.
 tmp(:,ny-1,:) = 0.
 tmpm(:,0   ) = 0.
 tmpm(:,ny-1) = 0.
 yadv1p = (-1.)*dim_avg( tmp*v3p )/dy
 yadv1m = (-1.)*tmpm*v3m/dy
 yadv1 = yadv1p + yadv1m

; w4 = runave(w3,3,-1)
; w3 = w4*1.5 - w3*0.5
; w4 = runave_n(w3,3,0,0)
; w3 = w4*1.5 - w3*0.5
 tmp = 0.
 tmpm = 0.
 do k=1, nz-2
   tmp(k,:,:) = u3p(k+1,:,:) - u3p(k-1,:,:)
   tmpm(k,:) = u3m(k+1,:) - u3m(k-1,:)
 end do
 tmp(0   ,:,:) = 0.
 tmp(nz-1,:,:) = 0.
 tmpm(0   ,:) = 0.
 tmpm(nz-1,:) = 0.
 zadv1p = (-1.)*dim_avg( tmp*w3p )/dz
 zadv1m = (-1.)*tmpm*w3m/dz
 zadv1 = zadv1p + zadv1m
; zzz = new(nz*2,float)
; zzz(:nz-1) = z-0.5*dzsml
; zzz(nz:  ) = z+0.5*dzsml
; u3tmp = new((/ny,nx,nz/),float)
; do k=0, nz-1
;   u3tmp(:,:,k) = u3(k,:,:)
; end do
; tmp1 = ftcurv(z,u3tmp,zzz)
; tmp = (/u3/)
; do k=0, nz-1
;   tmp(k,:,:) = (tmp1(:,:,k+nz) - tmp1(:,:,k))/dzsml
; end do
; zadv1 = (-1.)*dim_avg( tmp*w3 )

 sml1 = (-1.)*dim_avg(w3)*conform_dims(nd,f_y,1)
 sml2 = dim_avg(u3*v3)*conform_dims(nd,tanphi,1)/conform_dims(nd,r,0)
 sml3 = (-1.)*dim_avg(u3*w3)/conform_dims(nd,r,0)
 sml4 = sml1*0.
 tmp = 0.
 do i=1, nx-2
   tmp(:,:,i) = p3(:,:,i+1) - p3(:,:,i-1)
 end do
 tmp(:,:,0   ) = p3(:,:,1) - p3(:,:,nx-1)
 tmp(:,:,nx-1) = p3(:,:,0) - p3(:,:,nx-2)
 sml4 = (-1.)*dim_avg( tmp/r3 )/dx

 wgt = 86400./tofloat(n2-n1)
 if ( n .eq. n1 .or. n .eq. n2 ) then
   wgt = wgt*0.5
   if (n .eq. n1) then
     du = dim_avg(u3)
   end if
   if (n .eq. n2) then
     du = (dim_avg(u3) - du)/tofloat((n2-n1)/8)
   end if
 end if

 cor   = cor   + cor1*wgt
 yadv  = yadv  + yadv1*wgt
 zadv  = zadv  + zadv1*wgt
 yadvp = yadvp + yadv1p*wgt
 zadvp = zadvp + zadv1p*wgt
 yadvm = yadvm + yadv1m*wgt
 zadvm = zadvm + zadv1m*wgt
 sml01 = sml01 + sml1*wgt
 sml02 = sml02 + sml2*wgt
 sml03 = sml03 + sml3*wgt
 sml04 = sml04 + sml4*wgt

end do

 adv = yadv + zadv
 sml = sml01 + sml02 + sml03 + sml04
 tot0 = cor + adv + sml

 f0 = addfile("../dat/arch/uj/pe/uanuja.pe_195801.nc","r")
; f1 = addfile("../dat/arch/uj/pa/uanuja.pa_1958010100.nc","r")
; f2 = addfile("../dat/arch/uj/pa/uanuja.pa_1958020100.nc","r")
; f3 = addfile("../dat/tem/1958/uanuj.tem_yzt.1958.01.nc","r")
 f4 = addfile("../dat/arch/uj/pj/uanuja.pj_195801.nc","r")
um = dim_avg( f4->u(0,:,:,:) )

; f0 = addfile("/data7/kyh/arch/up/pe/uanupa.pe_194401.nc","r")
; f1 = addfile("/data7/kyh/arch/up/pa/uanupa.pa_1944010100.nc","r")
; f2 = addfile("/data7/kyh/arch/up/pa/uanupa.pa_1944020100.nc","r")
; f3 = addfile("../dat/tem_up2/1944/uanup.tem_yzt.1944.01.nc","r")
; f4 = addfile("/data7/kyh/arch/up/pj/uanupa.pj_194401.nc","r")

 y = lat
 z = z/1.e3

 var = new((/6,nz,ny/),float)

 var(0,:,:) = dim_avg(f0->du_gw    (0,:,:,:))*72.
 var(0,:,:) = var(0,:,:) + dim_avg(f0->dudt_gwconv(0,:,:,:))*86400.
 var(1,:,:) = dim_avg(f0->du_adv   (0,:,:,:))*72.
 var(2,:,:) = dim_avg(f0->du_bl    (0,:,:,:))*72.
 var(3,:,:) = dim_avg(f0->du_solver(0,:,:,:))*72.
 var(4,:,:) = dim_avg(f0->dudt_conv_p(0,:,:,:))*86400.

 var(5,:,:) = dim_sum_n(var(0:4,:,:),0)

 tot = tot0 + var(0,:,:) + var(2,:,:) + var(4,:,:)

 dif = tot-du
 dif(0   ,:) = 0.
 dif(nz-1,:) = 0.
;--------------------------------------- process 0
print(sum(abs(var(5,:,:)-du)))

 fvar = (/ du, tot, dif, \
           var(1,:,:)-cor, tot0-cor, var(0,:,:), \
           cor, adv, yadvp, \
           yadvm,zadvm,zadvp/)

;--------------------------------------- cn level
 cnlev = fspan(-0.5,0.5,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlev2= (/-16,-8,-4,-2,-1,1,2,4,8,16/)

;--------------------------------------- contour
 nl = 12

tit = (/"tend","F_total","Diff","MADV-COR","EUL-COR","GWD","COR","ADV","ADVY_P","ADVY_M","ADVZ_M","ADVZ_P"/)

 wi = 60
 he = 40
 xf = 34 + (/0,1,2,0,1,2,0,1,2,0,1,2/)*70
 yf = 10 + (/0,0,0,1,1,1,2,2,2,3,3,3/)*60

do il=0, nl-1

 contour1 = cnshadeplot(wks,y,z,fvar(il,:,:),cnlev,cnfil," ")
 cnaddline(contour1,wks,y,z,fvar(il,:,:),cnlev2,0,"D",0,"T",1,"C",1)

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-90,90,20,60,False,False)
 axis_str(contour1,tit(il),"","")

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end



