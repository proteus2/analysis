load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.4

cx = 10 ;500.e3 / 

begin

 wks = cr_wks("r1_z",600,0)
 gsn_define_colormap(wks,"gsdtol")
 icol_center = 11
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
 pi = 3.141592

 f0 = addfile("../res/uvn_background_p1_region.nc","r")

 z  = f0->height_m({20000:}) / 1.e3
 ub = f0->ub(:,{20000:},:,:)
 nb = f0->nb(:,{20000:},:,:)

 nz = dimsizes(z)

 u0 = new(nz,float,1.e32)
 n0 = new(nz,float,1.e32)
 do k=0, nz-1
   u0(k) = avg(ub(:,k,:,:))
   n0(k) = avg(nb(:,k,:,:))
 end do

 m0 = abs(n0/(cx-u0))
 m0i= 1./m0
 l0 = 2.*pi/m0

 ratio = 100
 nz2 = (nz-1)*ratio+1
 d2z = 1./ratio
 z2 = new(nz2,float,1.e32)
 do k=0, nz2-1
   z2(k) = z(0)+0.5*d2z*k
 end do

;;;
 m0 = 0.5*(9.+6.3)*1.e-4 - 0.5*(9.-6.3)*1.e-4*tanh(z-29)
 u0 = -20 ;0.5*(-17.-21.) - 0.5*(-17.+21.)*tanh(z-29)
 n0 = abs(m0*(cx-u0))
;;;

 s0 = ftcurv(z,n0*(cx-u0),z2)
 m02= ftcurv(z,m0,z2)

 r0r = new(nz2,float,1.e32)
 r0r(nz2-2:nz2-1) = 0.
 r0i = r0r
 do k=0, nz2-3
   kk = nz2-3-k
   r0r(kk) = r0r(kk+2)+d2z*( -0.5*log(s0(kk+2)/s0(kk))/d2z*(1.-r0r(kk+1)^2+r0i(kk+1)^2) - 2.*m02(kk+1)*r0i(kk+1) )
   r0i(kk) = r0i(kk+2)+d2z*( -0.5*log(s0(kk+2)/s0(kk))/d2z*(-2.*r0r(kk+1)*r0i(kk+1)) + 2.*m02(kk+1)*r0r(kk+1) )
 end do
 r0 = sqrt(r0r^2+r0i^2)

;--------------------------------------- line attribute
 cnlev1 = (/0.05,0.1,0.15,0.2,0.3,0.4,0.6/)
 cnfil1 = (/0,28,24,20,16,12,8,4/)
;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 12 + (/0,1,2/)*53
 yf = 10

 mstr = (/"","MFx_T_up","T variance"/)


 plot = xylineplot(wks,m02^2,z2,"D",0,"T",1.2,"C",1)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,1.e-7,2.e-6,-999,45,True,False)
 axis_str(plot,"","m:S:2:N: [rad:S:2:N: m:S:-2:N:]","z [km]")

 plotb = xylineplot(wks,r0,z2,"D",0,"T",1.2,"C",1)

 set_size(plotb,xf(1),yf,wi,he)
 axis_range(plotb,0,0.3,-999,45,False,False)
 axis_str(plotb,"","R","")


 draw(plot)
 draw(plotb)


 frame(wks)

end

