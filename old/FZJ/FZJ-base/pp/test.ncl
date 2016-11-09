load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

cx = 10 ;500.e3 / 

begin

 wks = cr_wks("test",600,0)
 gsn_define_colormap(wks,"gsdtol")
 icol_center = 11
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
 pi = 3.141592

 f0 = addfile("../res/uvn_background_p1_region.nc","r")
 f1 = addfile("../res/t_profile_ave_p1.nc","r")
 f2 = addfile("../res/uvn_background2_p1_region.nc","r")

 z  = f0->height_m({20000:}) / 1.e3
 ub = f0->ub(:,{20000:},:,:)
 nb = f0->nb(:,{20000:},:,:)

 nz = dimsizes(z)

 tbavg = f1->t_ave({20000:})
 ubavg = f2->ubavg({20000:})
 vbavg = f2->vbavg({20000:})

 u0 = new(nz,float,1.e32)
 n0 = new(nz,float,1.e32)
 usd   = new(nz,float,1.e32)
 polyu = new(nz*2,float,1.e32)
 polyuz= polyu
 do k=0, nz-1
   u0(k) = avg(ub(:,k,:,:))
   n0(k) = avg(nb(:,k,:,:))
   usd   (k) = sqrt( avg((ub(:,k,:,:)-u0(k))^2) )
   polyu (k       ) = u0(k) - usd(k)
   polyu (2*nz-1-k) = u0(k) + usd(k)
   polyuz(k       ) = z(k)
   polyuz(2*nz-1-k) = z(k)
 end do

 m0 = abs(n0/(cx-u0))
 m0i= 1./m0
 l0 = 2.*pi/m0

 mb = abs(nb/(cx-ub))
 mbi= 1./mb
 lb = 2.*pi/mb

 mb2 = mb^2

 mavg  = new(nz,float,1.e32)
 msd   = new(nz,float,1.e32)
 polym = new(nz*2,float,1.e32)
 polymz= polym
 do k=0, nz-1
   mavg  (k) = avg(mb2(:,k,:,:))
   msd   (k) = sqrt( avg((mb2(:,k,:,:)-mavg(k))^2) )
   polym (k       ) = mavg(k) - msd(k)
   polym (2*nz-1-k) = mavg(k) + msd(k)
   polymz(k       ) = z(k)
   polymz(2*nz-1-k) = z(k)
 end do

 cgz =  n0*2.*pi/500.e3/m0^2

 tlast = new(nz,float,1.e32)
 tlast(0) = 0.
 do k=1, nz-1
   tlast(k) = tlast(k-1) + 2./(cgz(k-1)+cgz(k)) * (z(k)-z(k-1))*1.e3 / 3600.
 end do

print(z(:)+"   "+mavg(:))

 wl = wavelet_default(dtrend_quadratic(vbavg,0),0)
; wl = wavelet_default(dtrend(ubavg,0),0)

 wla = wl(0,:,:)^2 + wl(1,:,:)^2

 m2t = (2.*pi/(wl@period*500.))^2
print(dtrend(vbavg,0))

;--------------------------------------- line attribute
 cnlev1 = (/0.05,0.1,0.15,0.2,0.3,0.4,0.6/)
 cnfil1 = (/0,28,24,20,16,12,8,4/)
;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 10 + (/0,1,2/)*50
 yf = 10

 mstr = (/"","MFx_T_up","T variance"/)


 plot = xylineplot(wks,u0,z,"D",0,"T",1.2,"C",1)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,-50,-10,-999,45,False,False)
 axis_str(plot,"","U [m s:S:-1:N:]","z [km]")
 setvalues plot
  "tmXTOn" : False
 end setvalues

 plot2 = xylineplot(wks,n0*100.,z,"D",1,"T",1.2,"C",1)

 set_size(plot2,xf(0),yf,wi,he)
 axis_range(plot2,1.5,3,-999,45,False,False)
 axis_str(plot2,"","N [10:S:-2:N: s:S:-1:N:]","")
 setvalues plot2
  "tmXBOn" : False
  "tmXTLabelsOn" : True
  "tiXAxisSide" : "Top"
 end setvalues

 plres = True
 plres@gsFillColor = colind
 dummy = gsn_add_polygon(wks,plot,polyu,polyuz,plres)


 plotb = xylineplot(wks,mavg,z,"D",0,"T",1.2,"C",1)

 set_size(plotb,xf(1),yf,wi,he)
 axis_range(plotb,1.e-7,2.e-6,-999,45,True,False)
 axis_str(plotb,"","m:S:2:N: [rad:S:2:N: m:S:-2:N:]","")

 plres = True
 plres@gsFillColor = colind
 dummyb = gsn_add_polygon(wks,plotb,polym,polymz,plres)

 draw(plot)
 draw(plot2)
 draw(plotb)


 contour = cnshadeplot(wks,z,m2t,wla, \  ; 2.e2
           cnlev1,cnfil1,0)

 set_size(contour,xf(2),yf,wi,he)
 axis_range(contour,-999,999,1.e-7,4.e-6,False,True)
 axis_str(contour,"","x [km]","")
print(wla)

 cnfinalize(contour,wks)

 draw(contour)


dr_txt(wks,"(a)",xf(0)-2,yf-7,25,8,0.)
dr_txt(wks,"(b)",xf(1)-2,yf-7,25,8,0.)


 frame(wks)

end

