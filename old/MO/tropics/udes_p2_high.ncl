load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+144
 lat0 = 5.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 m2 = 0*12
 m3 = 0*12
 nmon0 = nmon-4

 wks = cr_wks("fig/udes_p2_high-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,:,{-lat0:lat0},0) )
   lat1 = f1->latitude({-lat0:lat0})
   p1 = f1->p
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )

   ny = dimsizes(lat1)
   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
   u4 = dim_avg( f3->u(m2:m2+nmon0-1,:,{-lat0:lat0}) )
   u5 = dim_avg( f3->u(m3:m3+nmon0-1,:,{-lat0:lat0}) )
   lat = f3->latitude({-lat0:lat0})
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   ny3 = dimsizes(lat)
   nz3 = dimsizes(p)
   nt3 = nmon0

   u = new((/2,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:)/)
   end do

   ur = new((/3,nz3,nt3/),float)
   do n=0, nt3-1
     ur(0,:,n) = u3(n,:)
     ur(1,:,n) = u4(n,:)
     ur(2,:,n) = u5(n,:)
   end do

   t = fspan(1,nt,nt)
   t3 = fspan(1,nt3,nt3)

   uss = new((/2,nz,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

   uss0 = new((/3,nz3,12/),float)
   do n=0, 11
     uss0(:,:,n) = dim_avg(ur(:,:,ispan(n,nt3-1,12)))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,:,n) = ur(:,:,n) - uss0(:,:,mod(n,12))
   end do

   ifilt = nt/12-1
   ifilt0 = nt3/12-1

   fftu = ezfftf(udes)
   fftu(:,:,:,ifilt+1:) = 0.
   urun = ezfftb(fftu,0.)

   fftu0 = ezfftf(udes0)
   fftu0(:,:,:,ifilt0+1:) = 0.
   urun0 = ezfftb(fftu0,0.)

   uano  = udes  - urun
   uano0 = udes0 - urun0

   uu = (/udes,urun,uano/)
   uu@_FillValue = 1.e32
   uu0 = (/udes0,urun0,uano0/)
   uu0@_FillValue = 1.e32

   usd  = dim_stddev(uu)
   usd0 = dim_stddev(uu0)

;--------------------------------------- process 0

;--------------------------------------- line att.
 llbl = (/"CTL","CC","MERRA"/)
 lthc = (/1.,1.,2./)
 ldsh = (/14,0,0/)

;--------------------------------------- contour
 wi = 35
 he = 36
 xf = 12 + (/0,1,0/)*43
 yf = 10 + (/0,1,2,3,4/)*42

 mstr = (/"Total","Low Freq.","High Freq."/)
 xstr = "Std. Dev. (m s:S:-1:N:)"
 ystr = (/"log-p height (km)","",""/)
 xmax = (/28,28,15/)


do isc=2,2;0, 2

 plot = xylineplot(wks,usd(isc,:,:),zp1,"D",ldsh(0:1),"T",lthc(0:1),"C",1)
 xyaddline(plot,dim_avg_n(usd0(isc,:,:),0),zp,"D",ldsh(2),"T",lthc(2),"C",1)

 set_size(plot,xf(isc),yf,wi,he)
 axis_range(plot,0,xmax(isc),15,65,False,False)
 axis_str(plot,mstr(isc),xstr,ystr(0))

 draw(plot)

end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",lthc,"C",1,2)
 set_size_lg(lg,xf+wi*0.55,yf(0)+he*0.7,wi*0.25,he*0.3)
 draw(lg)


 frame(wks)

end



