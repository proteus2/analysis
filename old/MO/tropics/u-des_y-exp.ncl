load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 expn = "0"

 wks = cr_wks("u-des_y-"+expn,600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 30
 p     = (/10,30,50/)

 nyear = 32
 y1    = 0

 mstr1 = "U (UM-CTL; ERA-I; MERRA) ("+lat_r+"N-"+lat_r+"S)"

 mon_thr = 6

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/tropics/dat"
 file_dat1 = dir_dat+"/era-int.u_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat3 = dir_dat+"/uh-p_"+expn+".nc"
 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f3 = addfile(file_dat3,"r")

 lat1 = f1->latitude({-lat_r:lat_r})
 lat2 = f2->latitude({-lat_r:lat_r})
 lat3 = f3->latitude({-lat_r:lat_r})
 p1  = f1->levelist
 p2  = f2->levels
 p3  = f3->p

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1 

 u1i = doubletofloat( f1->u(m_1:m_2,:,{-lat_r:lat_r}) )
 u2i = f2->u(m_1:m_2,:,{-lat_r:lat_r})
 u3i = f3->u(:,:,{-lat_r:lat_r},0)
 u1i!1 = "p"
 u1i&p = p1

 ny1 = dimsizes(lat1)
 ny2 = dimsizes(lat2)
 ny3 = dimsizes(lat3)
 nz  = dimsizes(p)
 nmon_c = dimsizes(f3->t(:))

 u1 = new((/2,nz,ny1/),float,1.e32)
 u2 = new((/2,nz,ny2/),float,1.e32)
 u3 = new((/2,nz,ny3/),float,1.e32)
 u1 = 0.
 u2 = 0.
 u3 = 0.

 do n=0, nmon-1
   u1(0,:,:) = u1(0,:,:) + u1i(n,{p(:)},:)
   u2(0,:,:) = u2(0,:,:) + u2i(n,{p(:)},:)
 end do
 do n=0, nmon_c-1
   u3(0,:,:) = u3(0,:,:) + u3i(n,{p(:)},:)
 end do
 u1(0,:,:) = u1(0,:,:)/nmon
 u2(0,:,:) = u2(0,:,:)/nmon
 u3(0,:,:) = u3(0,:,:)/nmon_c

 u1ss = new((/12,nz,ny1/),float)
 u2ss = new((/12,nz,ny2/),float)
 u3ss = new((/12,nz,ny3/),float)
 do n=0, 11
 do k=0, nz-1
   do j=0, ny1-1
     u1ss(n,k,j) = avg( u1i(ispan(n,nmon-1,12),{p(k)},j) )
   end do
   do j=0, ny2-1
     u2ss(n,k,j) = avg( u2i(ispan(n,nmon-1,12),{p(k)},j) )
   end do
   do j=0, ny3-1
     u3ss(n,k,j) = avg( u3i(ispan(n,nmon_c-1,12),{p(k)},j) )
   end do
 end do
 end do

 u1des = new((/nmon,nz,ny1/),float)
 u2des = new((/nmon,nz,ny2/),float)
 u3des = new((/nmon_c,nz,ny3/),float)
 do n=0, nmon-1
   u1des(n,:,:) = u1i(n,{p(:)},:) - u1ss(mod(n,12),:,:)
   u2des(n,:,:) = u2i(n,{p(:)},:) - u2ss(mod(n,12),:,:)
 end do
 do n=0, nmon_c-1
   u3des(n,:,:) = u3i(n,{p(:)},:) - u3ss(mod(n,12),:,:)
 end do

 do n=0, nmon-1
   u1(1,:,:) = u1(1,:,:) + u1des(n,:,:)^2
   u2(1,:,:) = u2(1,:,:) + u2des(n,:,:)^2
 end do
 do n=0, nmon_c-1
   u3(1,:,:) = u3(1,:,:) + u3des(n,:,:)^2
 end do
 u1(1,:,:) = sqrt(u1(1,:,:)/nmon)
 u2(1,:,:) = sqrt(u2(1,:,:)/nmon)
 u3(1,:,:) = sqrt(u3(1,:,:)/nmon_c)

;--------------------------------------- process 0

;--------------------------------------- line attribute

;--------------------------------------- plot
 wi = 36
 he = 30
 xf = 20 + (/0,1,2/)*50
 yf = 10 + (/0,1,2/)*40

 do k=0, nz-1

 plot = xylineplot(wks,lat1,0.,"D",0,"T",0.5,"C",colind)
 xyaddline(plot,lat1,u1(0,k,:),"D",11,"T",0.8,"C",1)
 xyaddline(plot,lat2,u2(0,k,:),"D",0 ,"T",0.8,"C",1)
 xyaddline(plot,lat3,u3(0,k,:),"D",0 ,"T",2.0,"C",1)

 set_size(plot,xf(0),yf(k),wi,he)
 axis_range(plot,-30,30,-18,3,False,False)
 axis_str(plot,"","lat (deg)","mean (m s:S:-1:N:)")

 draw(plot)


 plot = xylineplot(wks,lat1,0.,"D",0,"T",0.5,"C",1)
 xyaddline(plot,lat1,u1(1,k,:),"D",11,"T",0.8,"C",1)
 xyaddline(plot,lat2,u2(1,k,:),"D",0 ,"T",0.8,"C",1)
 xyaddline(plot,lat3,u3(1,k,:),"D",0 ,"T",2.0,"C",1)

 set_size(plot,xf(1),yf(k),wi,he)
 axis_range(plot,-30,30,0,21,False,False)
 axis_str(plot,"","lat (deg)","amplitude (m s:S:-1:N:)")

 draw(plot)

 dr_txt(wks,p(k)+" hPa",xf(0)+6,yf(k)+he-3,22,6,0.)

 end do


 frame(wks)

end

