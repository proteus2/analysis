load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 wks = cr_wks("du_z",600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 10

 nyear = 32
 y1    = 0

 mstr1 = "Unresolved forcing (MERRA) ("+lat_r+"N-"+lat_r+"S)"

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/tropics/dat"
 file_dat1 = dir_dat+"/merra.dudtgwd_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/merra.dudtana_ypt."+year1+"-"+year2+".nc"
 file_dat0 = dir_dat+"/merra.h_ypt."+year1+"-"+year2+".nc"
 file_dat3 = dir_dat+"/ng_0.nc"
 file_dat4 = dir_dat+"/cg_n.nc"
 file_dat5 = dir_dat+"/cg_s.nc"

 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f0 = addfile(file_dat0,"r")

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 u1i = f1->dudtgwd(m_1:m_2,:,{-lat_r:lat_r}) + f2->dudtana(m_1:m_2,:,{-lat_r:lat_r})
 u0i = f0->h(m_1:m_2,:,{-lat_r:lat_r})/1.e3

 lat = f1->latitude({-lat_r:lat_r})
 p1  = f1->levels

 np1 = dimsizes(p1)
 nmon = dimsizes(f1->time(m_1:m_2))

 nz = 66
 z = fspan(0.,65.,nz)
 u1ii = new((/nz,nmon,dimsizes(lat)/),float)
 do n=0, nmon-1
 do j=0, dimsizes(lat)-1
   u1ii(:,n,j) = ftcurv(u0i(n,:,j),u1i(n,:,j),z)
 end do
 end do
 u1iii = dim_avg( u1ii ) 

 z!0 = "z"
 z&z = z

 u = new((/1,2,nz/),float,1.e32)
 u = 0.
 u!2 = "z"
 u&z = z

 do k=0, nz-1
 do n=0, nmon-1
   u(0,0,k) = u(0,0,k) + max( (/0.,u1iii(k,n)/) )
   u(0,1,k) = u(0,1,k) + min( (/0.,u1iii(k,n)/) )
 end do
 end do
 u = u/nmon * 86400.
 
;--------------------------------------- process 0

;--------------------------------------- line attribute
 ldsh = (/0,0,11/)
 ltck = (/2.0,0.8,0.8/)
 lcol = (/18,6/)

 ymin = (/40,15/)
 ymax = (/max(z),45/)
 xmin = (/-1.2,-0.2/)
 xmax = (/1.2,0.3/)
 mstr = (/mstr1,""/)
;--------------------------------------- plot
 wi = 50
 he = 48
 xf = 20 + (/0,1/)*65
 yf = 10 + (/0,1/)*60

 plot = new(2,graphic)

 do ip=0, 1

 plot(ip) = xylineplot(wks,0.,z,"D",0,"T",0.5,"C",1)
 do is=0, 1
 do il=0, 0;2
   xyaddline(plot(ip),u(il,is,:),z,"D",ldsh(il) ,"T",ltck(il),"C",lcol(is))
 end do
 end do

 set_size(plot(ip),xf(0),yf(ip),wi,he)
 axis_range(plot(ip),xmin(ip),xmax(ip),ymin(ip),ymax(ip),False,False)
 axis_str(plot(ip),mstr(ip),"zonal forcing (m s:S:-1:N: day:S:-1:N:)","log-p height (km)")

 end do

 draw(plot)


 frame(wks)

end

