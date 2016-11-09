load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 wks = cr_wks("du_p",600,0)
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

 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")

; lat = f1->lat({-lat_r:lat_r})
 p1  = f1->levels

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 u1i = dim_avg( f1->dudtgwd(m_1:m_2,:,{-lat_r:lat_r}) )
 u2i = dim_avg( f2->dudtana(m_1:m_2,:,{-lat_r:lat_r}) )
 u0i = u1i + u2i

 np1 = dimsizes(p1)
 nmon = dimsizes(f1->time(m_1:m_2))

 zp1 = 7.0*log(1.e3/p1)
 zp1!0 = "zp"
 zp1&zp = zp1

 u = new((/3,2,np1/),float,1.e32)
 u = 0.
 u!2 = "zp"
 u&zp = zp1

 do k=0, np1-1
 do n=0, nmon-1
   u(0,0,k) = u(0,0,k) + max( (/0.,u0i(n,k)/) )
   u(1,0,k) = u(1,0,k) + max( (/0.,u1i(n,k)/) )
   u(2,0,k) = u(2,0,k) + max( (/0.,u2i(n,k)/) )
   u(0,1,k) = u(0,1,k) + min( (/0.,u0i(n,k)/) )
   u(1,1,k) = u(1,1,k) + min( (/0.,u1i(n,k)/) )
   u(2,1,k) = u(2,1,k) + min( (/0.,u2i(n,k)/) )
 end do
 end do
 u = u/nmon * 86400.
 
;--------------------------------------- process 0

;--------------------------------------- line attribute
 ldsh = (/0,0,11/)
 ltck = (/2.0,0.8,0.8/)
 lcol = (/18,6/)

 ymin = (/40,15/)
 ymax = (/max(zp1),45/)
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

 plot(ip) = xylineplot(wks,0.,zp1,"D",0,"T",0.5,"C",1)
 do is=0, 1
 do il=0, 2
   xyaddline(plot(ip),u(il,is,:),zp1,"D",ldsh(il) ,"T",ltck(il),"C",lcol(is))
 end do
 end do

 set_size(plot(ip),xf(0),yf(ip),wi,he)
 axis_range(plot(ip),xmin(ip),xmax(ip),ymin(ip),ymax(ip),False,False)
 axis_str(plot(ip),mstr(ip),"zonal forcing (m s:S:-1:N: day:S:-1:N:)","log-p height (km)")

 end do

 draw(plot)


 frame(wks)

end

