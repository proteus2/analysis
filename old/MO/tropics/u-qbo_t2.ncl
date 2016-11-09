load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 wks = cr_wks("u-qbo_t2",600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 10

 nyear = 32
 y1    = 0

 p = (/10,30,50/)

 mstr1 = "U (ERA-I, MERRA, monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"

;--------------------------------------- read data
 dir_dat   = "../dat"
 file_dat1 = dir_dat+"/era-int.u_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"

 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")

; lat = f1->lat({-lat_r:lat_r})
 p1  = f1->levelist
 p2  = f2->levels

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 u1i = doubletofloat(dim_avg( f1->u(m_1:m_2,:,{-lat_r:lat_r}) ))
 u2i = dim_avg_Wrap( f2->u(m_1:m_2,:,{-lat_r:lat_r}) )
 u1i!1 = "p"
 u1i&p = p1

 nz = dimsizes(p)
 nmon = dimsizes(f1->time(m_1:m_2))

 mon = fspan(0,nmon-1,nmon)
 nl = 2
 u = new((/nz,nl,nmon/),float)
 do k=0, nz-1
   u(k,0,:) = u1i(:,{p(k)})
   u(k,1,:) = (/ u2i(:,{p(k)}) /)
 end do
 u0 = dim_avg(u)

;--------------------------------------- process 0

;--------------------------------------- line attribute
 ldsh = (/0,0,0/)
 ltck = (/1.2,1.2,1.2/)
 lcol = (/6,18/)

 ymin = (/-50,-40,-30/)
 ymax = (/30,20,20/)
 mstr = (/mstr1,"",""/)

;--------------------------------------- plot
 wi = 132
 he = 12
 xf = 20
 yf = 10 + (/0,1,2/)*25

  plot = new(nz,graphic)

 do ip=0, nz-1

 plot(ip) = xylineplot(wks,mon,0.,"D",0,"T",0.5,"C",1)
 do il=0, nl-1
   xyaddline(plot(ip),mon,u0(ip,il),"D",0,"T",0.5,"C",lcol(il))
   xyaddline(plot(ip),mon,u(ip,il,:),"D",ldsh(il) ,"T",ltck(il),"C",lcol(il))
 end do

 set_size(plot(ip),xf(0),yf(ip),wi,he)
 axis_range(plot(ip),0,mon(nmon-1),ymin(ip),ymax(ip),False,False)
 axis_str(plot(ip),mstr(ip),"time (month)","U (m s:S:-1:N:)")

 dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 end do

 draw(plot)


 frame(wks)

end

