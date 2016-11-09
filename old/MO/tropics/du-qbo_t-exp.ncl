load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 expn = "0"

 wks = cr_wks("du-qbo_t-"+expn,600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 5;10

 nyear = 6
 y1    = 0

 p = (/10,30,50/)
 zp = 7.0e3*log(1000./p)

 mstr1 = "Zonal forcing (UM-CTL, MERRA, monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/tropics/dat"
 file_dat0 = dir_dat+"/merra.h_ypt."+year1+"-"+year2+".nc"
 file_dat1 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat2a = dir_dat+"/merra.dudtgwd_ypt."+year1+"-"+year2+".nc"
 file_dat2b = dir_dat+"/merra.dudtana_ypt."+year1+"-"+year2+".nc"
 file_dat3 = dir_dat+"/u_0.nc"
 file_dat4 = dir_dat+"/ng_0.nc"
 file_dat5 = dir_dat+"/cg_n.nc"
 file_dat6 = dir_dat+"/cg_s.nc"

 f0 = addfile(file_dat0,"r")
 f1 = addfile(file_dat1,"r")
 f2a= addfile(file_dat2a,"r")
 f2b= addfile(file_dat2b,"r")
 f3 = addfile(file_dat3,"r")
 f4 = addfile(file_dat4,"r")
 f5 = addfile(file_dat5,"r")
 f6 = addfile(file_dat6,"r")

; lat = f1->lat({-lat_r:lat_r})
 p2  = f2a->levels

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 u0i = dim_avg_Wrap( f0->h       (m_1:m_2,:,{-lat_r:lat_r}) )
 u1i = dim_avg_Wrap( f1->u       (m_1:m_2,:,{-lat_r:lat_r}) )
 u2i = dim_avg_Wrap( f2a->dudtgwd(m_1:m_2,:,{-lat_r:lat_r}) + \
                     f2b->dudtana(m_1:m_2,:,{-lat_r:lat_r}) )
 u3i = dim_avg_Wrap( f3->u          (m_1+4:m_2+4,:,{-lat_r:lat_r},0) )
 u4i = dim_avg_Wrap( f4->field424   (m_1+4:m_2+4,:,{-lat_r:lat_r},0) )
 u5i = dim_avg_Wrap( f5->unspecified(m_1+4:m_2+4,:,{-lat_r:lat_r},0) )
 u6i = dim_avg_Wrap( f6->unspecified(m_1+4:m_2+4,:,{-lat_r:lat_r},0) )

 nz = dimsizes(p)
 nmon = dimsizes(f1->time(m_1:m_2))

 mon = fspan(0,nmon-1,nmon)
 nl = 6
 u = new((/nz,nl,nmon/),float)
 do k=0, nz-1
   u(k,0,:) = ftcurv(u0i,u1i,zp(k))
   u(k,1,:) = ftcurv(u0i,u2i,zp(k))*86400. *50.
   u(k,2,:) = (/ u3i(:,{zp(k)}) /)
   u(k,3,:) = (/ u4i(:,{zp(k)}) /)*86400. *50.
   u(k,4,:) = (/ u5i(:,{zp(k)}) /)*86400. *100.
   u(k,5,:) = (/ u6i(:,{zp(k)}) /)*86400. *100.
 end do
 u0 = dim_avg(u)

;--------------------------------------- process 0

;--------------------------------------- line attribute
 ldsh = (/1,0,1,0,0,0/)
 ltck = (/0.8,1.2,0.8,1.2,1.2,0.5/)
; lcol = (/1,1,6,6,18,18/)
 lcol = (/1,0,6,6,18,18/)

 ymin = (/-50,-40,-30/)
 ymax = (/30,20,20/)
 mstr = (/mstr1,"",""/)

;--------------------------------------- plot
 wi = 132
 he = 32
 xf = 20
 yf = 10 + (/0,1,2/)*45

  plot = new(nz,graphic)

 do ip=0, nz-1

 plot(ip) = xylineplot(wks,mon,0.,"D",0,"T",0.5,"C",1)
 do il=0, nl-1;-2
;   xyaddline(plot(ip),mon,u0(ip,il),"D",0,"T",0.5,"C",lcol(il))
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

