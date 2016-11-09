load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 expn = (/"i","f","fd"/)

 wks = cr_wks("u-qbo_t2-"+expn(0)+expn(1)+expn(2),600,0)
 gsn_define_colormap(wks,"amwg")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 5

 nyear = 10
 y1    = 0

 p = (/10,30,50/)

 mstr1 = "U (monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"

;--------------------------------------- read data
 dir_dat   = "../dat"
 file_dat1 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/uh-p_"+expn(0)+".nc"
 file_dat3 = dir_dat+"/uh-p_"+expn(1)+".nc"
 file_dat4 = dir_dat+"/uh-p_"+expn(2)+".nc"

 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f3 = addfile(file_dat3,"r")
 f4 = addfile(file_dat4,"r")

 p1  = f1->levels
 p2  = f2->p

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 u1i = dim_avg_Wrap( f1->u(m_1:m_2,:,{-lat_r:lat_r}) )
 nmon_c_all = (/ dimsizes( f2->t(4+m_1:) ), \
                 dimsizes( f3->t(4+m_1:) ), \
                 dimsizes( f4->t(4+m_1:) ) /)

 if (nmon .le. nmon_c_all(0)) then
   u2i = dim_avg_Wrap( f2->u(4+m_1:4+m_2,:,{-lat_r:lat_r},0) )
 else
   u2i0 = dim_avg_Wrap( f2->u(4+m_1:,:,{-lat_r:lat_r},0) )
   arrapp2 = new((/nmon-nmon_c_all(0),dimsizes(p2)/),float,1.e32)
   arrapp2 = 1.e32
   u2i = array_append_record(u2i0,arrapp2,0)
 end if
 if (nmon .le. nmon_c_all(1)) then
   u3i = dim_avg_Wrap( f3->u(4+m_1:4+m_2,:,{-lat_r:lat_r},0) )
 else
   u3i0 = dim_avg_Wrap( f3->u(4+m_1:,:,{-lat_r:lat_r},0) )
   arrapp3 = new((/nmon-nmon_c_all(1),dimsizes(p2)/),float,1.e32)
   arrapp3 = 1.e32
   u3i = array_append_record(u3i0,arrapp3,0)
 end if
 if (nmon .le. nmon_c_all(2)) then
   u4i = dim_avg_Wrap( f4->u(4+m_1:4+m_2,:,{-lat_r:lat_r},0) )
 else
   u4i0 = dim_avg_Wrap( f4->u(4+m_1:,:,{-lat_r:lat_r},0) )
   arrapp4 = new((/nmon-nmon_c_all(2),dimsizes(p2)/),float,1.e32)
   arrapp4 = 1.e32
   u4i = array_append_record(u4i0,arrapp4,0)
 end if

 nz = dimsizes(p)

 mon = fspan(0,nmon-1,nmon)
 nl = 4
 u = new((/nz,nl,nmon/),float,1.e32)
 do k=0, nz-1
   u(k,0,:) = (/ u1i(:,{p(k)}) /)
   u(k,1,:) = (/ u2i(:,{p(k)}) /)
   u(k,2,:) = (/ u3i(:,{p(k)}) /)
   u(k,3,:) = (/ u4i(:,{p(k)}) /)
 end do
 u0 = dim_avg(u)

;--------------------------------------- process 0

;--------------------------------------- line attribute
 ldsh = (/1,0,0,0/)
 ltck = (/1.5,1.2,1.2,1.2/)
 lcol = (/1,9,3,14/)

 ymin = (/-50,-40,-30/)
 ymax = (/30,30,30/)
 mstr = (/mstr1,"",""/)

;--------------------------------------- plot
 wi = 132
 he = 32
 xf = 20
 yf = 10 + (/0,1,2/)*45

  plot = new(nz,graphic)

 do ip=0, nz-1

 plot(ip) = xylineplot(wks,mon,0.,"D",0,"T",0.5,"C",1)
 do il=0, nl-1
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

