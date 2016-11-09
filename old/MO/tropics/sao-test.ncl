load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 expn = (/"0","v","f","d","u"/)  ; bl/ye/gr/re

 wks = cr_wks("sao-test",600,0)
 gsn_define_colormap(wks,"amwg")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 5
 p = (/0.3,1.,3./)

 nyear = 32
 y1    = 0

 mstr1 = "U (UM-CTL; ERA-I; MERRA) ("+lat_r+"N-"+lat_r+"S)"

;--------------------------------------- read data
 dir_dat   = "../dat"
 file_dat1 = dir_dat+"/era-int.u_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat3 = dir_dat+"/uh-p_"+expn+".nc"
 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f3 = addfile(file_dat3(0),"r")

 nexp = dimsizes(expn)

; lat = f1->lat({-lat_r:lat_r})
 p1  = f1->levelist
 p2  = f2->levels
 p3  = f3->p

 np = dimsizes(p)
 kp = new((/np,3/),integer)
 do k=0, np-1
   kp(k,0) = minind(abs(log(p(k)/p1)))
   kp(k,1) = minind(abs(log(p(k)/p2)))
   kp(k,2) = minind(abs(log(p(k)/p3)))
 end do
print(p1(kp(:,0)))
print(p2(kp(:,1)))
print(p3(kp(:,2)))

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1 

 u1i = doubletofloat(dim_avg( f1->u(m_1:m_2,:,{-lat_r:lat_r}) ))
 u2i = dim_avg( f2->u(m_1:m_2,:,{-lat_r:lat_r}) )

 np1 = dimsizes(p1)
 np2 = dimsizes(p2)
 np3 = dimsizes(p3)

 nmon_c = 36 ;dimsizes(f3->t(4:))
 u3i = new((/nexp,nmon_c,np3/),float)
 do ie=0, nexp-1
   f3 = addfile(file_dat3(ie),"r")
   u3i(ie,:,:) = dim_avg( f3->u(4:nmon_c+3,:,{-lat_r:lat_r},0) )
 end do

 zp1 = 7.0*log(1.e3/p1)
 zp2 = 7.0*log(1.e3/p2)
 zp3 = 7.0*log(1.e3/p3)

 u1ss = new((/14,np1/),float)
 u2ss = new((/14,np2/),float)
 u3ss = new((/nexp,14,np3/),float)
 do n=1, 12
   u1ss(n,:) = dim_avg_n( u1i(ispan(n-1,nmon-1,12),:),0 )
   u2ss(n,:) = dim_avg_n( u2i(ispan(n-1,nmon-1,12),:),0 )
   do ie=0, nexp-1
     u3ss(ie,n,:) = dim_avg_n( u3i(ie,ispan(n-1,nmon_c-1,12),:),0 )
   end do
 end do
 u1ss(0,:) = u1ss(12,:)
 u2ss(0,:) = u2ss(12,:)
 u3ss(:,0,:) = u3ss(:,12,:)
 u1ss(13,:) = u1ss(1,:)
 u2ss(13,:) = u2ss(1,:)
 u3ss(:,13,:) = u3ss(:,1,:)

 u1sd = new((/14,np1/),float)
 u2sd = new((/14,np2/),float)
 u3sd = new((/nexp,14,np3/),float)
 do n=1, 12
   u1sd(n,:) = dim_stddev_n( u1i(ispan(n-1,nmon-1,12),:),0 )
   u2sd(n,:) = dim_stddev_n( u2i(ispan(n-1,nmon-1,12),:),0 )
   do ie=0, nexp-1
   u3sd(ie,n,:) = dim_stddev_n( u3i(ie,ispan(n-1,nmon_c-1,12),:),0 )
   end do
 end do
 u1sd(0,:) = u1sd(12,:)
 u2sd(0,:) = u2sd(12,:)
 u3sd(:,0,:) = u3sd(:,12,:)
 u1sd(13,:) = u1sd(1,:)
 u2sd(13,:) = u2sd(1,:)
 u3sd(:,13,:) = u3sd(:,1,:)

;--------------------------------------- process 0

;--------------------------------------- line attribute
 lcol = (/1,4,12,8,17/)

;--------------------------------------- plot
 wi = 56
 he = 50
 xf = 20 + (/0,1,2/)*70
 yf = 10 + (/0,1,2/)*60

 mon = ispan(0,13,1)

 do k=0, np-1

   plot = xylineplot(wks,(/0,13/),0.,"D",0,"T",0.5,"C",colind)
   do ie=0, nexp-1
     xyaddline(plot,mon,u3ss(ie,:,kp(k,2)),"D",0 ,"T",2.0,"C",lcol(ie))
   end do
   if (p(k) .ge. 1.) then
     xyaddline(plot,mon,u1ss(:,kp(k,0)),"D",11,"T",0.8,"C",1)
   end if
   xyaddline(plot,mon,u2ss(:,kp(k,1)),"D",0 ,"T",0.8,"C",1)

   set_size(plot,xf(0),yf(k),wi,he)
   axis_range(plot,0.5,12.5,-60,40,False,False)
   axis_str(plot,"","month","U (m s:S:-1:N:)")
   setvalues plot
    "tmXBMode"   : "Explicit"
    "tmXBValues" : mon
    "tmXBLabels" : mon
   end setvalues

   draw(plot)

   plot = xylineplot(wks,(/0,13/),0.,"D",0,"T",0.5,"C",colind)
   do ie=0, nexp-1
     xyaddline(plot,mon,u3sd(ie,:,kp(k,2)),"D",0 ,"T",2.0,"C",lcol(ie))
   end do
   if (p(k) .ge. 1.) then
     xyaddline(plot,mon,u1sd(:,kp(k,0)),"D",11,"T",0.8,"C",1)
   end if
   xyaddline(plot,mon,u2sd(:,kp(k,1)),"D",0 ,"T",0.8,"C",1)

   set_size(plot,xf(1),yf(k),wi,he)
   axis_range(plot,0.5,12.5,0,25,False,False)
   axis_str(plot,"","month","SD_U (m s:S:-1:N:)")
   setvalues plot
    "tmXBMode"   : "Explicit"
    "tmXBValues" : mon
    "tmXBLabels" : mon
   end setvalues

   draw(plot)

 end do


 frame(wks)

end

