load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 opt_mean = 2   ;  0 - raw  ;  1 - anomaly  ;  2 - deseasonalized anomaly

 wks = cr_wks("u-des_p",600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 10

 nyear = 32
 y1    = 0

 mstr1 = "U (ERA-I; MERRA) ("+lat_r+"N-"+lat_r+"S)"

 mon_thr = 6

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/tropics/dat"
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
 u2i = dim_avg( f2->u(m_1:m_2,:,{-lat_r:lat_r}) )

 np1 = dimsizes(p1)
 np2 = dimsizes(p2)

 zp1 = 7.0*log(1.e3/p1)
 zp2 = 7.0*log(1.e3/p2)

 u1 = new((/5,np1/),float,1.e32)
 u2 = new((/5,np2/),float,1.e32)
 u1 = 0.
 u2 = 0.

 do n=0, nmon-1
   u1(0,:) = u1(0,:) + u1i(n,:)
   u2(0,:) = u2(0,:) + u2i(n,:)
 end do
 u1(0,:) = u1(0,:)/nmon
 u2(0,:) = u2(0,:)/nmon

 u1ss = new((/12,np1/),float)
 u2ss = new((/12,np2/),float)
 u1ss = 0.
 u2ss = 0.
 if (opt_mean .eq. 1) then
   u1ss(0,:) = dim_avg_n(u1i,0)
   u2ss(0,:) = dim_avg_n(u2i,0)
 end if
 if (opt_mean .eq. 2) then
   do n=0, 11
     u1ss(n,:) = dim_avg_n( u1i(ispan(n,nmon-1,12),:),0 )
     u2ss(n,:) = dim_avg_n( u2i(ispan(n,nmon-1,12),:),0 )
   end do
 end if

 u1pn = new((/2,nmon,np1/),float)
 u2pn = new((/2,nmon,np2/),float)
 do n=0, nmon-1
   u1pn(0,n,:) = u1i(n,:) - u1ss(mod(n,12),:)
   u2pn(0,n,:) = u2i(n,:) - u2ss(mod(n,12),:)
 end do
 u1pn(1,:,:) = u1pn(0,:,:)
 u2pn(1,:,:) = u2pn(0,:,:)

 do n=0, nmon-1
   do k=0, np1-1
     if (u1pn(0,n,k) .gt. 0.) then
       u1(1,k) = u1(1,k) + 1.
       u1pn(1,n,k) = 0.
     else
       u1pn(0,n,k) = 0.
     end if
   end do
   do k=0, np2-1
     if (u2pn(0,n,k) .gt. 0.) then
       u2(1,k) = u2(1,k) + 1.
       u2pn(1,n,k) = 0.
     else
       u2pn(0,n,k) = 0.
     end if
   end do
 end do
 u1(2,:) = nmon - u1(1,:)
 u2(2,:) = nmon - u2(1,:)

 do n=0, nmon-1
   u1(3:4,:) = u1(3:4,:) + u1pn(0:1,n,:)^2
   u2(3:4,:) = u2(3:4,:) + u2pn(0:1,n,:)^2
 end do
 do iv=3, 4
   do k=0, np1-1
     u1(iv,k) = sqrt( u1(iv,k)/max((/u1(iv-2,k),0.01/)) )
   end do
   do k=0, np2-1
     u2(iv,k) = sqrt( u2(iv,k)/max((/u2(iv-2,k),0.01/)) )
   end do
 end do

 u1(1:2,:) = 0.
 u2(1:2,:) = 0.
 stmp = (/"+","-"/)
 do j=0, 1
 do k=0, np1-1
   tmp1 = 0
   tmp2 = 0
   cnt = -999
   do n=0, nmon-1
     if (u1pn(j,n,k) .eq. 0.) then
       if (cnt .gt. mon_thr) then
         tmp1 = tmp1 + cnt
         tmp2 = tmp2 + 1
         if (p1(k) .eq. 30.) then
           print(stmp(j)+" phase duration at "+p1(k)+" hPa : "+cnt)
         end if
       end if
       cnt = 0
     else
       if (cnt .ge. 0) then
         cnt = cnt + 1
       end if
     end if
   end do
   if (tmp2 .gt. 0) then
     u1(j+1,k) = int2flt(tmp1)/tmp2
   else
     u1(j+1,k) = 1.e32
   end if
 end do
 end do
 do j=0, 1
 do k=0, np2-1
   tmp1 = 0
   tmp2 = 0
   cnt = -999
   do n=0, nmon-1
     if (u2pn(j,n,k) .eq. 0.) then
       if (cnt .gt. mon_thr) then
         tmp1 = tmp1 + cnt
         tmp2 = tmp2 + 1
         if (p2(k) .eq. 30.) then
           print(stmp(j)+" phase duration at "+p2(k)+" hPa : "+cnt)
         end if
       end if
       cnt = 0
     else
       if (cnt .ge. 0) then
         cnt = cnt + 1
       end if
     end if
   end do
   if (tmp2 .gt. 0) then
     u2(j+1,k) = int2flt(tmp1)/tmp2
   else
     u2(j+1,k) = 1.e32
   end if
 end do
 end do

;--------------------------------------- process 0

;--------------------------------------- line attribute
 cnlev = (/-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60/)*0.5
 cnfil = (/colind,3,4,5,6,7,8,9,0,0,14,15,16,17,18,19,20/)

;--------------------------------------- plot
 wi = 36
 he = 40
 xf = 20 + (/0,1,2/)*50
 yf = 10 + (/0,1/)*50

 plot = xylineplot(wks,0.,zp1,"D",0,"T",0.5,"C",colind)
 xyaddline(plot,u1(0,:),zp1,"D",11,"T",0.8,"C",1)
 xyaddline(plot,u2(0,:),zp2,"D",0 ,"T",0.8,"C",1)

 set_size(plot,xf(0),yf(0),wi,he)
 axis_range(plot,-20,5,15,47,False,False)
 axis_str(plot,"","mean (m s:S:-1:N:)","log-p height (km)")

 draw(plot)


 plot = xylineplot(wks,0.,zp1,"D",0,"T",0.5,"C",1)
 xyaddline(plot,u1(3,:),zp1,"D",11,"T",0.8,"C",18)
 xyaddline(plot,u2(3,:),zp2,"D",0 ,"T",0.8,"C",18)
 xyaddline(plot,u1(4,:),zp1,"D",11,"T",0.8,"C",6 )
 xyaddline(plot,u2(4,:),zp2,"D",0 ,"T",0.8,"C",6 )

 set_size(plot,xf(1),yf(0),wi,he)
 axis_range(plot,0,20,15,47,False,False)
 axis_str(plot,mstr1,"amplitude (m s:S:-1:N:)","log-p height (km)")

 draw(plot)


 plot = xylineplot(wks,mon_thr,zp1,"D",0,"T",0.5,"C",1)
 xyaddline(plot,u1(1,:),zp1,"D",11,"T",0.8,"C",18)
 xyaddline(plot,u2(1,:),zp2,"D",0 ,"T",0.8,"C",18)
 xyaddline(plot,u1(2,:),zp1,"D",11,"T",0.8,"C",6 )
 xyaddline(plot,u2(2,:),zp2,"D",0 ,"T",0.8,"C",6 )

 set_size(plot,xf(2),yf(0),wi,he)
 axis_range(plot,mon_thr,19,15,47,False,False)
 axis_str(plot,"","duration (month)","log-p height (km)")

 draw(plot)


 frame(wks)

end

