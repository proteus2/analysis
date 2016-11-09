load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expn = (/"fd","y"/)
 opt_mean = 0   ;  0 - raw  ;  1 - anomaly  ;  2 - deseasonalized anomaly

 nmon = 4+144
 lat0 = 5

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 m2 = 0*12
 m3 = 0*12
 nmon0 = 32*12;nmon-4

 wks = cr_wks("qbo-diag",600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 nyear = 32
 y1    = 0

;--------------------------------------- read data
 dir_dat   = "../dat"
 file_dat1 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat2 = dir_dat+"/uh-p_"+expn(1)+".nc"
 file_dat3 = dir_dat+"/uh-p_"+expn(0)+".nc"
 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f3 = addfile(file_dat3,"r")

 lat1 = f1->latitude({-lat0:lat0})
 lat2 = f2->latitude({-lat0:lat0})

 p1  = f1->levels
 p2  = f2->p
 zp1 = 7.0*log(1.e3/p1)
 zp2 = 7.0*log(1.e3/p2)

 ny1 = dimsizes(lat1)
 ny2 = dimsizes(lat2)
 nz1 = dimsizes(p1)
 nz2 = dimsizes(p2)
 nt1 = nmon0
 nt2 = nmon-4

 u1ia = dim_avg( f1->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
 u1ib = dim_avg( f1->u(m2:m2+nmon0-1,:,{-lat0:lat0}) )
 u1ic = dim_avg( f1->u(m3:m3+nmon0-1,:,{-lat0:lat0}) )
 u2i = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )
 u3i = dim_avg( f3->u(4:nmon-1,:,{-lat0:lat0},0) )

 ur = new((/3,nz1,nt1/),float)
 do n=0, nt1-1
   ur(:,:,n) = (/ u1ia(n,:), u1ib(n,:), u1ic(n,:) /)
 end do

 u = new((/2,nz2,nt2/),float)
 do n=0, nt2-1
   u(:,:,n) = (/ u2i(n,:), u3i(n,:) /)
 end do

 uss1 = new((/3,nz1,12/),float)
 do n=0, 11
   uss1(:,:,n) = dim_avg(ur(:,:,ispan(n,nt1-1,12)))
 end do
 udes1 = new(dimsizes(ur),float)
 do n=0, nt1-1
   udes1(:,:,n) = ur(:,:,n) ;- uss1(:,:,mod(n,12))
 end do

 uss2 = new((/2,nz2,12/),float)
 do n=0, 11
   uss2(:,:,n) = dim_avg(u(:,:,ispan(n,nt2-1,12)))
 end do
 udes2 = new(dimsizes(u),float)
 do n=0, nt2-1
   udes2(:,:,n) = u(:,:,n) ;- uss2(:,:,mod(n,12))
 end do

 ifilt1 = nt1/12-1
 ifilt2 = nt2/12-1

 fftu1 = ezfftf(udes1)
 fftu1(:,:,:,ifilt1+1:) = 0.
 bar1 = fftu1@xbar
 if (opt_mean .ne. 0) then
   bar1 = bar1*0.
 end if
 ulow1 = ezfftb(fftu1,bar1)

 fftu2 = ezfftf(udes2)
 fftu2(:,:,:,ifilt2+1:) = 0.
 bar2 = fftu2@xbar
 if (opt_mean .ne. 0) then
   bar2 = bar2*0.
 end if
 ulow2 = ezfftb(fftu2,bar2)

 u1 = new((/4,3,nz1/),float,1.e32)
 u2 = new((/4,2,nz2/),float,1.e32)
 u1 = 0.
 u2 = 0.

 u1pn = (/ulow1,ulow1/)
 u2pn = (/ulow2,ulow2/)

 do ie=0, 2
 do k=0, nz1-1
   do n=0, nt1-1
     if (u1pn(0,ie,k,n) .gt. 0.) then
       u1(0,ie,k) = u1(0,ie,k) + 1.
       u1pn(1,ie,k,n) = 0.
     else
       u1pn(0,ie,k,n) = 0.
     end if
   end do
 end do
 end do
 do ie=0, 1
 do k=0, nz2-1
   do n=0, nt2-1
     if (u2pn(0,ie,k,n) .gt. 0.) then
       u2(0,ie,k) = u2(0,ie,k) + 1.
       u2pn(1,ie,k,n) = 0.
     else
       u2pn(0,ie,k,n) = 0.
     end if
   end do
 end do
 end do
 u1(1,:,:) = nt1 - u1(0,:,:)
 u2(1,:,:) = nt2 - u2(0,:,:)

 u1(2:3,:,:) = sqrt( dim_sum( u1pn(0:1,:,:,:)^2 )/(u1(0:1,:,:)+0.0000001) )
 u2(2:3,:,:) = sqrt( dim_sum( u2pn(0:1,:,:,:)^2 )/(u2(0:1,:,:)+0.0000001) )

 u1(0:1,:,:) = 0.
 u2(0:1,:,:) = 0.
 stmp = (/"+","-"/)
 do j=0, 1
 do ie=0, 2
 do k=0, nz1-1
   tmp1 = 0
   tmp2 = 0
   cnt = -999
   do n=0, nt1-1
     if (u1pn(j,ie,k,n) .eq. 0.) then
       if (cnt .gt. 6) then
         tmp1 = tmp1 + cnt
         tmp2 = tmp2 + 1
         if (p1(k) .eq. 50.) then
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
     u1(j,ie,k) = int2flt(tmp1)/tmp2
   else
     u1(j,ie,k) = 1.e32
   end if
 end do
 end do
 end do
 do j=0, 1
 do ie=0, 1
 do k=0, nz2-1
   tmp1 = 0
   tmp2 = 0
   cnt = -999
   do n=0, nt2-1
     if (u2pn(j,ie,k,n) .eq. 0.) then
       if (cnt .gt. 6) then
         tmp1 = tmp1 + cnt
         tmp2 = tmp2 + 1
         if (abs(p2(k)-50.)/50. .le. 0.1) then
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
     u2(j,ie,k) = int2flt(tmp1)/tmp2
   else
     u2(j,ie,k) = 1.e32
   end if
 end do
 end do
 end do

;--------------------------------------- process 0

;--------------------------------------- line attribute
 llbl = (/"CTL","SC","MERRA"/)
 lthc = (/1.,1.,2./)
 ldsh = (/14,0,0/)

;--------------------------------------- plot
 wi = 35
 he = 36
 xf = 12 + (/0,1,2,3/)*43
 yf = 10 + (/0,1,2,3,4/)*42

 mstr = (/"Total","Low Freq.","High Freq.",""/)
 xstr = "Std. Dev. (m s:S:-1:N:)"
 ystr = (/"log-p height (km)","","",""/)
 xmin = (/9,9,0,0/)
 xmax = (/20,20,40,40/)

do iv=0, 3

 plot = xylineplot(wks,u2(iv,:,:),zp2,"D",ldsh(0:1),"T",lthc(0:1),"C",1)
 xyaddline(plot,dim_avg_n(u1(iv,:,:),0),zp1,"D",ldsh(2),"T",lthc(2),"C",1)

 set_size(plot,xf(iv),yf,wi,he)
 axis_range(plot,xmin(iv),xmax(iv),15,40,False,False)
 axis_str(plot,mstr(iv),xstr,ystr(iv))

 draw(plot)

end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",lthc,"C",1,2)
 set_size_lg(lg,xf+wi*0.55,yf(0),wi*0.25,he*0.3)
 draw(lg)


 frame(wks)

end

