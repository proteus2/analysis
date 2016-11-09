load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 30.
 kz  = (/9,15,21/) ;9,12,15,18,21
 kz0 = (/7,13,19/)  ;7,9,13,15,19

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 m2 = 0*12
 m3 = 0*12
 nmon0 = nmon-4

 wks = cr_wks("udes_y-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = f1->u(4:nmon-1,kz,{-lat0:lat0},0)
   lat1 = f1->latitude({-lat0:lat0})
   p1 = f1->p(kz)
print(p1)

   tmp = 10.^(floattoint(log10(p1)+100.)-100)
   p1 = round(p1/tmp,0)*tmp

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = f2->u(4:nmon-1,kz,{-lat0:lat0},0)

   ny = dimsizes(lat1)
   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = f3->u(m1:m1+nmon0-1,kz0,{-lat0:lat0})
   u4 = f3->u(m2:m2+nmon0-1,kz0,{-lat0:lat0})
   u5 = f3->u(m3:m3+nmon0-1,kz0,{-lat0:lat0})
   lat = f3->latitude({-lat0:lat0})
   p  = f3->levels(kz0)
print(p)

   ny3 = dimsizes(lat)
   nz3 = dimsizes(p)
   nt3 = nmon0

   u = new((/2,nz,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/u2(n,:,:),u1(n,:,:)/)
   end do

   ur = new((/3,nz3,ny3,nt3/),float)
   do n=0, nt3-1
     ur(0,:,:,n) = u3(n,:,:)
     ur(1,:,:,n) = u4(n,:,:)
     ur(2,:,:,n) = u5(n,:,:)
   end do

   t = fspan(1,nt,nt)
   t3 = fspan(1,nt3,nt3)

   uss = new((/2,nz,ny,12/),float)
   do n=0, 11
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n,nt-1,12)))
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,:,n) = u(:,:,:,n) - uss(:,:,:,mod(n,12))
   end do

   uss0 = new((/3,nz3,ny,12/),float)
   do n=0, 11
     uss0(:,:,:,n) = dim_avg(ur(:,:,:,ispan(n,nt3-1,12)))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,:,:,n) = ur(:,:,:,n) - uss0(:,:,:,mod(n,12))
   end do

   ifilt = nt/24-1
   ifilt0 = nt3/24-1

   fftu = ezfftf(udes)
   fftu(:,:,:,:,ifilt+1:) = 0.
   urun = ezfftb(fftu,0.)

   fftu0 = ezfftf(udes0)
   fftu0(:,:,:,:,ifilt0+1:) = 0.
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
 llbl = (/"CTL","SC","MERRA"/)
 lthc = (/1.,1.,2./)
 ldsh = (/14,0,0/)

;--------------------------------------- contour
 wi = 35
 he = 24
 xf = 12 + (/0,1,2/)*43
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/nz,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Total","Low Freq.","High Freq."/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "lat (deg)"
 ystr = (/"Std. Dev. (m s:S:-1:N:)","",""/)
 txtstr = new((/nz,3/),string)
 txtstr(:,:) = ""
 txtstr(:,0) = p1(:)+" hPa"
 ymax = (/20,12/)


do k=0, nz-1
do isc=0, 1;2

 plot = xylineplot(wks,lat1,usd(isc+1,:,k,:),"D",ldsh(0:1),"T",lthc(0:1),"C",1)
 xyaddline(plot,lat,dim_avg_n(usd0(isc+1,:,k,:),0),"D",ldsh(2),"T",lthc(2),"C",1)

 set_size(plot,xf(isc),yf(k),wi,he)
 axis_range(plot,-lat0,lat0,-999,999,False,False)
 axis_str(plot,mstr(k,isc+1),xstr(k),ystr(isc))
; getvalues plot
;  "trYMaxF" : ymax
; end getvalues
 setvalues plot
  "trYMaxF" : ymax(isc) ;ymax*1.4
 end setvalues

 draw(plot)

 dr_txt(wks,txtstr(k,isc),xf(isc)+6,yf(k)+3,22,6,0.)

end do
end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",lthc,"C",1,2)
 set_size_lg(lg,xf(1),yf(0),wi*0.15,he*0.4)
 draw(lg)


 frame(wks)

end



