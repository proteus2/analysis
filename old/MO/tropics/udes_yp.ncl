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

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 m2 = 0*12
 m3 = 0*12
 nmon0 = nmon-4

 wks = cr_wks("udes_yp-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = f1->u(4:nmon-1,:,{-lat0:lat0},0)
   lat1 = f1->latitude({-lat0:lat0})
   p1 = f1->p
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = f2->u(4:nmon-1,:,{-lat0:lat0},0)

   ny = dimsizes(lat1)
   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = f3->u(m1:m1+nmon0-1,:,{-lat0:lat0})
   u4 = f3->u(m2:m2+nmon0-1,:,{-lat0:lat0})
   u5 = f3->u(m3:m3+nmon0-1,:,{-lat0:lat0})
   lat = f3->latitude({-lat0:lat0})
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

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

   urun  = runave(udes ,9,0)
   urun@_FillValue = 1.e32
   urun0 = runave(udes0,9,0)
   urun0@_FillValue = 1.e32

   uano  = udes  - urun
   uano0 = udes0 - urun0

   uu = (/udes,urun,uano/)
   uu@_FillValue = 1.e32
   uu0 = (/udes0,urun0,uano0/)
   uu0@_FillValue = 1.e32

   usd  = dim_stddev(uu)
   usd0 = dim_stddev(uu0)

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(2,24,12)
 cnfil = (/0,10,9,8,7,6,5,4,3,2,26,25,24,23,22,21,19/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38 + (/0,0,1/)*8
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = (/"","","lat (deg)"/)
 ystr = (/"log-p height (km)","","log-p height (km)"/)
 lbopt = new((/3,3/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(m s:S:-1:N:)","",""/)


do isc=0, 2
do iv=0, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,lat1,zp1,usd(isc,iv,:,:),cnlev,cnfil,lbopt(isc,iv))
 else
   contour1 = cnshadeplot(wks,lat,zp,dim_avg_n(usd0(isc,:,:,:),0),cnlev,cnfil,lbopt(isc,iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.0
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(isc),wi,he)
 axis_range(contour1,-lat0,lat0,15,65,False,False)
 axis_str(contour1,mstr(isc,iv),xstr(isc),ystr(iv))

 cnfinalize(contour1,wks)


 draw(contour1)

end do
end do


 frame(wks)

end



