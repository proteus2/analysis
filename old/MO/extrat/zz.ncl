load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 expn = (/"0","fa","fd"/)

 month = (/4,5/)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 nyear = 30
 y1    = 1

 mstr = (/"CTL","CTL - MERRA", \
          "CGW2_1 - CTL","CGW2_1 - MERRA", \
          "CGW2_2 - CTL","CGW2_2 - MERRA"/)

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

 lat1 = f1->latitude
 lat2 = f2->latitude
 p1  = f1->levels
 p2  = f2->p

 zp1 = 7.0*log(1.e3/p1)
 zp2 = 7.0*log(1.e3/p2)

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 
 nmon_c = min( (/ dimsizes(f3->t(4:)), dimsizes(f4->t(4:)) /) ) - 1

 u1i = (f1->u(m_1:m_2,:,:) + f1->u(m_1+1:m_2+1,:,:) )/2.
 u2i = (f2->u(4:nmon_c+3,:,:,0) + f2->u(5:nmon_c+4,:,:,0))/2.
 u3i = (f3->u(4:nmon_c+3,:,:,0) + f3->u(5:nmon_c+4,:,:,0))/2.
 u4i = (f4->u(4:nmon_c+3,:,:,0) + f4->u(5:nmon_c+4,:,:,0))/2.
 u1i!1 = "p"
 u1i&p = p1

 ny = dimsizes(lat1)
 nz = dimsizes(p1)
 ny_c = dimsizes(lat2)
 nz_c = dimsizes(p2)

 nyr = new(12,integer)
 um0 = new((/12,nz,ny/),float,1.e32)
 um1 = new((/12,nz_c,ny_c/),float,1.e32)
 um2 = new((/12,nz_c,ny_c/),float,1.e32)
 um3 = new((/12,nz_c,ny_c/),float,1.e32)
 do n=0, 11
   um0(n,:,:) = dim_avg_n( u1i(ispan(n,nmon-1,12),:,:),0 )
   nyr(n) = dimsizes(ispan(n,nmon_c-1,12))
   if (nyr(n) .gt. 1) then
     um1(n,:,:) = dim_avg_n( u2i(ispan(n,nmon_c-1,12),:,:),0 )
     um2(n,:,:) = dim_avg_n( u3i(ispan(n,nmon_c-1,12),:,:),0 )
     um3(n,:,:) = dim_avg_n( u4i(ispan(n,nmon_c-1,12),:,:),0 )
   else
     um1(n,:,:) = u2i(n,:,:)
     um2(n,:,:) = u3i(n,:,:)
     um3(n,:,:) = u4i(n,:,:)
   end if
 end do

;--------------------------------------- process 0
 do im=0, dimsizes(month)-1

 wks = cr_wks("fig/zz-"+month(im),600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 print(nyr(month(im)-1)+"")

 uf0 = um1(month(im)-1,:,:)
 uf2 = um2(month(im)-1,:,:) - uf0
 uf4 = um3(month(im)-1,:,:) - uf0

 um0i = new((/nz_c,ny_c/),float,1.e32)
 do j=0, ny_c-1
   um0i(:,j) = ftcurv(zp1,um0(month(im)-1,:,j),zp2)
 end do
 uf1 = um1(month(im)-1,:,:) - um0i
 uf3 = um2(month(im)-1,:,:) - um0i
 uf5 = um3(month(im)-1,:,:) - um0i

;--------------------------------------- line attribute
 cnlev = (/-80,-40,-20,-10,-5,0,5,10,20,40,80/)
 cnfil = (/4,5,7,8,9,0,0,14,15,17,18,21/)

 cnlev2 = (/-30,-20,-10,10,20,30/)*0.1
 cnfil2 = (/6,8,9,0,14,15,17/)

;--------------------------------------- plot
 wi = 58
 he = 32
 xf = 20 + (/0,1/)*72
 yf = 10 + (/0,1,2/)*56

 contour1 = cnshadeplot(wks,lat2,zp2,uf0,cnlev,cnfil,1)
 cnaddline(contour1,wks,lat2,zp2,uf0,fspan(-100,100,41),0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",1)
; cnline_att(contour1,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,-90,90,10,64,False,False)
 axis_str(contour1,mstr(0),"lat","log-p height (km)")

 cnfinalize(contour1,wks)

 draw(contour1)


 contour2 = cnshadeplot(wks,lat2,zp2,uf1,cnlev2,cnfil2,1)
 cnaddline(contour2,wks,lat2,zp2,uf1,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour2,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour2,xf(1),yf(0),wi,he)
 axis_range(contour2,-90,90,10,64,False,False)
 axis_str(contour2,mstr(1),"lat","log-p height (km)")

 cnfinalize(contour2,wks)

 draw(contour2)


 contour3 = cnshadeplot(wks,lat2,zp2,uf2,cnlev2,cnfil2,1)
 cnaddline(contour3,wks,lat2,zp2,uf2,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour3,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour3,xf(0),yf(1),wi,he)
 axis_range(contour3,-90,90,10,64,False,False)
 axis_str(contour3,mstr(2),"lat","log-p height (km)")

 cnfinalize(contour3,wks)

 draw(contour3)


 contour4 = cnshadeplot(wks,lat2,zp2,uf3,cnlev2,cnfil2,1)
 cnaddline(contour4,wks,lat2,zp2,uf3,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour4,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour4,xf(1),yf(1),wi,he)
 axis_range(contour4,-90,90,10,64,False,False)
 axis_str(contour4,mstr(3),"lat","log-p height (km)")

 cnfinalize(contour4,wks)

 draw(contour4)


 contour5 = cnshadeplot(wks,lat2,zp2,uf4,cnlev2,cnfil2,1)
 cnaddline(contour5,wks,lat2,zp2,uf4,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour5,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour5,xf(0),yf(2),wi,he)
 axis_range(contour5,-90,90,10,64,False,False)
 axis_str(contour5,mstr(4),"lat","log-p height (km)")

 cnfinalize(contour5,wks)

 draw(contour5)


 contour6 = cnshadeplot(wks,lat2,zp2,uf5,cnlev2,cnfil2,1)
 cnaddline(contour6,wks,lat2,zp2,uf5,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour6,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour6,xf(1),yf(2),wi,he)
 axis_range(contour6,-90,90,10,64,False,False)
 axis_str(contour6,mstr(5),"lat","log-p height (km)")

 cnfinalize(contour6,wks)

 draw(contour6)


 frame(wks)


 end do

end

