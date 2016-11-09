load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 expn = (/"f","0"/)

 month = (/1,7/)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 nyear = 32
 y1    = 0

 mstr = (/"CGW","CGW - CTL","CTL - MERRA","CGW - MERRA"/)

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/dat"
; file_dat1 = dir_dat+"/era-int.u_ypt."+year1+"-"+year2+".nc"
 file_dat1 = dir_dat+"/merra.u_ypt."+year1+"-"+year2+".nc"
 file_dat3 = dir_dat+"/uh-p_"+expn+".nc"

 f1 = addfile(file_dat1,"r")
 f3 = addfile(file_dat3(0),"r")
 f4 = addfile(file_dat3(1),"r")

 lat1 = f1->latitude
 lat3 = f3->latitude
; p1  = f1->levelist
 p1  = f1->levels
 p3  = f3->p

 zp1 = 7.0*log(1.e3/p1)
 zp3 = 7.0*log(1.e3/p3)

 nmon = nyear*12
 m_1 = y1*12
 m_2 = m_1+nmon-1

 nmon_c = dimsizes(f3->t(4:))

; u1i = doubletofloat( f1->u(m_1:m_2,:,:) )
 u1i = f1->u(m_1:m_2,:,:)
 u3i = f3->u(4:,:,:,0)
 u4i = f4->u(4:,:,:,0)
 u1i!1 = "p"
 u1i&p = p1

 ny = dimsizes(lat1)
 nz = dimsizes(p1)
 ny_c = dimsizes(lat3)
 nz_c = dimsizes(p3)

 nyr = new(12,integer)
 um0 = new((/12,nz,ny/),float,1.e32)
 um1 = new((/12,nz_c,ny_c/),float,1.e32)
 um2 = new((/12,nz_c,ny_c/),float,1.e32)
 do n=0, 11
   um0(n,:,:) = dim_avg_n( u1i(ispan(n,nmon-1,12),:,:),0 )
   nyr(n) = dimsizes(ispan(n,nmon_c-1,12))
   if (nyr(n) .gt. 1) then
     um1(n,:,:) = dim_avg_n( u3i(ispan(n,nmon_c-1,12),:,:),0 )
     um2(n,:,:) = dim_avg_n( u4i(ispan(n,nmon_c-1,12),:,:),0 )
   else
     um1(n,:,:) = u3i(n,:,:)
     um2(n,:,:) = u4i(n,:,:)
   end if
 end do

;--------------------------------------- process 0
 do im=0, dimsizes(month)-1

 wks = cr_wks("u_yp-"+month(im)+"-"+expn(0),600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 print(nyr(month(im)-1)+"")

 uf0 = um1(month(im)-1,:,:)
 uf1 = uf0 - um2(month(im)-1,:,:)

 um0i = new((/nz_c,ny_c/),float,1.e32)
 do j=0, ny_c-1
   um0i(:,j) = ftcurv(zp1,um0(month(im)-1,:,j),zp3)
 end do
 uf2 = um2(month(im)-1,:,:) - um0i
 uf3 = uf0 - um0i

;--------------------------------------- line attribute
 cnlev = (/-80,-40,-20,-10,-5,0,5,10,20,40,80/)
 cnfil = (/4,5,7,8,9,0,0,14,15,17,18,21/)

 cnlev2 = (/-30,-20,-10,10,20,30/)
 cnfil2 = (/6,8,9,0,14,15,17/)

;--------------------------------------- plot
 wi = 58
 he = 32
 xf = 20 + (/0,1/)*72
 yf = 10 + (/0,1/)*56

 contour1 = cnshadeplot(wks,lat3,zp3,uf0,cnlev,cnfil,1)
 cnaddline(contour1,wks,lat3,zp3,uf0,fspan(-100,100,41),0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",1)
; cnline_att(contour1,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,-90,90,10,64,False,False)
 axis_str(contour1,mstr(0),"lat","log-p height (km)")

 cnfinalize(contour1,wks)

 draw(contour1)


 contour2 = cnshadeplot(wks,lat3,zp3,uf1,cnlev2,cnfil2,1)
 cnaddline(contour2,wks,lat3,zp3,uf1,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour2,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour2,xf(1),yf(0),wi,he)
 axis_range(contour2,-90,90,10,64,False,False)
 axis_str(contour2,mstr(1),"lat","log-p height (km)")

 cnfinalize(contour2,wks)

 draw(contour2)


 contour3 = cnshadeplot(wks,lat3,zp3,uf2,cnlev2,cnfil2,1)
 cnaddline(contour3,wks,lat3,zp3,uf2,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour3,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour3,xf(0),yf(1),wi,he)
 axis_range(contour3,-90,90,10,64,False,False)
 axis_str(contour3,mstr(2),"lat","log-p height (km)")

 cnfinalize(contour3,wks)

 draw(contour3)


 contour4 = cnshadeplot(wks,lat3,zp3,uf3,cnlev2,cnfil2,1)
 cnaddline(contour4,wks,lat3,zp3,uf3,fspan(-100,100,41),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour4,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour4,xf(1),yf(1),wi,he)
 axis_range(contour4,-90,90,10,64,False,False)
 axis_str(contour4,mstr(3),"lat","log-p height (km)")

 cnfinalize(contour4,wks)

 draw(contour4)


 frame(wks)


 end do

end

