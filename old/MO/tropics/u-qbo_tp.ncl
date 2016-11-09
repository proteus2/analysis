load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 wks = cr_wks("u-qbo_tp",600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 lat_r = 10

 nyear = 12
 y1    = 0

 mstr1 = "U (ERA-I, monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"
 mstr2 = "U (MERRA, monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"

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
 nmon = dimsizes(f1->time(m_1:m_2))

 zp1 = 7.0*log(1.e3/p1)
 zp2 = 7.0*log(1.e3/p2)

 mon = fspan(0,nmon-1,nmon)
 u1 = new((/np1,nmon/),float)
 u2 = new((/np2,nmon/),float)
 do i=0, nmon-1
   u1(:,i) = u1i(i,:)
   u2(:,i) = u2i(i,:)
 end do

;--------------------------------------- process 0

;--------------------------------------- line attribute
 cnlev = (/-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60/)*0.5
 cnfil = (/colind,3,4,5,6,7,8,9,0,0,14,15,16,17,18,19,20/)

;--------------------------------------- plot
 wi = 132
 he = 32
 xf = 20
 yf = 10 + (/0,1/)*50

 contour1 = cnshadeplot(wks,mon,zp1,u1,cnlev,cnfil,0)
 cnaddline(contour1,wks,mon,zp1,u1,(/0/),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)
; cnline_att(contour1,"D",(/0,0,1/),"T",(/1.0,0.5,1.0/),"C",1)

 set_size(contour1,xf,yf(0),wi,he)
 axis_range(contour1,-999,999,15,47,False,False)
 axis_str(contour1,mstr1,"month","log-p height (km)")

 cnfinalize(contour1,wks)

 contour2 = cnshadeplot(wks,mon,zp2,u2,cnlev,cnfil,1)
 cnaddline(contour2,wks,mon,zp2,u2,(/0/),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)
; cnline_att(contour2,"D",(/0,0,1/),"T",(/1.0,0.5,1.0/),"C",1)

 set_size(contour2,xf,yf(1),wi,he)
 axis_range(contour2,-999,999,15,47,False,False)
 axis_str(contour2,mstr2,"month","log-p height (km)")

 cnfinalize(contour2,wks)


 draw(contour1)
 draw(contour2)


 frame(wks)

end

