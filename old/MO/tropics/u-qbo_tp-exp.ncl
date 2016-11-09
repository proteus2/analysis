load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 expn = "fd"

 wks = cr_wks("u-qbo_tp-"+expn,600,0)
 gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 lat_r = 5;10

 nyear = 8;12
 y1    = 0

 mstr1 = "U (UM-CTL, monthly, "+lat_r+"N-"+lat_r+"S, zonal mean)"

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/tropics/dat"
 file_dat1 = dir_dat+"/uh-p_"+expn+".nc"

 f1 = addfile(file_dat1,"r")

; lat = f1->lat({-lat_r:lat_r})
 p1  = f1->p

 nmon = nyear*12
 m_1 = 4 + y1*12
 m_2 = m_1+nmon-1

 u1i = dim_avg( f1->u(m_1:m_2,:,{-lat_r:lat_r},0) )

 np1 = dimsizes(p1)
 nmon = dimsizes(f1->t(m_1:m_2))

 zp1 = 7.0*log(1.e3/p1)

 mon = fspan(0,nmon-1,nmon)
 u1 = new((/np1,nmon/),float)
 do i=0, nmon-1
   u1(:,i) = u1i(i,:)
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

 contour1 = cnshadeplot(wks,mon,zp1,u1,cnlev,cnfil,1)
 cnaddline(contour1,wks,mon,zp1,u1,(/0/),0, \
           "D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)
; cnline_att(contour1,"D",(/0,0,1/),"T",(/1.0,0.5,1.0/),"C",1)

 set_size(contour1,xf,yf(0),wi,he)
; axis_range(contour1,-999,999,15,47,False,False)
 axis_range(contour1,0,143,15,47,False,False)
 axis_str(contour1,mstr1,"month","log-p height (km)")

 cnfinalize(contour1,wks)


 draw(contour1)


 frame(wks)

end

