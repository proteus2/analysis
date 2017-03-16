load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fname0 = "dat/O3_Yonsei_SPARC_1952_2009q_L60_yzt.nc"
 y0 = 1940
 m0 = 11

 y1 = 1979
; y2 = 2006
 y2 = 2008
;y1=1994
;y2=2005
 latb = 5.


 wks = cr_wks("fig/o3prof",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 n1 = (y1 - y0)*12
 if (m0 .ne. 1) then
   n1 = n1 - (m0 - 1)
 end if
 n2 = n1 + (y2 - y1 + 1)*12 - 1
 
 ; to include half-month edges
 n1 = n1 - 1
 n2 = n2 + 1

 f0 = addfile(fname0,"r")

; m2v = 28.97/48.00
 m2v = 1./(48.*0.62198/18.)  ; following Andrew (based on the UM values)
 oz = f0->O3(n1:n2,:,{-latb:latb}) * 1.e6 * m2v
 z = f0->hybrid_ht * 1.e-3
 lat0 = f0->latitude({latb})

 ndim = dimsizes(oz)
 ny = ndim(2)
 nt = ndim(0)

 wgt_t = new(nt,float)
 wgt_t = 1.
 wgt_t((/0,nt-1/)) = 0.5

 wgt_y = new(ny,float)
 wgt_y = 1.
; if (lat0 .eq. latb) then
;   wgt_y((/0,ny-1/)) = 0.5
; else
;   print("Rewrite the script considering weighting.")
; end if

 ozc = dim_avg_wgt_n(dim_avg_wgt(oz,wgt_y,0),wgt_t,0,0)

;print(z+"   "+ozc)

 pm = asciiread("dat/pm_t",-1,"float")

print(z+"   "+ozc+"   "+pm+"   "+(7.0*log(1.e3/pm)))

;--------------------------------------- cn level

;--------------------------------------- ctr
 wi = 40
 he = 60
 xf = 20
 yf = 10

 xstr = "concentration [ppm (vol)]"
 ystr = "z [km]"

 plot1 = xylineplot(wks,ozc,z,"D",0,"T",1,"C",1)
 xyaddmarker(plot1,1,"S",3.0,"T",0.6,"C",1)

 set_size(plot1,xf,yf,wi,he)
 axis_range(plot1,0,12,0,80,False,False)
 axis_str(plot1,"Annual Equatorial (5N-5S) Ozone",xstr,ystr)

 draw(plot1)

 frame(wks)

end

