load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.0

load "read_dims.ncl"

begin

;--------------------------------------- parameters
 models = (/"CMAM","MIROC-ESM","MRI"/)
; models = (/"CMAM","MRI"/)
 y1 = 1993
 y2 = 2007
 mmi = "05"
 dd = "01"
 hh = "06"
 lata = 5.

 nm = dimsizes(models)

do im=0, nm-1

 tag = 0
 do yr=y1, y2

   sysout = systemfunc("env_qboi.sh 5 ua 6hr XYPT r1i1p1 "+yr+mmi)

   i0 = ind( sysout .eq. "start_"+models(im) ) + 2
   nfile = stringtoint(sysout(i0-1))
   files = sysout(i0:i0+nfile-1)

   f0 = addfile(files(0),"r")

   ; time: n1
   t = f0->$f0->ua!0$(0:4)
   date = cd_calendar(t,-3)
   n1 = ind( date+"" .eq. yr+mmi+dd+hh )
   if ( ismissing(n1) ) then
     print(" MISSING: "+models(im)+" for "+(yr+mmi+dd+hh))
     exit
   end if

   if (tag .eq. 0) then
     tag = 1

     lat = getlat_q(models(im),f0,"ua")
     lat_bnd1 = lat@bnd1
     lat_bnd2 = lat@bnd2
     dlat = lat_bnd2 - lat_bnd1

     lev = getlev_q(models(im),f0,"ua")
     pscale = 1.
     levp = tofloat( lev )
     if ( isatt(lev, "units") ) then
     if (lev@units .eq. "Pa") then
       pscale = 100.
     end if
     end if
     levp = levp/pscale

     if (lat(0) .lt. lat(1)) then
       j1 = min( ind( lat .ge. -lata ) )
       j2 = max( ind( lat .le. lata  ) )
     else
       j1 = min( ind( lat .le. lata  ) )
       j2 = max( ind( lat .ge. -lata ) )
     end if
 
     if ( abs(lat(j1)) .ne. lata .or.  \
          abs(lat(j2)) .ne. lata ) then
       tag = 2
     end if
   end if

   u0 = dim_avg_Wrap( f0->ua(n1,:,:,:) )

   u = dim_sum_Wrap( u0(:,j1:j2)*conform_dims(dimsizes(u0(:,j1:j2)),  \
                                 dlat(j1:j2),1) )
   if (tag .eq. 2) then
print("U needs to be updated")
   end if

 end do

print(lat_bnd1+"   "+lat_bnd2)

 delete( files )
 delete( [/ u, u0 /] )
 delete( t )
 delete( [/ lat , lev, levp /] )
 delete( [/ lat_bnd1, lat_bnd2, dlat /] )

end do

exit

 latb = 5.


;#ncks -h -d t,$n1,$n2 $f0 -O zz.nc
;#ncwa -h -a latitude -d t,$n1,$n2 -d latitude,-5.0,5.0 $f0 -O zz.nc
;ncwa -h -a latitude -w $wgt  -d t,$n1,$n2 -d latitude,68,76 $f0 -O zz2.nc

 wks = cr_wks("fig/u_init",600,0)
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

