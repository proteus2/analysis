load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "uj"
 varname = "u"

 y1 = 1953
 y2 = 2006
 lat0 = 15

 iz = 0

 z = (/18.6,23.3/)  ; (/18.6,23.3/)  ; (/20,30,50,70/) ~ (/26.9,23.3,21.0,18.6/)
 p = (/70,30/)

 wks = cr_wks("fig/ws-y"+toint(lat0)+"-z"+toint(z(iz))+"-"+expname,600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt0 = (y2-y1+1)*12
 t0 = y1 + (fspan(1,nt0,nt0)-0.5)/12.

 fname  = "$DATD/L60CGW/tem/uan"+expname+".tem2_yzt.nc"
 fname2 = "/data18/kyh/dat/L60CGW/xm.uan"+expname+"a.pj.nc"
 f1 = addfile(fname ,"r")
 f2 = addfile(fname2,"r")

 nt = dimsizes( f1->time )
 t = t0(0:nt-1)

 nz = dimsizes(z)
 print( f1->z({z*1.e3})/1.e3 )

 w0   = f1->w_res(:,{z*1.e3},{-lat0:lat0})
 rho0 = f2->rho  (:,{z*1.e3},{-lat0:lat0})
 w = transpose( dim_avg(w0)*1.e3 )
; w = transpose( dim_avg(w0*1.e3*rho0) ) / \
;     conform_dims((/nz,nt/), dim_avg_n(dim_avg(rho0),0), 0)

;--------------------------------------- process 0
 wss = new((/nz,14/),float)
 do n=1, 12
   wss(:,n) = dim_avg(w(:,ispan(n-1,nt-1,12)))
 end do
 wss(:,0 ) = wss(:,12)
 wss(:,13) = wss(:,1 )

 wdes = new(dimsizes(w),float)
 do n=0, nt-1
   wdes(:,n) = w(:,n) - wss(:,mod(n+1,12))
 end do

 freq = fspan(0,nt/2,nt/2+1)/tofloat(y2-y1+1)
 psdw = new((/nz,nt/2+1/),float)
 wfft = ezfftf(wdes)
 psdw(:,0) = 0.
 psdw(:,1:nt/2) = 0.5*( wfft(0,:,:)^2 + wfft(1,:,:)^2 )/freq(1)
 psdw(:,nt/2) = psdw(:,nt/2)*2.

;--------------------------------------- cn level
 lcol = (/2,1,3,4/)

;--------------------------------------- contour
 nl = 3

 wi = (/126,40,60/)
 he = 25
 xf = (/20,20,86,20/)
 yf = 10. + (he+20.)*(/0.,1.,1./)

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = (/cexp,"",""/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "W:S:*:N: (mm s:S:-1:N:)"

 plot1 = xylineplot(wks,t,wdes(iz,:),"D",0,"T",1.0,"C",1)
 xyaddline(plot1,t,0,"D",0,"T",0.6,"C",colind)
 set_size(plot1,xf(0),yf(0),wi(0),he)
 axis_range(plot1,y1,y2+1,-999,999,False,False)
 axis_str(plot1,"Deseasonalized anomaly","year","W:S:*:N: (mm s:S:-1:N:)")
 draw(plot1)

 plot2 = xylineplot(wks,ispan(0,13,1),wss(iz,:),"D",0,"T",1.0,"C",1)
 set_size(plot2,xf(1),yf(1),wi(1),he)
 axis_range(plot2,0.5,12.5,-999,999,False,False)
 axis_str(plot2,"Seasonal cycle","month","W:S:*:N: (mm s:S:-1:N:)")
 setvalues plot2
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot2)

 plot3 = xylineplot(wks,freq,psdw(iz,:),"D",0,"T",1.0,"C",1)
 set_size(plot3,xf(2),yf(2),wi(2),he)
 axis_range(plot3,0,2.5,-999,999,False,False)
 axis_str(plot3,"PSD of deseasonalized anomaly","frequency (cyc yr:S:-1:N:)","PSD (mm:S:2:N: s:S:-2:N: / cyc yr:S:-1:N:)")
 setvalues plot3
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 0.5
  "tmXBMinorPerMajor" : 4
 end setvalues
 draw(plot3)

 frame(wks)

end



