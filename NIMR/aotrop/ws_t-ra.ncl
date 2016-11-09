load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "merra"

 y0 = 1979
 y1 = 1979
 y2 = 2015
 lat0 = 15

; z = (/18.6,23.3,26.9,31.0/)  ; (/20,30,50,70/) ~ (/26.9,23.3,21.0,18.6/)
 p = (/70,30,20,10/)

 wks = cr_wks("fig/ws_t-y"+toint(lat0)+"-"+expname,600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt0 = (y2-y1+1)*12
 n0 = (y1-y0)*12
 t0 = y1 + (fspan(1,nt0,nt0)-0.5)/12.

 fname  = "/data18/kyh/dat/qbo/"+expname+".tem_ypt.1979-2006.nc"
 f1 = addfile(fname ,"r")

 nt = dimsizes( f1->time ) - n0
 t = t0(0:nt-1)

 nz = dimsizes(p)
 print( f1->p({p}) )

 w0   = f1->w_res(n0:,{p},{-lat0:lat0})
 w = transpose( dim_avg(w0)*1.e3 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 lcol = (/2,4,3,5/)

;--------------------------------------- contour
 nl = 2
 ys = (/1979,1997/)
 yl = 18
 yls = 2

 wi = 126
 he = 20
 xf = 34
 yf = 10 + (/0,1,2/)*28

 wi2 = tofloat(wi)/yl*yls

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = (/expname,""/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "W:S:*:N: (mm s:S:-1:N:)"

do il=0, nl-1

 plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
 do k=0, nz-1
   xyaddline(plot,t,w(k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,ys(il),ys(il)+yl,0,0.6,False,False)
 axis_str(plot,"",xstr(il),ystr(il))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

if (il .gt. 0) then
  setvalues plot
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  plot2 = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
  do k=0, nz-1
    xyaddline(plot2,t,w(k,:),"D",0,"T",1.0,"C",lcol(k))
  end do
  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),0,0.6,False,False)
  axis_str(plot2,"","",ystr(il))
  setvalues plot2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  draw(plot2)
end if

 draw(plot)

end do


 frame(wks)

end



