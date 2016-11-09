load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 jobname = (/"tk","tc"/)
 expname = (/"AOL60CTL","AOL60CGW"/)
 ppid    = "j"
 varname = "T_1.5m"

 lat1 = 30
 lat2 = 70
 y0 = 1951
 y1 = 1951
 y2 = 2009
 mms = "djf"

 wks = cr_wks("fig/t_t15-"+mms+"_"+lat1+"-"+lat2+"-"+jobname(1),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(jobname)
 nyr = y2 - y1 + 1

 mm3s = (/"djf","mam","jja","son"/)
 mm3i = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
 do is=0, dimsizes(mm3s)-1
   if ( mms .eq. mm3s(is) ) then
     mm = mm3i(is,:)
   end if
 end do
 nm = dimsizes(mm)

do ie=0, ne-1

 fname = "../dat/"+expname(ie)+"/xm.uan"+jobname(ie)+"a.p"+ppid+".nc"
 f0 = addfile(fname,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp1 = f0->$varname$(n0:n0+nyr*12-1,0,{lat1:lat2})
 tmp2 = f0->p(n0:n0+nyr*12-1,0,{lat1:lat2})

 tmp = dim_avg( tmp1 * (tmp2/1.e5)^0.286 )

 if (ie .eq. 0) then
   vara = new((/ne,nm,nyr/),float)
 end if
 do im=0, nm-1
   vara(ie,im,:) = tmp(ispan(mm(im)-1,nyr*12-1,12))
 end do

 delete(tmp)

end do  ; ie

 t = ispan(y1,y2,1)

;--------------------------------------- process 0
 var0 = dim_avg_n(vara,1) - 273.15

 var1 = runave(var0,9,0)

 var = array_append_record( var0, var1, 0 )

;--------------------------------------- cn level
 lcol = (/1,2,1,2/)
 ltck = (/0.5,0.5,2,2/)
 ldsh = (/2,2,0,0/)

;--------------------------------------- contour
 wi = 64
 he = 24
 xf = 16 + (/0,1,2/)*40
 yf = 10

 ystr  = new(ne+1,string)

 mstr = "" ;"T at 1.5 m"
 xstr = "time"
 ystr = "T (:S:o:N:C)"

 plot = xylineplot(wks,t,var,"D",ldsh,"T",ltck,"C",lcol)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,y1,y2,-999,999,False,False)
 axis_str(plot,mstr,xstr,ystr)

 draw(plot)


 frame(wks)

end

