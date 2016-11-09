load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 jobname = (/"tk","tc"/)
; jobname = (/"ui","uj"/)
 expname = "L60"+(/"CTL","CGW"/)
; expname = "L60"+(/"CTL","CGW-t"/)
 expname = "AO"+expname
 year1 = 1951;+2
 year2 = 1980;+26

 wks = cr_wks("fig/easmi-"+jobname(1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- data reading
 nt = year2 - year1 + 1
 time = fspan(year1,year2,nt)

 easmi = new((/2,nt/),float)
 easmi = 0.

 n = 0
 do iy=year1, year2
   do im=6, 8
;     file_uv = "/data15/HG2AMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
     file_uv = "/data18/HG2CMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
               iy+sprinti("%0.2i",im)+".nc"
     fu0 = addfile(file_uv(0),"r")
     fu1 = addfile(file_uv(1),"r")
     if ( .not. isdefined("um") ) then
       lon = fu0->u&longitude_1
       lat = fu0->u&latitude_1
     end if
     ; calculate EASMI
     u = (/fu0->u_p_uv(0,{850},:,:),fu1->u_p_uv(0,{850},:,:)/)
     v = (/fu0->v_p_uv(0,{850},:,:),fu1->v_p_uv(0,{850},:,:)/)
     u!1 = "lat"
     u!2 = "lon"
     u&lat = lat
     u&lon = lon
     do ie=0, 1
       easmi(ie,n) = easmi(ie,n) + ( avg(u(ie,{22.5:32.5},{110:140})) - \
                                     avg(u(ie,{ 5.0:15.0},{ 90:130})) )/3.
     end do
   end do
   n = n + 1
 end do
print(dim_avg(easmi))

; easmi = easmi - avg(easmi)
 easmi = easmi + 3.96

;--------------------------------------- process 0

;--------------------------------------- cn level

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 plot = xylineplot(wks,time,easmi,"D",0,"T",1.2,"C",(/24,4/))
 xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
 xyaddline(plot,(/year1-6,year2+5/),0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1-1.0,year2+0.5,-999,999,False,False)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickStartF"    : 1950
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues
 axis_str(plot,"East Asian Summer Monsoon Index","year","(m s:S:-1:N:)")
 draw(plot)

 frame(wks)

 status_exit(10)

end

