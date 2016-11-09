load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 jobname = (/"tk","tc","ui","uj","merra","era"/)
 expname = (/"AO","AO","",""/)+"L60"+(/"CTL","CGW","CTL","CGW-t"/)
 exptype = (/"CMIP","CMIP","AMIP","AMIP"/)
 year1 = (/1951,1951,1953,1953,1979,1958/)
 year2 = (/1980,1980,2006,2006,2006,2002/)

 wks = cr_wks("fig/easmi_all-"+jobname(1),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- data reading
 ne = dimsizes(jobname)
 year1m = min(year1)
 year2m = max(year2)
 nt = year2m - year1m + 1
 time = fspan(year1m,year2m,nt)

 easmi = new((/ne,nt/),float)
 easmi = 0.

do ie=0, ne-3
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=6, 8
     file_uv = "/data18/HG2"+exptype(ie)+"/"+expname(ie)+"/pj/uan"+ \
               jobname(ie)+"a.pj_"+iy+sprinti("%0.2i",im)+".nc"
     fu0 = addfile(file_uv,"r")
     ; calculate EASMI
     u = fu0->u_p_uv(0,{850},:,:)
     v = fu0->v_p_uv(0,{850},:,:)
     easmi(ie,n) = easmi(ie,n) + ( avg(u({22.5:32.5},{110:140})) - \
                                   avg(u({ 5.0:15.0},{ 90:130})) )/3.
   end do
   n = n + 1
 end do
end do

 ie = ne - 2
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=6, 8
     file_u = "diag_easmi/dat/"+iy+"/merra.u_m.assm.pl"+850+"."+ \
              iy+sprinti("%0.2i",im)+".nc"
     file_v = "diag_easmi/dat/"+iy+"/merra.v_m.assm.pl"+850+"."+ \
              iy+sprinti("%0.2i",im)+".nc"
     fur0 = addfile(file_u,"r")
     fvr0 = addfile(file_v,"r")
     ; calculate EASMI
     ur = fur0->u(0,0,:,:)
     vr = fvr0->v(0,0,:,:)
     easmi(ie,n) = easmi(ie,n) + ( avg(ur({22.5:32.5},{110:140})) - \
                                   avg(ur({ 5.0:15.0},{ 90:130})) )/3.
   end do
   n = n + 1
 end do

 ie = ne - 1
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=6, 8
     file_u = "u8_"+sprinti("%0.2i",im)+".nc"
     fur0 = addfile(file_u,"r")
     ; calculate EASMI
     ur2 = fur0->u(iy-1958,:,:)
     easmi(ie,n) = easmi(ie,n) + tofloat( avg(ur2({22.5:32.5},{110:140})) - \
                                          avg(ur2({ 5.0:15.0},{ 90:130})) )/3.*tofloat(ur2@scale_factor)
   end do
   n = n + 1
 end do

 easmi@_FillValue = 0.

; easmi = easmi - avg(easmi)
 easmi = easmi - avg(easmi(ne-2,:))
print(dim_avg(easmi))

;--------------------------------------- process 0
 m1 = max(year1(2:4)) - year1m
 m2 = nt - 1 - (year2m - min(year2(2:4)))
 cor2 = sum((easmi(2,m1:m2)-avg(easmi(2,m1:m2)))*easmi(4,m1:m2))/ \
        sqrt( sum((easmi(2,m1:m2)-avg(easmi(2,m1:m2)))^2)* \
              sum(easmi(4,m1:m2)^2) )
 cor3 = sum((easmi(3,m1:m2)-avg(easmi(3,m1:m2)))*easmi(4,m1:m2))/ \
        sqrt( sum((easmi(3,m1:m2)-avg(easmi(3,m1:m2)))^2)* \
              sum(easmi(4,m1:m2)^2) )
 print(cor2+"  "+cor3)

 m1 = max(year1((/0,1,5/))) - year1m
 m2 = nt - 1 - (year2m - min(year2((/0,1,5/))))
 cor0 = sum((easmi(0,m1:m2)-avg(easmi(0,m1:m2)))*easmi(5,m1:m2))/ \
        sqrt( sum((easmi(0,m1:m2)-avg(easmi(0,m1:m2)))^2)* \
              sum(easmi(5,m1:m2)^2) )
 cor1 = sum((easmi(1,m1:m2)-avg(easmi(1,m1:m2)))*easmi(5,m1:m2))/ \
        sqrt( sum((easmi(1,m1:m2)-avg(easmi(1,m1:m2)))^2)* \
              sum(easmi(5,m1:m2)^2) )
 print(cor0+"  "+cor1)

;--------------------------------------- cn level

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/30,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 plot = xylineplot(wks,time,easmi,"D",(/0,0,2,2,0,0/),"T",(/1.6,1.6,0.8,0.8,1.0,1.0/),"C",(/4,2,4,2,1,3/))
 xyaddmarker(plot,16,"S",0.25,"T",1.0,"C",1)
 xyaddline(plot,(/year1m-6,year2m+5/),0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1m-1.0,year2m+0.5,-999,999,False,False)
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

