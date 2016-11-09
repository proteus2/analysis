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
 year1 = (/1951,1951,1953,1953,1979,1957/)
 year2 = (/1979,1979,2005,2005,2005,2001/)

 wks = cr_wks("fig/eawmi_all-"+jobname(1),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- data reading
 ne = dimsizes(jobname)
 year1m = min(year1)
 year2m = max(year2)
 nt = year2m - year1m + 1
 time = fspan(year1m,year2m,nt)

 eawmi = new((/ne,nt/),float)
 eawmi = 0.

 mon3 = (/12,1,2/)

do ie=0, ne-3
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=0, 2
     file_uv = "/data18/HG2"+exptype(ie)+"/"+expname(ie)+"/pj/uan"+ \
               jobname(ie)+"a.pj_"+iy+sprinti("%0.2i",mon3(im))+".nc"
     if (mon3(im) .le. 2) then
       file_uv = "/data18/HG2"+exptype(ie)+"/"+expname(ie)+"/pj/uan"+ \
                 jobname(ie)+"a.pj_"+(iy+1)+sprinti("%0.2i",mon3(im))+".nc"
     end if
     fu0 = addfile(file_uv,"r")
     ; calculate EAWMI
     u = fu0->u_p_uv(0,{300},:,:)
     v = fu0->v_p_uv(0,{300},:,:)
     eawmi(ie,n) = eawmi(ie,n) + ( avg(u({27.5:37.5},{110:170})) - \
                                   avg(u({50.0:60.0},{ 80:140})) )/3.
   end do
   n = n + 1
 end do
end do

 ie = ne - 2
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=0, 2
     file_u = "diag_eawmi/dat/"+iy+"/merra.u_m.assm.pl"+300+"."+ \
              iy+sprinti("%0.2i",mon3(im))+".nc"
     file_v = "diag_eawmi/dat/"+iy+"/merra.v_m.assm.pl"+300+"."+ \
              iy+sprinti("%0.2i",mon3(im))+".nc"
     if (mon3(im) .le. 2) then
       file_u = "diag_eawmi/dat/"+(iy+1)+"/merra.u_m.assm.pl"+300+"."+ \
                (iy+1)+sprinti("%0.2i",mon3(im))+".nc"
       file_v = "diag_eawmi/dat/"+(iy+1)+"/merra.v_m.assm.pl"+300+"."+ \
                (iy+1)+sprinti("%0.2i",mon3(im))+".nc"
     end if
     fur0 = addfile(file_u,"r")
     fvr0 = addfile(file_v,"r")
     ; calculate EAWMI
     ur = fur0->u(0,0,:,:)
     vr = fvr0->v(0,0,:,:)
     eawmi(ie,n) = eawmi(ie,n) + ( avg(ur({27.5:37.5},{110:170})) - \
                                   avg(ur({50.0:60.0},{ 80:140})) )/3.
   end do
   n = n + 1
 end do

 ie = ne - 1
 n = year1(ie) - year1m
 do iy=year1(ie), year2(ie)
   do im=0, 2
     file_u = "u3_"+sprinti("%0.2i",mon3(im))+".nc"
     fur0 = addfile(file_u,"r")
     ; calculate EAWMI
     ur2 = fur0->u(iy-1957,:,:)
     eawmi(ie,n) = eawmi(ie,n) + tofloat( avg(ur2({27.5:37.5},{110:170})) - \
                                          avg(ur2({50.0:60.0},{ 80:140})) )/3.*tofloat(ur2@scale_factor)
   end do
   n = n + 1
 end do

 eawmi@_FillValue = 0.

; eawmi = eawmi - avg(eawmi)
 eawmi = eawmi - avg(eawmi(ne-2,:))
print(dim_avg(eawmi))

;--------------------------------------- process 0
 m1 = max(year1(2:4)) - year1m
 m2 = nt - 1 - (year2m - min(year2(2:4)))
 cor2 = sum((eawmi(2,m1:m2)-avg(eawmi(2,m1:m2)))*eawmi(4,m1:m2))/ \
        sqrt( sum((eawmi(2,m1:m2)-avg(eawmi(2,m1:m2)))^2)* \
              sum(eawmi(4,m1:m2)^2) )
 cor3 = sum((eawmi(3,m1:m2)-avg(eawmi(3,m1:m2)))*eawmi(4,m1:m2))/ \
        sqrt( sum((eawmi(3,m1:m2)-avg(eawmi(3,m1:m2)))^2)* \
              sum(eawmi(4,m1:m2)^2) )
 print(cor2+"  "+cor3)

 m1 = max(year1((/0,1,5/))) - year1m
 m2 = nt - 1 - (year2m - min(year2((/0,1,5/))))
 cor0 = sum((eawmi(0,m1:m2)-avg(eawmi(0,m1:m2)))*eawmi(5,m1:m2))/ \
        sqrt( sum((eawmi(0,m1:m2)-avg(eawmi(0,m1:m2)))^2)* \
              sum(eawmi(5,m1:m2)^2) )
 cor1 = sum((eawmi(1,m1:m2)-avg(eawmi(1,m1:m2)))*eawmi(5,m1:m2))/ \
        sqrt( sum((eawmi(1,m1:m2)-avg(eawmi(1,m1:m2)))^2)* \
              sum(eawmi(5,m1:m2)^2) )
 print(cor0+"  "+cor1)

;--------------------------------------- cn level

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/30,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 plot = xylineplot(wks,time,eawmi,"D",(/0,0,2,2,0,0/),"T",(/1.6,1.6,0.8,0.8,1.0,1.0/),"C",(/4,2,4,2,1,3/))
 xyaddmarker(plot,16,"S",0.25,"T",1.0,"C",1)
 xyaddline(plot,(/year1m-6,year2m+5/),0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1m-1.0,year2m+0.5,-999,999,False,False)
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(1940,2100,5)
  "tmXBLabels" : sprinti("%2.2i",mod(ispan(1940,2100,5),100))+"/"+ \
                 sprinti("%2.2i",mod(ispan(1941,2101,5),100))
  "tmXBMinorValues" : ispan(1940,2100,1)
 end setvalues
 axis_str(plot,"East Asian Winter Monsoon Index","year","(m s:S:-1:N:)")
 draw(plot)

 frame(wks)

 status_exit(10)

end

