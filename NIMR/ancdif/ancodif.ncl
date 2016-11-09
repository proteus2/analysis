load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.5

begin

;--------------------------------------- parameters
 expname = (/"l38","l60"/)

 wks = cr_wks("ancodif",600,0)
 gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "waterfix_"+expname+".nc"

   f10 = addfile(fnam1(0),"r")
   f11 = addfile(fnam1(1),"r")

   pe_38 = f10->field672(0,0,:,:)
   pe_60 = f11->field672(0,0,:,:)
   var = (/pe_38,pe_60/)

   lon = f10->longitude

   vara = dim_avg(var)

   c3 = 2.75
   ny = dimsizes(f10->latitude)
   lat = new(ny,float)
   caa = 30./(c3-1.)
   cbb = 47.5/log(c3)
   lat(0:60) = fspan(-90.,-30.,61)
   lat(ny-61:ny-1) = fspan(30.,90.,61)
   lat(ny/2:ny-62) = caa*(exp(fspan(0.5,46.5,47)/cbb)-1.)
   lat(61:ny/2-1) = -caa*(exp(fspan(46.5,0.5,47)/cbb)-1.)

;--------------------------------------- process 0
 cnlev = fspan(1,20,20)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = ispan(36,201,8)

;--------------------------------------- plot
 xmax = (/2.5/)
 wi = (/50,50,30/)
 he = 28
 xf = 15 + (/0,1,2/)*60
 yf = 12 + (/0,1,2,3,4/)*42 + (/0,1,1,1,1/)*15

 xstr = (/"lat (deg)"/)
 ystr = (/"height (km)","",""/)
 xstr2 = (/"P-E flux"/)
 unit2 = " "+(/"(10:S:-6:N: kg m:S:-2:N: s:S:-1:N:)"/)
 lbopt = (/"","(10:S:-7:N: kg m:S:-2:N: s:S:-1:N:)"/)

 do iv=0, 0

   plot = new(3,graphic)

   do i=0, 1
     plot(i) = cnshadeplot(wks,lon,lat,var(i,:,:)/1.e-7,cnlev,cnfil,lbopt(i))
   end do
   plot(2) = xylineplot(wks,vara/1.e-6,lat,"D",0,"T",1.5,"C",(/1,170/))

   do i=0, 2
     set_size(plot(i),xf(i),yf(iv),wi(i),he)
     if (i .le. 1) then
       axis_range(plot(i),0.,360.,-90.,90.,False,False)
       axis_str(plot(i),"",xstr(iv),ystr(i))
       setvalues plot(i)
         
       end setvalues
       cnfinalize(plot(i),wks)
     else
       axis_range(plot(i),0.,xmax(iv),-90.,90.,False,False)
       axis_str(plot(i),"",xstr2(iv)+unit2(iv),ystr(i))
     end if
   end do

   draw(plot)

 end do

 dr_txt(wks,"Flux correction field ::  Salinity flux (P - E)",xf(1)+10,yf(0)-8,22,10,0.)

 frame(wks)

end



