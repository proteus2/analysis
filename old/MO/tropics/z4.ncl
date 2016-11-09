load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 30
 kz = (/15,18,21/)

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 32*12
 kz0 = (/13,15,19/)

 wks = cr_wks("z4_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   u1 = f1->u(4:nmon-1,kz,:,0)
   lat = f1->latitude
   p1  = f1->p(kz)
print(p1)

   tmp = 10.^(floattoint(log10(p1)+100.)-100)
   p1 = round(p1/tmp,0)*tmp

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = f2->u(4:nmon-1,kz,:,0)

   nz = dimsizes(kz)
   ny = dimsizes(lat)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3  = addfile(fnam0,"r")
   u3 = f3->u(m1:m1+nmon0-1,kz0,:)
   lat3= f3->latitude
   p   = f3->levels(kz0)

   ny3 = dimsizes(lat3)
   nt3 = nmon0

   u = new((/2,nz,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/u2(n,:,:),u1(n,:,:)/)
   end do

   ur = new((/nz,ny3,nt3/),float)
   do n=0, nt3-1
     ur(:,:,n) = u3(n,:,:)
   end do

   t = fspan(0,13,14)

   uss = new((/3,nz,ny,14/),float)
   do n=1, 12
     uss(0:1,:,:,n) = dim_avg(u(0:1,:,:,ispan(n-1,nt-1,12)))
;     uss(2,:,:,n) = linint1(lat3,dim_avg(ur(:,:,ispan(n-1,nt-1,12))),False,lat,0)
   end do
   uss(2,:,:,:) = uss(1,:,:,:) - uss(0,:,:,:)
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,50,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38 + (/0,0,1/)*8
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/nz,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "month"
 ystr = (/"lat (deg)","","lat (deg)"/)
 lbopt = new((/nz,3/),string)
 lbopt(:,:) = ""
 lbopt(nz-1,0) = (/"(m s:S:-1:N:)"/)
 txtstr = new((/nz,3/),string)
 txtstr(:,0) = p1(:)+" hPa"
 txtstr(:,1) = p1(:)+" hPa"
 txtstr(:,2) = p(:)+" hPa"


do k=0, nz-1
do iv=0, 2

 contour1 = cnshadeplot(wks,t,lat,uss(iv,k,:,:),cnlev,cnfil,lbopt(k,iv))
 cnaddline(contour1,wks,t,lat,uss(iv,k,:,:), \
           (/-80,-70,-60,-50,50,60,70,80/),0, \
           "D",0,"T",0.8,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.0
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(k),wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr(k,iv),xstr(k),ystr(iv))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
;  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 0
;  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt(wks,txtstr(k,iv),xf(iv)+6,yf(k)+3,22,6,0.)

end do
end do


 frame(wks)

end



