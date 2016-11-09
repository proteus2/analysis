load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 30
 p1 = (/0.3, 1/)
 kz = (/66 ,64/)
;  100,50,30,10, 3, 1,0.3
;   47,53,55,59,62,64,66

 tail = (/"-"+expname(0),""/)

do ie=0, 1

 wks = cr_wks("z3_ty"+tail(ie),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   u1 = f1->unspecified(4:nmon-1,kz,:,0)
   lat = f1->latitude
   z1  = f1->hybrid_ht(kz)/1.e3
;   p1 = 1.e3*exp(-z1/7.)
;print(z1)

;   tmp = 10.^(floattoint(log10(p1)+100.)-100)
;   p1 = round(p1/tmp,0)*tmp

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = f2->field424(4:nmon-1,kz,:,0)
   nz = dimsizes(kz)
   ny = dimsizes(lat)
   nt = nmon-4

   if (ie .eq. 1) then
     u2 = u2+u1
     fnam2 = "../dat/ng_"+expname(1)+".nc"
     f2 = addfile(fnam2,"r")
     delete(u1)
     u1 = f2->field424(4:nmon-1,kz,:,0)
   end if

   u = new((/2,nz,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/u1(n,:,:),u2(n,:,:)/)*86400.
   end do

   t = fspan(0,13,14)

   uss = new((/3,nz,ny,14/),float)
   do n=1, 12
     uss(0:1,:,:,n) = dim_avg(u(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(2,:,:,:) = uss(1,:,:,:) - uss(0,:,:,:)
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/-9,-8,-7,-6,-5,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6,7,8,9/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/2,nz,3/),string)
 mstr(:,:,:) = ""
 mstr(0,0,:) = (/"CGW","NGW",""/)
 mstr(1,0,:) = (/"CTL","SC",""/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "month"
 ystr = (/"lat (deg)","",""/)
 lbopt = new((/nz,3/),string)
 lbopt(:,:) = ""
 lbopt(nz-1,0) = (/"(m s:S:-1:N: day:S:-1:N:)"/)
 txtstr = new((/nz,3/),string)
 txtstr(:,:) = ""
 txtstr(:,0) = p1(:)+" hPa"


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
 axis_str(contour1,mstr(ie,k,iv),xstr(k),ystr(iv))
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

end do

end



