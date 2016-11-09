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
 lat0 = 30;15.
 p1 = (/30, 3, 1,0.3/)
 kz = (/55,62,64,66 /)
;  100,50,30,10, 3, 1,0.3
;   47,53,55,59,62,64,66

 tail = (/"-"+expname(0),""/)

do ie=0, 1

 wks = cr_wks("gwddes_ty"+tail(ie),600,0)
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

   t = fspan(1,nt,nt)

   uss = new((/2,nz,ny,12/),float)
   do n=0, 11
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n,nt-1,12)))
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,:,n) = u(:,:,:,n) - uss(:,:,:,mod(n,12))
   end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/-80,-60,-40,-20,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,20,40,60,80/)*0.1
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 70
 he = 25
 xf = 12 + (/0,1,2,3/)*80
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/2,nz,2/),string)
 mstr(:,:,:) = ""
 mstr(0,0,:) = (/"CGW","NGW"/)
 mstr(1,0,:) = (/"CTL","SC"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "time (mon)"
 ystr = (/"lat (deg)",""/)
 lbopt = new((/nz,2/),string)
 lbopt(:,:) = ""
 lbopt(nz-1,0) = (/"(m s:S:-1:N: day:S:-1:N:)"/)
 txtstr = new((/nz,2/),string)
 txtstr(:,:) = ""
 txtstr(:,0) = p1(:)+" hPa"


do k=0, nz-1
do iv=0, 1

 contour1 = cnshadeplot(wks,t,lat,udes(iv,k,:,:),cnlev,cnfil,lbopt(k,iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.2
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.02
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.025
 end setvalues

 set_size(contour1,xf(iv),yf(k),wi,he)
 axis_range(contour1,-999,999,-lat0,lat0,False,False)
 axis_str(contour1,mstr(ie,k,iv),xstr(k),ystr(iv))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 12
  "tmXBTickStartF"    : 0
  "tmXBMinorPerMajor" : 11
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt(wks,txtstr(k,iv),xf(iv)+6,yf(k)+3,22,6,0.)

end do
end do


 frame(wks)

end do

end



