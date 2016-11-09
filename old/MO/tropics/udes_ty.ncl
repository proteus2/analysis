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
 kz = (/15,18,21/)

 wks = cr_wks("udes_ty-"+expname(0),600,0)
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

   u = new((/2,nz,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/u2(n,:,:),u1(n,:,:)/)
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
 cnlev = fspan(-50,50,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 70
 he = 25
 xf = 12 + (/0,1,2,3/)*80
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/nz,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "time (mon)"
 ystr = (/"lat (deg)",""/)
 lbopt = new((/nz,2/),string)
 lbopt(:,:) = ""
 lbopt(nz-1,0) = (/"(m s:S:-1:N:)"/)
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
 axis_str(contour1,mstr(k,iv),xstr(k),ystr(iv))
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

end



