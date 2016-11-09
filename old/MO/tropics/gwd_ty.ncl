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
 lat0 = 60;15.
 kz = (/49,54,58,63,64/)+1
; 18,23,30,45
 kzl = 21

 wks = cr_wks("gwd_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   mfe1 = f1->unspecified(4:nmon-1,kz,:,0)
   lat = f1->latitude
   z   = f1->hybrid_ht/1.e3
print(z(kz))

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->field424(4:nmon-1,kz,:,0)

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfe3 = f3->field424(4:nmon-1,kz,:,0)

   nz = dimsizes(kz)
   ny = dimsizes(lat)
   nt = nmon-4

mfe2 = mfe3
mfe1 = mfe1+mfe2
   mf = new((/3,nz,ny,nt/),float)
   do n=0, nt-1
     mf(:,:,:,n) = (/mfe1(n,:,:),mfe2(n,:,:),mfe3(n,:,:)/)*86400.
   end do

   t = fspan(1,nt,nt)

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/-80,-60,-40,-20,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,20,40,60,80/)*0.1
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 70
 he = 25
 xf = 12 + (/0,1,2,3/)*80
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/nz,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CGW","NGW"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "time (mon)"
 ystr = (/"lat (deg)",""/)
 lbopt = new((/nz,2/),string)
 lbopt(:,:) = ""
 lbopt(nz-1,0) = (/"(m s:S:-1:N: day:S:-1:N:)"/)


do k=0, nz-1
do iv=0, 1;2

 contour1 = cnshadeplot(wks,t,lat,mf(iv,k,:,:),cnlev,cnfil,lbopt(k,iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
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

end do
end do


 frame(wks)

end



