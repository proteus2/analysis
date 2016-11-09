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
 kz = (/45,49,54,58,63/)
; 15,18,23,30,45
 kzl = 21

 wks = cr_wks("mom_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   mfe1 = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw1 = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   lat = f1->latitude
   z   = f1->hybrid_ht_1/1.e3
print(z(kz))

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->field420(4:nmon-1,kz,:,0)
   mfw2 = f2->field422(4:nmon-1,kz,:,0)

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfe3 = f3->field420(4:nmon-1,kz,:,0)
   mfw3 = f3->field422(4:nmon-1,kz,:,0)

   nz = dimsizes(kz)
   ny = dimsizes(lat)
   nt = nmon-4

mfe2 = mfe3
mfw2 = mfw3
;mfe1 = mfe1+mfe2
;mfw1 = mfw1+mfw2
   mf = new((/3,2,nz,ny,nt/),float)
   do n=0, nt-1
     mf(:,:,:,:,n) = \ ;dim_avg_n_Wrap( \
                     (/ (/mfe1(n,:,:), mfw1(n,:,:)/), \
                        (/mfe2(n,:,:), mfw2(n,:,:)/), \
                        (/mfe3(n,:,:), mfw3(n,:,:)/) /)*1.e3
   end do

   t = fspan(1,nt,nt)

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300/)*0.01
 cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18 ,17,0/)

;--------------------------------------- contour
 wi = 70
 he = 25
 xf = 12 + (/0,1,2,3/)*80
 yf = 10 + (/0,1,2,3,4/)*30

 lbopt = (/ (/"","","",""/), (/"","(mPa)","",""/) /)
 mstr = new((/nz,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CGW","NGW"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "time (mon)"
 ystr = (/"lat (deg)",""/)

do idr=0, 0;1
do k=0, nz-1
do iv=0, 1;2

 contour1 = cnshadeplot(wks,t,lat,mf(iv,idr,k,:,:),cnlev,cnfil,lbopt(idr,iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarParallelPosF"   : 0.7
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.2
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : -0.02
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
end do


 frame(wks)

end



