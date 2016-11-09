load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+12*12
 lat0 = 30.
 kz = 47

 wks = cr_wks("mf100ss_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   mfe1 = f1->unspecified_1(4:nmon-1,kz,{-lat0:lat0},0)
   mfw1 = f1->unspecified_2(4:nmon-1,kz,{-lat0:lat0},0)*(-1.)
   lat1 = f1->latitude({-lat0:lat0})

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->field420(4:nmon-1,kz,{-lat0:lat0},0)
   mfw2 = f2->field422(4:nmon-1,kz,{-lat0:lat0},0)
   mfe2 = mfe2 + mfe1
   mfw2 = mfw2 + mfw1

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfe3 = f3->field420(4:nmon-1,kz,{-lat0:lat0},0)
   mfw3 = f3->field422(4:nmon-1,kz,{-lat0:lat0},0)

   ny = dimsizes(lat1)
   nt = nmon-4

   u = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/ (/mfe1(n,:),mfw1(n,:)/), \
                     (/mfe2(n,:),mfw2(n,:)/), \
                     (/mfe3(n,:),mfw3(n,:)/) /)*1.e3
   end do

   t = fspan(0,13,14)

   uss = new((/3,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
; cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1/)*38
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-MF","Westward-MF"/)
 xstr = (/"","","month"/)
 ystr = (/"lat (deg)",""/)
 txtstr = (/"CGW", "CGW + NGW", "NGW (CTL)"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(mPa)",""/)


do iv=0, 2
do idr=0, 1

 contour1 = cnshadeplot(wks,t,lat1,uss(iv,idr,:,:),cnlev,cnfil,lbopt(iv,idr))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.18
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.4
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))
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

  dr_txt2(wks,txtstr(iv),xf(0)+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do
end do

  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,8,0.)
  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end



