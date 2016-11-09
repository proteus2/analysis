load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fa","0"/)

 nmon = 4+24

 kz = 47
 wks = cr_wks("mom_yt-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   mfect1 = f1->unspecified_1(3:nmon,kz,:,0)
   mfwct1 = f1->unspecified_2(3:nmon,kz,:,0) *(-1.)
   lat = f1->latitude
   print(f1->hybrid_ht_1(kz)/1.e3+" km")

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   mfeun1 = f2->field420(3:nmon,kz,:,0)
   mfwun1 = f2->field422(3:nmon,kz,:,0)

   mfesm1 = mfeun1+mfect1
   mfwsm1 = mfwun1+mfwct1

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfeun2 = f3->field420(3:nmon,kz,:,0)
   mfwun2 = f3->field422(3:nmon,kz,:,0)

   mf_avg = (/ (/mfect1, mfwct1/), \
               (/mfeun1, mfwun1/), \
               (/mfesm1, mfwsm1/), \
               (/mfeun2, mfwun2/) /)*1.e3

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220/)*0.01
 cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,19/)

;--------------------------------------- contour
 wi = 32
 he = 32
 xf = 12 + (/0,1,3,2/)*40
 yf = 10 + (/0,1/)*46

 mstr = (/ (/"CGW MF","NGW MF","Sum","NGW MF (CTL)"/), (/"","","",""/) /)
 lbopt = (/ (/"","","",""/), (/"","(mPa)","",""/) /)
 ystr = (/"time (month)","","",""/)

 t = fspan(0,nmon-3,nmon-2)

do id=0, 1
do iv=0, 3

 contour1 = cnshadeplot(wks,lat,t,mf_avg(iv,id,:,:),cnlev,cnfil,lbopt(id,iv))
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
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(id),wi,he)
 axis_range(contour1,-999,999,0.5,nmon-4+0.5,False,False)
 axis_str(contour1,mstr(id,iv),"lat (deg)",ystr(iv))

 cnfinalize(contour1,wks)

if (iv .ne. 2) then
 draw(contour1)
end if

end do
end do


 frame(wks)

end



