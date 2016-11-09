load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"f","0"/)

 nmon = 4+144

 wks = cr_wks("momann_yz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1  = addfile(fnam1,"r")
   mfect1 = f1->unspecified_1(4:nmon-1,:,:,0)
   mfwct1 = f1->unspecified_2(4:nmon-1,:,:,0) *(-1.)
   lat = f1->latitude
   z   = f1->hybrid_ht_1/1.e3

   fnam2 = "../dat/ng_"+expname(0)+".nc"

   f2 = addfile(fnam2,"r")
   mfeun1 = f2->field420(4:nmon-1,:,:,0)
   mfwun1 = f2->field422(4:nmon-1,:,:,0)

   mfesm1 = mfeun1+mfect1
   mfwsm1 = mfwun1+mfwct1

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfeun2 = f3->field420(4:nmon-1,:,:,0)
   mfwun2 = f3->field422(4:nmon-1,:,:,0)

   mf_avg = dim_avg_n_Wrap( \
             (/ (/mfect1, mfwct1/), \
                (/mfeun1, mfwun1/), \
                (/mfesm1, mfwsm1/), \
                (/mfeun2, mfwun2/) /), 2 )*1.e3

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220/)*0.01
 cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,19/)

;--------------------------------------- contour
 wi = 32
 he = 42
 xf = 12 + (/0,1,2,3/)*40
 yf = 10 + (/0,1/)*56

 mstr = (/ (/"CGW MF","NGW MF","Sum","NGW MF (CTL)"/), (/"","","",""/) /)
 lbopt = (/ (/"","","",""/), (/"","(mPa)","",""/) /)
 ystr = (/"model height (km)","","",""/)

do id=0, 1
do iv=0, 3

 contour1 = cnshadeplot(wks,lat,z,mf_avg(iv,id,:,:),cnlev,cnfil,lbopt(id,iv))
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
 axis_range(contour1,-999,999,0,80,False,False)
 axis_str(contour1,mstr(id,iv),"lat (deg)",ystr(iv))

 cnfinalize(contour1,wks)


 draw(contour1)

end do
end do


 frame(wks)

end



