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
 nmonc= 4+12*12
 lat0 = 30.
 kz = 47

 wks = cr_wks("mf100_fy-"+expname(0),600,0)
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
   mfe3 = f3->field420(4:nmonc-1,kz,{-lat0:lat0},0)
   mfw3 = f3->field422(4:nmonc-1,kz,{-lat0:lat0},0)

   ny = dimsizes(lat1)
   nt = nmon-4
   ntc= nmonc-4

   u = new((/2,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/ (/mfe1(n,:),mfw1(n,:)/), \
                     (/mfe2(n,:),mfw2(n,:)/) /)*1.e3
   end do
   uc = new((/2,ny,ntc/),float)
   do n=0, ntc-1
     uc(:,:,n) = (/mfe3(n,:),mfw3(n,:)/)*1.e3
   end do

   nfc = ntc/2
   nf  = nt/2
   frc = fspan(1,nfc,nfc)/(int2flt(ntc)/12.)
   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
print(frc(0)+"  "+fr(0))

   fftuc = ezfftf(uc)
   fftu  = ezfftf(u )

   pwsc = ( fftuc(0,:,:,:)^2 + fftuc(1,:,:,:)^2 )*0.5*ntc/12.
   pws  = ( fftu (0,:,:,:,:)^2 + fftu (1,:,:,:,:)^2 )*0.5*nt /12.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/2,5,10,20,50,100,200,500/)*0.01
 cnfil = (/0,9,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1/)*38
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-MF","Westward-MF"/)
 xstr = (/"","","frequency (cycle yr:S:-1:N:)"/)
 ystr = (/"lat (deg)",""/)
 txtstr = (/"CGW", "CGW + NGW", "NGW (CTL)"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(mPa:S:2:N: / cycle yr:S:-1:N:)",""/)


do iv=0, 2
do idr=0, 1

 if (iv .lt. 2) then
   contour1 = cnshadeplot(wks,fr,lat1,pws(iv,idr,:,:),cnlev,cnfil,lbopt(iv,idr))
 else
   contour1 = cnshadeplot(wks,frc,lat1,pwsc(idr,:,:),cnlev,cnfil,lbopt(iv,idr))
 end if
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.0
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,-999,2.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))

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



