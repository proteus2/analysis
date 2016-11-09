load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+12*12
 nmonc= 4+12*12;-12
 lat0 = 15.

 wks = cr_wks("fig/gwdvar_fz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/cg_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->unspecified(4:nmon-1,:,{-lat0:lat0},0) )*86400.
   lat1 = f1->latitude({-lat0:lat0})
   z1 = f1->hybrid_ht/1.e3

   fnam1 = "../dat/ng_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u2 = dim_avg( f1->field424(4:nmon-1,:,{-lat0:lat0},0) )*86400.

   u3 = u1 + u2

   fnam2 = "../dat/ng_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u0 = dim_avg( f2->field424(4:nmonc-1,:,{-lat0:lat0},0) )*86400.

   ny = dimsizes(lat1)
   nz = dimsizes(z1)
   nt = nmon-4
   ntc= nmonc-4

   uc = new((/nz,ntc/),float)
   do n=0, ntc-1
     uc(:,n) = u0(n,:)
   end do
   u = new((/3,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u1(n,:),u2(n,:),u3(n,:)/)
   end do

 ss = new((/4,nz,12/),float)
 do n=0, 11
   ss(0:2,:,n) = dim_avg(u(:,:,ispan(n,nt-1 ,12)))
   ss(3,:,n) = dim_avg(uc(:,ispan(n,ntc-1,12)))
 end do

 do n=0, nt-1
   u(0:2,:,n) = u(0:2,:,n) - ss(0:2,:,mod(n,12))
 end do
 do n=0, ntc-1
   uc(:,n) = uc(:,n) - ss(3,:,mod(n,12))
 end do

   nfc = ntc/2
   nf  = nt/2
   frc = fspan(1,nfc,nfc)/(int2flt(ntc)/12.)
   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
print(frc(0)+"  "+fr(0))

   fftuc = ezfftf(uc)
   fftu  = ezfftf(u )

   pwsc = ( fftuc(0,:,:)^2 + fftuc(1,:,:)^2 )*0.5*ntc/12.
   pws  = ( fftu (0,:,:,:)^2 + fftu (1,:,:,:)^2 )*0.5*nt /12.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/1,2,5,10,20,50,100,200/)*0.1
 cnfil = (/0,11,9,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 42
 he = 24
 xf = 12 + (/0,1,0,1/)*50 ;+ (/0,0,0,1/)*8
 yf = 10 + (/0,0,1,1/)*45

 mstr = new((/3,4/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CGW forcing","NGW forcing","Sum","NGW forcing (CTL)"/)
 xstr = (/"","","","frequency (cycle yr:S:-1:N:)"/)
 ystr = (/"model height (km)","","model height (km)",""/)
 lbopt = new((/3,4/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"","","(10:S:-1:N: (m s:S:-1:N: day:S:-1:N:):S:2:N: / cycle yr:S:-1:N:)",""/)


do iv=0, 3

 if (iv .lt. 3) then
   contour1 = cnshadeplot(wks,fr,z1,pws(iv,:,:)*10.,cnlev,cnfil,lbopt(2,iv))
 else
   contour1 = cnshadeplot(wks,frc,z1,pwsc*10.,cnlev,cnfil,lbopt(2,iv))
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
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,0,4,15,65,False,False)
 axis_str(contour1,mstr(0,iv),xstr(3),ystr(iv))

 cnfinalize(contour1,wks)


 draw(contour1)

end do


 frame(wks)

end



