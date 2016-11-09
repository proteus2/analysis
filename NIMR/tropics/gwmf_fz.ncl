load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 lat0 = 5.

 wks = cr_wks("fig/gwmf_fz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/mmxm.uanuja.pd.nc"
   fnam2 = "../dat/mmgwusspz_xm.uanuja.pg.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat = f1->latitude({-lat0:lat0})
   z  = f1->hybrid_ht / 1.e3

   u10 = dim_avg( f1->F_e_gwconv(0:215,:,{-lat0:lat0}) )*1.e3
   u11 = dim_avg( f1->F_e_gwconv(216:431,:,{-lat0:lat0}) )*1.e3
   u12 = dim_avg( f1->F_e_gwconv(432:,:,{-lat0:lat0}) )*1.e3

   u20 = dim_avg( f1->F_w_gwconv(0:215,:,{-lat0:lat0}) )*1.e3 *(-1.)
   u21 = dim_avg( f1->F_w_gwconv(216:431,:,{-lat0:lat0}) )*1.e3 *(-1.)
   u22 = dim_avg( f1->F_w_gwconv(432:,:,{-lat0:lat0}) )*1.e3 *(-1.)

   u10b = dim_avg( f2->Fp_e_gwussp(0:215,:,{-lat0:lat0}) )*1.e3
   u11b = dim_avg( f2->Fp_e_gwussp(216:431,:,{-lat0:lat0}) )*1.e3
   u12b = dim_avg( f2->Fp_e_gwussp(432:,:,{-lat0:lat0}) )*1.e3

   u20b = dim_avg( f2->Fp_w_gwussp(0:215,:,{-lat0:lat0}) )*1.e3 *(-1.)
   u21b = dim_avg( f2->Fp_w_gwussp(216:431,:,{-lat0:lat0}) )*1.e3 *(-1.)
   u22b = dim_avg( f2->Fp_w_gwussp(432:,:,{-lat0:lat0}) )*1.e3 *(-1.)

   u10 = u10b + u10
   u11 = u11b + u11
   u12 = u12b + u12
   u20 = u20b + u20
   u21 = u21b + u21
   u22 = u22b + u22

   ny = dimsizes(lat)
   nz = dimsizes(z)
   nt = dimsizes(f1->t)/3

   u = new((/2,3,2,nz,nt/),float)
   do k=0, nz-1
     u(:,:,:,k,:) = (/ (/ (/u10(:,k),u20(:,k)/), (/u11(:,k),u21(:,k)/), (/u12(:,k),u22(:,k)/) /), \
                       (/ (/u10b(:,k),u20b(:,k)/), (/u11b(:,k),u21b(:,k)/), (/u12b(:,k),u22b(:,k)/) /) /)
   end do

;--------------------------------------- process 0
 u = dtrend(u ,False)

 win = sin(3.141592*fspan(0,1,nt))        ; cosine window
win = 1.

 mswin = avg((win)^2)

 do n=0, nt-1
   u(:,:,:,:,n) = u(:,:,:,:,n)*win(n)
 end do

   nf  = nt/2
   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
;print(fr(0)+"  "+fr3(0))

   fftu  = ezfftf(u )

   pws0 = ( fftu (0,:,:,:,:,:)^2 + fftu (1,:,:,:,:,:)^2 )*0.5*nt /12. / mswin

   pws = dim_avg_n(pws0, 1)

 pws = sqrt(pws*fr(0))/fr(0)

;--------------------------------------- cn level
 cnlev = (/0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6/)
 cnfil = (/0,6,4,2,26,24,22,21,19,17/)

;--------------------------------------- contour
 wi = 45
 he = 24
 xf = 12 + (/0,1,2/)*53 + (/0,0,0/)*6.6
 yf = 10 + (/0,1,2,3,4/)*30

 varn = "";(/"Small-scale GW"/)
 mstr = (/ (/"Eastward-GW MF","Westward-GW MF"/), (/"",""/) /)
 xstr = (/"frequency (cycle yr:S:-1:N:)",""/)
 ystr = (/"height (km)","",""/) ;"log-p height (km)"/)
 lbopt = new((/2,2/),string)
 lbopt(:,:) = ""
 lbopt(0,:) = (/"(mPa / cycle yr:S:-1:N:)",""/) ;(/"(m s:S:-1:N:)","",""/)


do ie=0, 0
do iv=0, 1

 contour1 = cnshadeplot(wks,fr,z,pws(ie,iv,:,:),cnlev,cnfil,lbopt(ie,iv))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.2 ;1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
  "pmLabelBarHeightF"        : 0.2
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.020
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.25
  "lbTitleFontHeightF" : 0.024
 end setvalues

 set_size(contour1,xf(iv),yf(ie),wi,he)
 axis_range(contour1,-999,2.5,10,50,False,False)
 axis_str(contour1,mstr(ie,iv),xstr(ie),ystr(iv))
 setvalues contour1
  "tmXBMaxTicks" : 5
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt2(wks,varn(ie),xf(iv)+wi-1,yf(ie)+2,22,6,0.,"CenterRight",1)

end do
end do

 frame(wks)

end



