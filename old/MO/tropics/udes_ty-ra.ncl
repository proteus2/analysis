load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 year1 = 1979
 year2 = 2010

 m1 = 0*12
 nmon = 16*12
 lat0 = 30;15.
 kz = (/7,9,13,15,19/)

 wks = cr_wks("udes_ty",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"

   f1  = addfile(fnam1,"r")
   u1 = f1->u(m1:m1+nmon-1,kz,:)
   lat = f1->latitude
   p   = f1->levels(kz)
print(p)

   nz = dimsizes(kz)
   ny = dimsizes(lat)
   nt = nmon

   u = new((/nz,ny,nt/),float)
   do n=0, nt-1
     u(:,:,n) = u1(n,:,:)
   end do

   t = fspan(1,nt,nt)

   uss = new((/nz,ny,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,50,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 140
 he = 25*0.66
 xf = 12 + (/0,1,2,3/)*80
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new(nz,string)
 mstr(:) = ""
 mstr(0) = "MERRA"
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "time (mon)"
 ystr = "lat (deg)"
 lbopt = new(nz,string)
 lbopt(:) = ""
 lbopt(nz-1) = "(m s:S:-1:N:)"
 txtstr = p+" hPa"


do k=0, nz-1

 contour1 = cnshadeplot(wks,t,lat,udes(k,:,:),cnlev,cnfil,lbopt(k))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarParallelPosF"   : 0.5
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.01
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.012
 end setvalues

 set_size(contour1,xf,yf(k),wi,he)
 axis_range(contour1,-999,999,-lat0,lat0,False,False)
 axis_str(contour1,mstr(k),xstr(k),ystr)
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 12
  "tmXBTickStartF"    : 0
  "tmXBMinorPerMajor" : 11
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt(wks,txtstr(k),xf+10,yf(k)+3,22,6,0.)

end do


 frame(wks)

end



