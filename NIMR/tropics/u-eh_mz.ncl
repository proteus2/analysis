load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 lat0 = (/5,10,15/)

 wks = cr_wks("fig/u-eh_mz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "$DATD/L60CGW-t/xm-eh.uan"+expname(0)+"a.pj.nc"
   f0 = addfile(fnam0,"r")
   u = f0->u(:,{1.5e3:20.e3},{-30:30})
   z = u&$u!1$/1.e3
   nz = dimsizes(z)
   nt = dimsizes(u&$u!0$)
   nrng = dimsizes(lat0)

   u3 = new((/nrng,nt,nz/),float)
   do j=0, nrng-1
     u3(j,:,:) = dim_avg(u(:,:,{-lat0(j):lat0(j)}))
   end do

   t = fspan(0,13,14)

   uss = new((/nrng,nz,14/),float)
   do n=1, 12
     uss(:,:,n) = dim_avg_n( u3(:,ispan(n-1,nt-1,12),:), 1 )
   end do
   uss(:,:,0 ) = uss(:,:,12)
   uss(:,:,13) = uss(:,:,1 )
   uss!1 = "z"
   uss&z = z

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = ispan(-16,0,2)
 cnfil = (/22,24,25,26,17,18,19,20,21,0/)
 cnfil = (/22,23,24,25,26,13,14,15,16,0/)

;--------------------------------------- contour
 wi = 30
 he = 22
 xf = 12 + (/0,1,2,3/)*38
 yf = 10 + (/0,0,0,0/)*28

 ystr  = new(nrng,string,"")
 lbopt = new(nrng,string)
 xstr = "month"
 ystr(0) = "z (km)"
 lbopt = ""
 lbopt(nrng-2) = "(mPa)"

 do j=0, nrng-1

 contour1 = cnshadeplot(wks,t,z,uss(j,:,:),cnlev,cnfil,lbopt(j))
 cnaddline(contour1,wks,t,z,uss(j,:,:),ispan(-20,20,1),0,"D",0,"T",0.5,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.01
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(j),yf(j),wi,he)
 axis_range(contour1,0.5,12.5,2,19,False,False)
 axis_str(contour1,"",xstr,ystr(j))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

; dr_txt2(wks,txtstr(j),xf(j)+0.6,yf(j)+he-2,22,6,0.,"CenterLeft",1)

end do

;  txtstrb = ":F22:"+lat0avg+":S:o:N:N:NF29:-:F22:"+lat0avg+":S:o:N:S"
;  dr_txt(wks,txtstrb,xf(3)+wi-8,yf(3)+2,22,6,0.)

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,6,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(3)-2,25,6,0.)


 frame(wks)

end



