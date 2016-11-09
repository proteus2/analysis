load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expn = (/"ERA-I","era-int","ome","W_GDS0_ISBL"/)
; expn = (/"MERRA","merra","omega","omega"/)
 ftem = "$DATD/tem/"+expn(0)+"/"+expn(1)+".tem_ypt.nc"
 p = 70
 y0 = 1979
 y1 = 1980
 y2 = 2009
 h_scale = 6.6
 lat_r = 85  ; 52
 opt_w = 0 ;1
 fw = "$DATD/clim/"+expn(0)+"/"+expn(1)+"."+expn(2)+"_ypt.nc"
 str_w = expn(3) ;"W_GDS0_ISBL"
 
 wks = cr_wks("fig/ws_my-p"+p+"_"+y1+"-"+y2+"-"+expn(0),600,0)
 gsn_define_colormap(wks,"MPL_BrBG")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 n1 = (y1 - y0)*12
 n2 = (y2 - y0)*12 + 11

 f1 = addfile(ftem,"r")
 ws = f1->w_res(:,{p},{-lat_r:lat_r})

 lat = f1->$ws!1$({-lat_r:lat_r})

 nv = 1
 if (opt_w .ne. 0) then
   nv = 2
   f2 = addfile(fw,"r")
   w = f2->$str_w$(:,{p},{-lat_r:lat_r})
   w = w/(p*100.)*(h_scale*1.e3)*(-1.)
 end if

;--------------------------------------- process 0
 mon = ispan(0,13,1)

 varp = new((/nv,dimsizes(lat),14/),float)
 do n=1, 12
   varp(0,:,n) = dim_avg_n( ws(ispan(n1+n-1,n2,12),:), 0 )
 end do
 if (opt_w .ne. 0) then
   do n=1, 12
     varp(1,:,n) = dim_avg_n( w(ispan(n1+n-1,n2,12),:), 0 )
   end do
 end if

 varp(:,:,0 ) = varp(:,:,12)
 varp(:,:,13) = varp(:,:,1 )

;--------------------------------------- cn level
 cnlev = ispan(-10,10,1)*0.1
 cnfil = toint(fspan(129,2,22))
 cnlev2 = ispan(-50,50,1)*0.1

;--------------------------------------- contour
 he = 40
 wi = 32
 xf = 10 + (/0,1/)*(wi+6)
 yf = 10 + (/0,0/)*(he+6)

 mstr = (/"W*","W"/)
 xstr = "month"
 ystr = (/"lat [deg]",""/)
 lbopt = (/"[mm s:S:-1:N:]",""/)

do il=0, nv-1

 contour1 = cnshadeplot(wks,mon,lat,varp(il,:,:)*1.e3,cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,mon,lat,varp(il,:,:)*1.e3,cnlev2,0,"D",0,"T",(/0.1,2.5,0.1/),"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
  "pmLabelBarHeightF"        : 0.08
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.024
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.028
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0.5,12.5,-999,999,False,False)
 axis_str(contour1,mstr(il),xstr,ystr(il))
 setvalues contour1
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

end do

 frame(wks)

end



