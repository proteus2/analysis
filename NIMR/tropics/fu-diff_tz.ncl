load "$UTIL/ncl/header.ncl"

xspace = 80   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = (/"uj","tc"/)
 datname = (/"L60CGW-t","AOL60CGW"/)
 uname = (/"u","u"/)

 ie = 1
 y0 = 1953
 y1 = 2005 - 15
 y2 = 2010 - 15
 lat0 = 5 ;10
 ztop = 35

 wks = cr_wks("fig/fu-diff_tz-"+y1+"-"+y2+"-ya"+lat0+"-"+expname(ie),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12

 fnam0 = "$DATD/"+datname(ie)+"/xm.uan"+expname(ie)+"a.pj.nc"
 fnam3 = "$DATD/"+datname(ie)+"/tem/uan"+expname(ie)+".tem-wc_yzt.nc"
 fnamg = "$DATD/"+datname(ie)+"/xm.uan"+expname(ie)+"a.pe.nc"

 latrng = (/-lat0,lat0/)

   f0 = addfile(fnam0,"r")
   f3 = addfile(fnam3,"r")
   fg = addfile(fnamg,"r")

   lat = f3->lat({-lat0:lat0})
   dy05 = 0.5*(lat(1) - lat(0))
   wgt = where( lat .ge. latrng(0)+dy05 .and. lat .le. latrng(1)-dy05, 1., 0.)
   wgt = where( abs(lat-latrng(0)) .lt. dy05, 0.5*(1.+(lat-latrng(0))/dy05), wgt)
   wgt = where( abs(lat-latrng(1)) .lt. dy05, 0.5*(1.-(lat-latrng(1))/dy05), wgt)

   z = f3->z({15.e3:}) / 1.e3
   nz = dimsizes(z)

   print(expname(ie)+"  "+abs(lat(0))+"  "+wgt(0))

   zname = f0->$uname(ie)$!1
   u0 = transpose( dim_avg( f0->$uname(ie)$(n0:n0+nt-1,{15.e3:},{-lat0:lat0}) ) )
   zu = f0->hybrid_ht({15.e3:}) / 1.e3
;   za = f3->z / 1.e3

   ut   = transpose( dim_avg_wgt_Wrap( f3->u_tend (n0:n0+nt-1,{15.e3:},{-lat0:lat0}), wgt,1) )
   uf   = transpose( dim_avg_wgt_Wrap( f3->u_force(n0:n0+nt-1,{15.e3:},{-lat0:lat0}), wgt,1) )

   gw = transpose( dim_avg( fg->du_gw      (n0:n0+nt-1,{15.e3:},{-lat0:lat0}) )*72. + \
                   dim_avg( fg->dudt_gwconv(n0:n0+nt-1,{15.e3:},{-lat0:lat0}) )*86400. )

   u1 = transpose( dim_avg_wgt_Wrap( f0->$uname(ie)$(n0:n0+nt-1,{15.e3:},{-lat0:lat0}), wgt, 1) )
   dz = 1.2e3

   u4z = new((/nz,nt/),float)
   do k=2, nz-3
     u4z(k,:) = (u1(k-2,:)-4.*u1(k-1,:)+6.*u1(k,:)-4.*u1(k+1,:)+u1(k+2,:)) / dz^4
   end do
   u4z = -u4z*86400. * 1.5e6

;   ureq = (/ ut - advz - advy /)

;   if (ie .eq. 0) then
     fu = new((/2,nz,nt/),float)
     t = y1 + (fspan(1,nt,nt)-0.5)/12.
;     tt = y1 + (fspan(1,nt-1,nt-1))/12.
;   end if
   fu = (/ ut - uf - gw, u4z /)

; end do  ; ie

;--------------------------------------- process 0
 fu = fu*30.

;--------------------------------------- cn level
 cscale = (/1,1/)

 cnlev = fspan(-15,15,21)

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-50,50,11)

;--------------------------------------- contour
 nv = 2

 wi = 60 * tofloat(y2-y1+1)/8.
 he = 15
 xf = 10
 yf = 10 + ispan(0,nv-1,1)*(he+8)

 mstr  = new(nv,string,"")
 xstr  = new(nv,string,"")
 lbopt = new(nv,string)

 mstr(0) = "";expname(ie)
; xstr(nv-1) = "year"
 xstr(1) = "year"
 ystr = "height (km)"
 lbopt(:) = ""
; lbopt(nv-1) = "(m s:S:-1:N: month:S:-1:N:)"
 lbopt(1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"D",":F18:-:F22:K:B:4:N::F18:s:F22::S:4:N:U/:F18:s:F22:z:S:4:N:"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)
 do iv=0, nv-1
   if (cscale(iv) .ne. 1) then
     txtstr(iv) = txtstr(iv)+" (:F34:4:F22: "+cscale(iv)+")"
   end if
 end do

do iv=0, 1 ;nv-1

 contour1 = cnshadeplot(wks,t,z,fu(iv,:,:)*cscale(iv),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12 ;+ 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.65
  "pmLabelBarHeightF"        : 0.25
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.016
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05 + 0.05
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf,yf(iv),wi,he)
 axis_range(contour1,y1,y2+1,19,ztop,False,False)
 axis_str(contour1,mstr(iv),xstr(iv),ystr)
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",(/0.6,1.5,0.6/),"C",1)

 set_size(contouru,xf,yf(iv),wi,he)
 axis_range(contouru,y1,y2+1,19,ztop,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(iv),xf+1,yf(iv)-2,22,5,0.,"CenterLeft",1)
; dr_txt(wks,subl(iv),xf(0)-8,yf(iv),25,8,0.)

 delete(contour1)
 delete(contouru)

 ps = (/100,70,50,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf+wi+2,yf(iv),1,he)
 axis_range(plotx,-999,999,19,ztop,False,False)
 axis_str(plotx,"","","p (hPa)")
 setvalues plotx
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : True
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tiYAxisSide" : "Right"
  "tiYAxisAngleF" : 270
  "tmYRLabelsOn" : True
  "tmYLMode" : "Explicit"
  "tmYLValues" : zs
  "tmYLLabels" : ps
 end setvalues

 draw(plotx)

end do


 frame(wks)

end

