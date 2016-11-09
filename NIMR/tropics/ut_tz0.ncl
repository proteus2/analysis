load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"

 y1 = 1953
 y2 = 2006
 lat0 = 5.

 wks = cr_wks("fig/ut_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   fnam1 = "../dat/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"

   f1 = addfile(fnam1,"r")
   zu = f1->hybrid_ht / 1.e3

   nzu = dimsizes(zu)

   u0 = transpose( dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) ) )

   ut0 = new((/nzu,nt-1/),float)
   ut0(:,:) = (u0(:,1:nt-1) - u0(:,0:nt-2))/30.

   t = y1 + (fspan(1,nt,nt)-0.5)/12.
   tt = y1 + (fspan(1,nt-1,nt-1))/12.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-20,20,21)*0.03
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 3
 ys = (/1953,1971,1989/)
 yl = 18
 yls = 2

 wi = 126
 he = 20
 xf = 34
 yf = 10 + (/0,1,2/)*28

 wi2 = tofloat(wi)/yl*yls

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:) = ""
 mstr(0) = cexp(0)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "z (km)"
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: day:S:-1:N:)"

do il=0, nl-1

 contour1 = cnshadeplot(wks,tt,zu,ut0,cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50+0.05
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.4
  "pmLabelBarHeightF"        : 0.15
  "lbLabelFontHeightF" : 0.009
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.012
 end setvalues

 set_size(contour1,xf,yf(il),wi,he)
 axis_range(contour1,ys(il),ys(il)+yl,16,45,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf,yf(il),wi,he)
 axis_range(contouru,ys(il),ys(il)+yl,16,45,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

if (il .gt. 0) then

  setvalues contour1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contour2 = cnshadeplot(wks,tt,zu,ut0,cnlev,cnfil,"")
  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),16,45,False,False)
  axis_str(contour2,"","",ystr(il))
  setvalues contour2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  cnfinalize(contour2,wks)
  draw(contour2)

  setvalues contouru
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contouru2 = cnlineplot(wks,t,zu,u0,cnlevu,0)
  cnline_att(contouru2,"D",(/0,0,1/),"T",0.6,"C",1)
  set_size(contouru2,xf-wi2,yf(il),wi2,he)
  axis_range(contouru2,ys(il)-yls,ys(il),16,45,False,False)
  axis_str(contouru2,"","","")
  setvalues contouru2
   "tmXBOn" : False
   "tmXTOn" : False
  end setvalues
  cnfinalize(contouru2,wks)
  draw(contouru2)

end if

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

end do


 frame(wks)

end



