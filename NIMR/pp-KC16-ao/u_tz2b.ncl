load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin
;--------------------------------------- parameters
 cexp    = "CGW"
 expname = "tc"
 ppid    = "j"
 varname = "u"

 y1 = 1955
 y2 = y1 - 1 + 40
 lat0 = 5.

 wks = cr_wks("fig/u_tz2b-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir  = "$DATD/AOL60"+cexp
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$
 z = tmp&$tmp!1$ / 1.e3
 nz = dimsizes(z)
 t = tmp&$tmp!0$
 t_char = tochar(t@units)
 if ( t_char(0:10) .ne. "days since " ) then
   print(t@units)
   print("Check this. It should be 'days since ...'.")
   exit
 end if
 y0i = toint(tostring(t_char(11:14)))
 m0i = toint(tostring(t_char(16:17)))
 d0i = toint(tostring(t_char(19:20)))
 print("date0 : "+y0i+" "+m0i+" "+d0i)

 t = y0i + (m0i-1)/12. + (d0i-1)/360. + t/360.
 nt = dimsizes(t)
 print("data : "+t(0)+" - "+t(nt-1))

 u = transpose( dim_avg( tmp(:,:,{-lat0:lat0}) ) )

 delete(tmp)

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnlev = fspan(-40,28,18)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 yl = 13
 yls = 1;2

 wi = 64
 he = 15
 xf = 20
 yf = 10 + ispan(0,5,1)*(he+5)

 nfig = 2
 dt = 85

do ifig=0, nfig-1

 if (ifig .gt. 0) then
   y1 = y1 + dt
   y2 = y2 + dt
   yf = yf + (he+5)*3 + 4
 end if

 ys = ispan(y1+yls,y2,yl)
 ys(0) = ys(0) - yls

 nl = dimsizes(ys)

 wi2 = tofloat(wi)/yl*yls

 xf0 = xf - wi2
 wi0 = wi + wi2

 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = ":F25:"+(/"(a)","(b)"/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = "height [km]"
 lbopt(:) = ""
 if (ifig .eq. nfig-1) then
   lbopt(nl-1) = "[m s:S:-1:N:]"
 end if


do il=0, nl-1

 contour1 = cnshadeplot(wks,t,z,u,cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,t,z,u,(/0./),0,"D",0,"T",0.8,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.18 ;0.12
  "pmLabelBarParallelPosF"   : 0.50
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.015
  "lbTitleFontHeightF" : 0.018
 end setvalues

 ye = ys(il)+yl
 xf1 = tofloat(xf)
 wi1 = tofloat(wi)
 if (il .eq. 0) then
   ye = ye + yls
   xf1 = xf0
   wi1 = wi0
 end if
 if (ye .gt. y2+1) then
   ye = y2+1
   wi1 = wi*(tofloat(y2+1-ys(il))/yl)
 end if

 set_size(contour1,xf1,yf(il),wi1,he)
 axis_range(contour1,ys(il),ye,17,39,False,False)
 axis_str(contour1,"",xstr(il),ystr)

 ps = (/100,70,50,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf1+wi1+2,yf(il),1,he)
 axis_range(plotx,-999,999,17,39,False,False)
 axis_str(plotx,"","","p [hPa]")
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

 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues


 if (il .gt. 0) then

  setvalues contour1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contour2 = cnshadeplot(wks,t,z,u,cnlev,cnfil,"")
  cnaddline(contour2,wks,t,z,u,(/0./),0,"D",0,"T",0.8,"C",1)

  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),17,39,False,False)
  axis_str(contour2,"","",ystr)

  setvalues contour2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues

  cnfinalize(contour2,wks)

  draw(contour2)

 end if


 cnfinalize(contour1,wks)

 draw(contour1)
 draw(plotx)

end do  ; il

 if (nfig .gt. 1) then
   dr_txt2(wks,mstr(ifig),xf-wi2-11,yf(0)-1,22,6,0.,"CenterLeft",1)
 end if

end do  ; ifig

 frame(wks)

end

