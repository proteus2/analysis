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
 y2 = 2015
 lat0 = 15.
 z0 = 27.e3;17.41e3

 wks = cr_wks("fig/fz-kw_ty0-z"+toint(z0/1.e3)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   fnam1 = "../dat/L60CGW/epf/uan"+expname(0)+".epfsa0_yz2_recon0.nc"
   fnam2 = "../dat/L60CGW/epf/uan"+expname(0)+".epf_yz2_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   y = f1->lat({-lat0:lat0})

   ny = dimsizes(y)
   nt = 648
   nw = dimsizes(f1->wg) + 2

   rw_i = new((/647,nw,ny/),float)
   rw_i(:,:nw-3,:) = f1->f_z(:,:,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw_i(:,nw-2:,:) = f2->f_z(:,8:9,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw0 = new((/nt,nw,ny/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
   rw0(647,:,:) = rw0(646,:,:)

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,7:8,:)/)

   rw(:,:,0,:) = rw(:,:,0,:)*(-1.)

   u = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = rw(:,n,:,:)
   end do

   t = y1 + (fspan(1,nt,nt)-0.5)/12.

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

   fnam1 = "../dat/L60CGW/xm.uan"+expname(0)+"a.pj.nc"

   f1 = addfile(fnam1,"r")
   yu = f1->latitude({-lat0:lat0})

   nyu = dimsizes(yu)

   u0 = transpose( f1->u(:nt-1,{z0},{-lat0:lat0}) )

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
; cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)

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
 ystr(:) = "lat (deg)"
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: day:S:-1:N:)"

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,y,u(0,0,:,:),cnlev,cnfil,lbopt(il))
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
 axis_range(contour1,ys(il),ys(il)+yl,-999,999,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 contouru = cnlineplot(wks,t,yu,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.3,"C",1)

 set_size(contouru,xf,yf(il),wi,he)
 axis_range(contouru,ys(il),ys(il)+yl,-999,999,False,False)
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

  contour2 = cnshadeplot(wks,t,y,u(0,0,:,:),cnlev,cnfil,"")
  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),-999,999,False,False)
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

  contouru2 = cnlineplot(wks,t,yu,u0,cnlevu,0)
  cnline_att(contouru2,"D",(/0,0,1/),"T",0.3,"C",1)
  set_size(contouru2,xf-wi2,yf(il),wi2,he)
  axis_range(contouru2,ys(il)-yls,ys(il),-999,999,False,False)
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



