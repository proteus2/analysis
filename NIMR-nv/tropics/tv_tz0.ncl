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

 y1 = 2002
 y2 = 2006
 lat0 = 14.

 wks = cr_wks("fig/tv_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"rainbow+gray")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 n1 = (y1-1953)*12
 n2 = (y2-1953)*12 + 11
 nt = n2 - n1 + 1

 fname = "../dat/reconstr/uanuj.var_wave_yz_recon0.nc"
 fname1= "../dat/xm.uanuja.pj.nc"
 f0 = addfile(fname ,"r")
 f1 = addfile(fname1,"r")

 z = f0->z({10.e3:})
 nz = dimsizes(z)

 vv0 = ( dim_sum( f0->var_theta(n1:n2-1,:,{10.e3:},{0:lat0}) )*2. - \
         f0->var_theta(n1:n2-1,:,{10.e3:},{0}) )/ \
       (dimsizes(f0->lat({0:lat0}))*2-1)

 vv0a = new((/1,7,nz/),float)
 vv0a = vv0(nt-2,:,:)
 vv1 = array_append_record(vv0,vv0a,0)

 vk  = transpose( vv1(:,2,:) )
 vr1 = transpose( vv1(:,3,:) )
 vrge= transpose( vv1(:,4,:) )
 vrgw= transpose( vv1(:,5,:) )

 vv = (/ vk, vr1, vrge, vrgw /)

 fname2 = "$DATD/L60CGW/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 k0 = minind(abs( f2->hybrid_ht - 10.e3 ))
 k1 = k0 + 1
 if ( f2->hybrid_ht(k0) .gt. 10.e3 ) then
   k1 = k0 - 1
 end if
k1 = k0  ; temp
 p0 = dim_avg(f2->p(k0:k0+nz-1,{-lat0:lat0}))
 p1 = dim_avg(f2->p(k1:k1+nz-1,{-lat0:lat0}))
 pi = exp(0.5*(log(p0) + log(p1)))

 vv = vv * conform_dims( dimsizes(vv), ((pi/1.e5)^0.287)^2, 1)

 u = transpose( dim_avg( f1->u(n1:n2,{10.e3:},{-lat0:lat0}) ) )

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev0a = fspan(0.1,1.9,19)
 cnlev0b = fspan(0.025,0.475,19)
 cnlev1 = (/cnlev0a,cnlev0b,cnlev0b,cnlev0b/)
 cnlev2 = fspan(-50,50,11)
 cnfil = ispan(1,211,11)
print(dimsizes(cnfil))

;--------------------------------------- contour
 nl = 4

 wi = 40
 he = 30
 xf = 34 + (/0,0,1,1/)*60
 yf = 10 + (/0,1,0,1/)*37.5

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,z/1.e3,vv(il,:,:),cnlev1(il,:),cnfil,0)
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

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,y1,y2,20,50,False,False)
 axis_str(contour1,"","","altitude (km)")
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)


 contour2 = cnlineplot(wks,t,z/1.e3,u,cnlev2,0)
 cnline_att(contour2,"D",(/0,0,14/),"T",(/0.8,2.5,0.8/),"C",0)
 setvalues contour2
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

 set_size(contour2,xf(il),yf(il),wi,he)
 axis_range(contour2,y1,y2,20,50,False,False)
 axis_str(contour1,"","","altitude (km)")
 setvalues contour2
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour2,wks)

 draw(contour2)

   label = create "lab" labelBarClass wks
;     "lbAutoManage" : False
     "lbOrientation" : "Vertical"
     "lbBoxCount"   : dimsizes(cnfil)
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
;     "lbLabelsOn" : False ;True
     "lbLabelStrings" : cnlev1(il,:)
     "lbLabelStride" : 4
;     "lbAutoManage" : False
     "lbLabelAlignment" : "InteriorEdges"
     "lbPerimOn" : False
;     "lbTitleString" : "[ K:S:2:N: / cyc day:S:-1:N: ]"
;     "lbTitlePosition" : "Bottom"
   end create
   set_size_lb(label,xf(il)+wi+2,yf(il),5,he)
   draw(label)

end do

 frame(wks)

end



