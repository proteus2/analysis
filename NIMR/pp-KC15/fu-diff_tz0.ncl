load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"

 y0 = 1953
 y1 = 1980
 y2 = 1991
 lat0 = 5.
 lat0r = lat0 + 1.0

 wks = cr_wks("fig/fu-diff_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   nt = 648

   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")
   z = fg->hybrid_ht_1 / 1.e3
   nz = dimsizes(z)

   gw = new((/nt,nz/),float)
   gw(:,:) = dim_avg( fg->du_gw      (:,:,{-lat0:lat0}) )*72. + \
             dim_avg( fg->dudt_gwconv(:,:,{-lat0:lat0}) )*86400.

   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"

   f1 = addfile(fnam1,"r")

   u0 = transpose( dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) ) )
   u0y = f1->u(:nt-1,:,{-lat0-2.5:lat0+2.5})
   u0z = dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) )

   ny = dimsizes(f1->latitude({-lat0:lat0}))
   dy = 1.25/180.*3.141592*6371.e3
   dz = 1.2e3

   u4y0 = new((/nt,nz,ny/),float)
   do j=0, ny-1
     u4y0(:,:,j) = (u0y(:,:,j)-4.*u0y(:,:,j+1)+6.*u0y(:,:,j+2)-4.*u0y(:,:,j+3)+u0y(:,:,j+4)) / dy^4
   end do
   u4y = -dim_avg(u4y0)*86400. * 2.e16

   u4z = new((/nt,nz/),float)
   do k=2, nz-3
     u4z(:,k) = (u0z(:,k-2)-4.*u0z(:,k-1)+6.*u0z(:,k)-4.*u0z(:,k+1)+u0z(:,k+2)) / dz^4
   end do
   u4z = -u4z*86400. * 1.5e6

   u2z = new((/nt,nz/),float)
   do k=1, nz-2
     u2z(:,k) = (u0z(:,k-1)-2.*u0z(:,k)+u0z(:,k+1)) / dz^2
   end do
   u2z0 = u2z*86400. * 0.8
   u2z = u2z*86400. * 0.5

   fnam2 = "$DATD/L60CGW-t/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   f2 = addfile(fnam2,"r")
   za = f2->z / 1.e3

   ut = transpose( dim_avg_Wrap( f2->u_tend(:,:,{-lat0:lat0}) ) )

   advz = transpose( dim_avg_Wrap( f2->uadv_z(:,:,{-lat0:lat0}) ) )

   uf = transpose( dim_avg_Wrap( f2->u_force(:,:,{-lat0:lat0}) ) )

   ue = transpose( dim_avg_Wrap( f2->epd(:,:,{-lat0:lat0}) ) )

   advy = transpose( dim_avg_Wrap( f2->uadv_y(:,:,{-lat0:lat0}) ) )
   advy = advy + transpose( dim_avg_Wrap( f2->cor(:,:,{-lat0:lat0}) ) )

   ureq = ut

   ureq = (/ ut - advz - advy /)

   wf = (/ gw,u4y,u4z,u2z0,u2z, \
           gw /)

   fu = new((/6,nz,nt/),float)
   do n=0, nt-1
     fu(:,:,n) = wf(:,n,:)
   end do
print(dimsizes(fu))
print(dimsizes(ureq))

   fu!1 = "z"
   fu&z = z*1.e3

   fu(0,{12.e3:},:) = ureq({12.e3:},:) - ue({12.e3:},:) - fu(0,{12.e3:},:)

   t = y0 + (fspan(1,nt,nt)-0.5)/12.

;   tt = y0 + (fspan(1,nt-1,nt-1))/12.
;   ttr = 1979. + (fspan(1,ntr-1,ntr-1))/12.

;--------------------------------------- process 0
 fu = fu*30.
 ut = ut*30.
 advz = advz*30.
 ureq = ureq*30.

 u3 = (/ut, advz, ureq/)

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*1.
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 5

 wi = 64 + (64./13.*2.)
 he = 15
 xf = 22 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,nl-1,1)*(he+5)

 xstr  = new(nl,string)
 lbopt = new(nl,string)

 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"z (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"D:B:u:H-10V15:_:N:","Uyyyy","Uzzzz","Uzz","Uzz2"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,z,fu(il,:,:),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12 + 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.55
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05 + 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues

 set_size(contour1,xf(0),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,"",xstr(il),ystr(0))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
;  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 contouru = cnlineplot(wks,t,z,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(0),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,19,33,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(il),xf(0)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)
 dr_txt(wks,subl(il),xf(0)-8,yf(il),25,8,0.)

 delete(contour1)
 delete(contouru)

end do


 frame(wks)

end

