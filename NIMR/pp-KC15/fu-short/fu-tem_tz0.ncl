load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
; cexp    = (/"HadGEM2-CGW","MERRA"/)
; expname = (/"uj","merra"/)
; varname = (/"u","u"/)
 cexp    = (/"HadGEM2-CGW","ERA-inter"/)
 expname = (/"uj","era-int"/)
 varname = (/"u","U_GDS0_ISBL"/)
 ppid    = (/"j" /)

 iplot1 = 0

 y0 = 1953
 y1 = 1979
 y2 = 1990
 lat0 = 5.
 lat0r = lat0 + 1.0

 wks = cr_wks("fig/fu-tem_tz0-"+expname(iplot1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   nt = 648
   ntr= 336

   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"
   fnam1r = "$DATD/clim/"+cexp(1)+"/"+expname(1)+".u_ypt.1979-2006.nc"
   f1 = addfile(fnam1,"r")
   f1r = addfile(fnam1r,"r")
   zu = f1->hybrid_ht / 1.e3
   zur = 18.6 + 6.6*log( 70./ f1r->$f1r->$varname(1)$!1$ )

   u0 = transpose( dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) ) )
   u0r = transpose( dim_avg( f1r->$varname(1)$(:ntr-1,:,{-lat0r:lat0r}) ) )

   fnam2 = "$DATD/L60CGW-t/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   fnam2r = "$DATD/tem/"+cexp(1)+"/"+expname(1)+".tem_ypt.1979-2006.nc"
   f2 = addfile(fnam2,"r")
   f2r = addfile(fnam2r,"r")
   za = f2->z / 1.e3
   zar = 18.6 + 6.6*log( 70./ f2r->p )

   ut = transpose( dim_avg( f2->u_tend(:,:,{-lat0:lat0}) ) )
   utr = transpose( dim_avg( f2r->u_tend(:,:,{-lat0r:lat0r}) ) )

   ue = transpose( dim_avg( f2->epd(:,:,{-lat0:lat0}) ) )
   uer = transpose( dim_avg( f2r->epd(:,:,{-lat0r:lat0r}) ) )

   advz = transpose( dim_avg( f2->uadv_z(:,:,{-lat0:lat0}) ) )
   advzr = transpose( dim_avg( f2r->uadv_z(:,:,{-lat0r:lat0r}) ) )

   advy = transpose( dim_avg( f2->uadv_y(:,:,{-lat0:lat0}) ) )
   advyr = transpose( dim_avg( f2r->uadv_y(:,:,{-lat0r:lat0r}) ) )
   advy = advy + transpose( dim_avg( f2->cor(:,:,{-lat0:lat0}) ) )
   advyr = advyr + transpose( dim_avg( f2r->cor(:,:,{-lat0r:lat0r}) ) )


   z1 = za(0)*1.e3
   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")

   gw = transpose( dim_avg( fg->du_gw      (:,{z1:},{-lat0:lat0}) )*72. + \
                   dim_avg( fg->dudt_gwconv(:,{z1:},{-lat0:lat0}) )*86400. )


   t = y0 + (fspan(1,nt,nt)-0.5)/12.
   tr = 1979. + (fspan(1,ntr,ntr)-0.5)/12.

 ureq = ut - advz - advy
 ureqr = utr - advzr - advyr

 tend = (/ureq,ue,ureq-ue,ureq-ue-gw/)
 tendr = (/ureqr,uer,ureqr-uer/)

;--------------------------------------- process 0
 tend = tend*30.
 tendr = tendr*30.

;--------------------------------------- cn level
 cscale = (/1,2,1,1/)

 cnlev = fspan(-10,10,21)*2.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 4
 if (iplot1 .eq. 1) then
   nl = 3
 end if

 wi = ( 64 + (64./13.*2.) )*12./15
 he = 15
 xf = 17 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,7,1)*(he+5)

 mstr  = new((/ne,nl/),string)
 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:,:) = ""
 mstr(:,0) = cexp
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"z (km)","z (km)"/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"F:B:U:N:_req","EPD","G+D","D"/)
 do il=0, nl-1
   if (cscale(il) .ne. 1) then
     txtstr(il) = txtstr(il)+" (:F34:4:F22: "+cscale(il)+")"
   end if
 end do

;do ie=0, ne-1
ie = iplot1
do il=0, nl-1

 if (ie .eq. 0) then
   contour1 = cnshadeplot(wks,t,za,tend(il,:,:)*cscale(il),cnlev,cnfil,lbopt(il))
 else
   contour1 = cnshadeplot(wks,tr,zar,tendr(il,:,:)*cscale(il),cnlev,cnfil,lbopt(il))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues

 set_size(contour1,xf(ie),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,17,36,False,False)
 axis_str(contour1,mstr(ie,il),xstr(il),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 if (ie .eq. 0) then
   contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 else
   contouru = cnlineplot(wks,tr,zur,u0r,cnlevu,0)
 end if
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(ie),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,17,36,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(il),xf(ie)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)

 delete(contour1)
 delete(contouru)

end do
;end do


 frame(wks)

end



