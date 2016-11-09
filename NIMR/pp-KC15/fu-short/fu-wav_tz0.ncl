load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"HadGEM2-CGW","MERRA"/)
 expname = (/"uj","merra"/)
 varname = (/"u","u"/)
; cexp    = (/"HadGEM2-CGW","ERA-inter"/)
; expname = (/"uj","era-int"/)
; varname = (/"u","U_GDS0_ISBL"/)
 ppid    = (/"j" /)

 iplot1 = 0

 y0 = 1953
 y1 = 1979
 y2 = 1990
 lat0 = 5.
 lat0r = lat0 + 1.0

 wks = cr_wks("fig/fu-wav_tz0-"+expname(iplot1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   fnam1 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"
   fnam1r = "$DATD/tem/"+cexp(1)+"/epf/"+expname(1)+".epfsa_yz_recon0.nc"
   fnam2r = "$DATD/tem/"+cexp(1)+"/epf/"+expname(1)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f1r = addfile(fnam1r,"r")
   f2r = addfile(fnam2r,"r")
   z  = f1->z / 1.e3
   zr = 18.6 + 6.6*log( 70./ f1r->p )

   nz = dimsizes(z)
   nt = 648
   nzr= dimsizes(zr)
   ntr= 336
   nw = 8

   rw_i = new((/nt-1,nw,nz/),float)
   rw_i(:,:,:) = dim_avg( f1->epd(:,:nw-1,:,{-lat0:lat0}) )
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     dim_avg( f2->epd(:,8:9,:,{-lat0:lat0}) )
   rw_ir = new((/ntr-1,nw,nzr/),float)
   rw_ir(:,:,:) = dim_avg( f1r->epd(:,:nw-1,:,{-lat0r:lat0r}) )
   rw_ir(:,nw-2:,:) = rw_ir(:,nw-2:,:) + \
                     dim_avg( f2r->epd(:,8:9,:,{-lat0r:lat0r}) )

; migrating tide
;   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
;                    ( dim_avg( f2->epd(:,2,:,{-lat0:lat0}) ) + \
;                      dim_avg( f2->epd(:,3,:,{-lat0:lat0}) ) )
;   rw_ir(:,nw-1,:) = rw_ir(:,nw-1,:) + \
;                    ( dim_avg( f2r->epd(:,2,:,{-lat0r:lat0r}) ) + \
;                      dim_avg( f2r->epd(:,3,:,{-lat0r:lat0r}) ) )

   rw_i(:,3,:) = rw_i(:,3,:) + dim_avg( f1->epd(:,13,:,{-lat0:lat0}) )
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + dim_avg( f1->epd(:,16:17,:,{-lat0:lat0}) )
   rw_ir(:,3,:) = rw_ir(:,3,:) + dim_avg( f1r->epd(:,13,:,{-lat0r:lat0r}) )
   rw_ir(:,6:7,:) = rw_ir(:,6:7,:) + dim_avg( f1r->epd(:,16:17,:,{-lat0r:lat0r}) )

   rw0 = new((/nt,nw,nz/),float)
   rw0(:nt-2,:,:) = rw_i(:nt-2,:,:)
   rw0(nt-1,:,:) = rw0(nt-2,:,:)
   rw0r = new((/ntr,nw,nzr/),float)
   rw0r(1:,:,:) = rw_ir(:,:,:)
   rw0r(0,:,:) = rw0r(1,:,:)


   z1 = f1->z(0)
   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")
   zg = fg->hybrid_ht_1 / 1.e3

   gw = new((/nt,nz/),float)
   gw(:,:) = dim_avg( fg->du_gw      (:,{z1:},{-lat0:lat0}) )*72. + \
             dim_avg( fg->dudt_gwconv(:,{z1:},{-lat0:lat0}) )*86400.


   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"
   fnam1r = "$DATD/clim/"+cexp(1)+"/"+expname(1)+".u_ypt.1979-2006.nc"

   f1 = addfile(fnam1,"r")
   f1r = addfile(fnam1r,"r")
   zu = f1->hybrid_ht / 1.e3
   zur = 18.6 + 6.6*log( 70./ f1r->$f1r->$varname(1)$!1$ )

   nzu = dimsizes(zu)
   nzur = dimsizes(zur)

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

   advz = transpose( dim_avg( f2->uadv_z(:,:,{-lat0:lat0}) ) )
   advzr = transpose( dim_avg( f2r->uadv_z(:,:,{-lat0r:lat0r}) ) )

   advy = transpose( dim_avg( f2->uadv_y(:,:,{-lat0:lat0}) ) )
   advyr = transpose( dim_avg( f2r->uadv_y(:,:,{-lat0r:lat0r}) ) )
   advy = advy + transpose( dim_avg( f2->cor(:,:,{-lat0:lat0}) ) )
   advyr = advyr + transpose( dim_avg( f2r->cor(:,:,{-lat0r:lat0r}) ) )

   wf = (/ rw0(:,2,:), rw0(:,3,:), rw0(:,4,:)+rw0(:,5,:), \
           rw0(:,6,:)+rw0(:,7,:), gw, \
           rw0(:,0,:)+rw0(:,1,:) /)

   wfr = (/ rw0r(:,2,:), rw0r(:,3,:), rw0r(:,4,:)+rw0r(:,5,:), \
            rw0r(:,6,:)+rw0r(:,7,:), rw0r(:,0,:)*0., \
            rw0r(:,0,:)+rw0r(:,1,:) /)

   fu = new((/6,nz,nt/),float)
   do n=0, nt-1
     fu(:,:,n) = wf(:,n,:)
   end do
   fur = new((/6,nzr,ntr/),float)
   do n=0, ntr-1
     fur(:,:,n) = wfr(:,n,:)
   end do

   t = y0 + (fspan(1,nt,nt)-0.5)/12.
   tr = 1979. + (fspan(1,ntr,ntr)-0.5)/12.

;--------------------------------------- process 0
 fu = fu*30.
 fur = fur*30.
 fr = (ut - advz - advy)*30.
 frr = (utr - advzr - advyr)*30.

;--------------------------------------- cn level
 cscale = (/4,4,8,8,1/)

 cnlev = fspan(-10,10,21)*2.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 6
 if (iplot1 .eq. 1) then
   nl = 5
 end if

 wi = ( 64 + (64./13.*2.) )*12./15.
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
 txtstr = (/"F:B:U:N:_req","K","R","MRG","IG","G"/)
 do il=1, nl-1
   if (cscale(il-1) .ne. 1) then
     txtstr(il) = txtstr(il)+" (:F34:4:F22: "+cscale(il-1)+")"
   end if
 end do

;do ie=0, ne-1
ie = iplot1
do il=0, nl-1

 if ( ie .eq. 1 .and. il .eq. 6 ) then
   continue
 end if

 if (il .eq. 0) then
   if (ie .eq. 0) then
     contour1 = cnshadeplot(wks,t,za,fr,cnlev,cnfil,lbopt(il))
   else
     contour1 = cnshadeplot(wks,tr,zar,frr,cnlev,cnfil,lbopt(il))
   end if
 else
   if (ie .eq. 0) then
     contour1 = cnshadeplot(wks,t,z,fu(il-1,:,:)*cscale(il-1),cnlev,cnfil,lbopt(il))
   else
     contour1 = cnshadeplot(wks,tr,zr,fur(il-1,:,:)*cscale(il-1),cnlev,cnfil,lbopt(il))
   end if
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



