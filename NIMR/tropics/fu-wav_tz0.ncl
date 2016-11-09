load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"HadGEM2-CGW","MERRA","ERA-inter"/)
 expname = (/"uj","merra","era-int"/)
 varname = (/"u","u","U_GDS0_ISBL"/)
 ppid    = (/"j" /)

 iplot1 = 0

 y0 = 1953
 y1 = 1980
 y2 = 1991
 lat0 = 5.
 lat0r = lat0 + 1.

 wks = cr_wks("fig/fu-wav_tz0-"+expname(iplot1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ie = iplot1

 lat00 = lat0
 if (ie .ne. 0) then
   lat00 = lat0r
 end if

 if (ie .eq. 0) then
   nt = 648
   t = y0 + (fspan(1,nt,nt)-0.5)/12.
 else
   nt = 336
   t = 1979. + (fspan(1,nt,nt)-0.5)/12.
 end if

 if (ie .eq. 0) then
   fnam1 = "$DATD/L60CGW-t/epf-wc/uan"+expname(ie)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc/uan"+expname(ie)+".epf_yz_recon0.nc"
 else
   fnam1 = "$DATD/tem/"+cexp(ie)+"/epf/"+expname(ie)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+cexp(ie)+"/epf/"+expname(ie)+".epf_yz_recon0.nc"
 end if
 f1 = addfile(fnam1,"r")
 f2 = addfile(fnam2,"r")

 if (ie .eq. 0) then
   z = f1->z / 1.e3
 else
   z = 18.6 + 6.6*log( 70./ f1->p )
 end if
 nz = dimsizes(z)
 nw = 8

   rw_i = new((/nt-1,nw,nz/),float)
   rw_i(:,:,:) = dim_avg( f1->epd(:,:nw-1,:,{-lat00:lat00}) )
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     dim_avg( f2->epd(:,8:9,:,{-lat00:lat00}) )

; migrating tide
;   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
;                    ( dim_avg( f2->epd(:,2,:,{-lat00:lat00}) ) + \
;                      dim_avg( f2->epd(:,3,:,{-lat00:lat00}) ) )

   rw_i(:,3,:) = rw_i(:,3,:) + dim_avg( f1->epd(:,13,:,{-lat00:lat00}) )
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + dim_avg( f1->epd(:,16:17,:,{-lat00:lat00}) )

   rw0 = new((/nt,nw,nz/),float)
   if (ie .eq. 0) then
     rw0(:nt-2,:,:) = rw_i(:,:,:)
     rw0(nt-1,:,:) = rw0(nt-2,:,:)
   else
     rw0(1:,:,:) = rw_i(:,:,:)
     rw0(0,:,:) = rw0(1,:,:)
   end if

 if (ie .eq. 0) then
   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(ie)+"a.pj.nc"
   fnam2 = "$DATD/L60CGW-t/tem/uan"+expname(ie)+".tem-wc_yzt.nc"
 else
   fnam1 = "$DATD/clim/"+cexp(ie)+"/"+expname(ie)+".u_ypt.1979-2006.nc"
   fnam2 = "$DATD/tem/"+cexp(ie)+"/"+expname(ie)+".tem_ypt.1979-2006.nc"
 end if
 f1 = addfile(fnam1,"r")
 f2 = addfile(fnam2,"r")

 if (ie .eq. 0) then
   zu = f1->hybrid_ht / 1.e3
   za = f2->z / 1.e3
 else
   zu = 18.6 + 6.6*log( 70./ f1->$f1->$varname(ie)$!1$ )
   za = 18.6 + 6.6*log( 70./ f2->p )
 end if

 u0 = transpose( dim_avg( f1->$varname(ie)$(:nt-1,:,{-lat00:lat00}) ) )

 ut   = transpose( dim_avg( f2->u_tend(:,:,{-lat00:lat00}) ) )
 advz = transpose( dim_avg( f2->uadv_z(:,:,{-lat00:lat00}) ) )
 advy = transpose( dim_avg( f2->uadv_y(:,:,{-lat00:lat00}) ) )
 advy = advy + transpose( dim_avg( f2->cor(:,:,{-lat00:lat00}) ) )

 if (ie .eq. 0) then
   z1 = z(0)*1.e3
   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(ie)+"a.pe.nc"
   fg = addfile(fnamg,"r")

   gw = new((/nt,nz/),float)
   gw(:,:) = dim_avg( fg->du_gw      (:,{z1:},{-lat00:lat00}) )*72. + \
             dim_avg( fg->dudt_gwconv(:,{z1:},{-lat00:lat00}) )*86400.
 end if

 wf = (/ rw0(:,2,:), rw0(:,3,:), rw0(:,4,:)+rw0(:,5,:), \
         rw0(:,6,:)+rw0(:,7,:), rw0(:,0,:)*0., \
         rw0(:,0,:)+rw0(:,1,:) /)
 if (ie .eq. 0) then
   wf(4,:,:) = gw
 end if

 fu = new((/6,nz,nt/),float)
 do n=0, nt-1
   fu(:,:,n) = wf(:,n,:)
 end do

;--------------------------------------- process 0
 fu = fu*30.
 fr = (ut - advz - advy)*30.

;--------------------------------------- cn level
 cscale = (/4,4,8,8,1/)

 cnlev = fspan(-10,10,21)*2.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 6
 if (ie .ne. 0) then
   nl = 5
 end if

 wi = 64 + (64./13.*2.)
 he = 15
 xf = 17 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,7,1)*(he+5)

 mstr  = new((/dimsizes(cexp),nl/),string)
 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:,:) = ""
 mstr(:,0) = cexp
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"z (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"REQ","Kelvin","Rossby","MRG","IG","CGW+BGW"/)
 do il=1, nl-1
   if (cscale(il-1) .ne. 1) then
     txtstr(il) = txtstr(il)+" (:F34:4:F22: "+cscale(il-1)+")"
   end if
 end do

do il=0, nl-1

 if (il .eq. 0) then
   contour1 = cnshadeplot(wks,t,za,fr,cnlev,cnfil,lbopt(il))
 else
   contour1 = cnshadeplot(wks,t,z,fu(il-1,:,:)*cscale(il-1),cnlev,cnfil,lbopt(il))
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

 set_size(contour1,xf(0),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,mstr(ie,il),xstr(il),ystr(0))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
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

 delete(contour1)
 delete(contouru)

end do


 frame(wks)

end



