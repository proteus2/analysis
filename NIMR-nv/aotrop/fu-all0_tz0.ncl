load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"tc"/);,"merra"/)
 ppid    = (/"j" /)
 varname = "u"

 y0 = 1951
 y1 = 1951
 y2 = 1965
 lat0 = 5.

 wks = cr_wks("fig/fu-all0_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   fnam1 = "../dat/AOL60CGW/epf/uan"+expname(0)+".epfsa0_yz_recon0.nc"
   fnam2 = "../dat/AOL60CGW/epf/uan"+expname(0)+".epf_yz_recon0.nc"
;   fnam1r = "/data4/kyh/dat/MERRA/epf/"+expname(1)+".epfsa0_yz_recon0.nc"
;   fnam2r = "/data4/kyh/dat/MERRA/epf/"+expname(1)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
;   f1r = addfile(fnam1r,"r")
;   f2r = addfile(fnam2r,"r")
   z  = f1->z / 1.e3
;   zr = 18.6 + 6.6*log( 70./ f1r->p )

   nz = dimsizes(z)
   nt = 240
;   nzr= dimsizes(zr)
;   ntr= 336
   nw = dimsizes(f1->wg) + 2

   dati1 = f1->epd(:,:,:,{-lat0:lat0})
   dati2 = f2->epd(:,8:9,:,{-lat0:lat0})
;   dati1r = f1r->epd(:,:,:,{-lat0-1.5:lat0+1.5})
;   dati2r = f2r->epd(:,8:9,:,{-lat0-1.5:lat0+1.5})
;   dati2(:,1,:,:) = dati2(:,1,:,:) + f2->epd(:,2,:,{-lat0:lat0}) + \
;                                     f2->epd(:,3,:,{-lat0:lat0})
;   dati2r(:,1,:,:) = dati2r(:,1,:,:) + f2r->epd(:,2,:,{-lat0-1.5:lat0+1.5}) + \
;                                       f2r->epd(:,3,:,{-lat0-1.5:lat0+1.5})

   rw_i = new((/nt-1,nw,nz/),float)
   rw_i(:,:nw-3,:) = dim_avg( dati1 )
   rw_i(:,nw-2:,:) = dim_avg( dati2 )
   rw0 = new((/nt,nw,nz/),float)
   rw0(:nt-2,:,:) = rw_i(:nt-2,:,:)
   rw0(nt-1,:,:) = rw0(nt-2,:,:)
;   rw_ir = new((/ntr-1,nw,nzr/),float)
;   rw_ir(:,:nw-3,:) = dim_avg( dati1r )
;   rw_ir(:,nw-2:,:) = dim_avg( dati2r )
;   rw0r = new((/ntr,nw,nzr/),float)
;   rw0r(1:,:,:) = rw_ir(:,:,:)
;   rw0r(0,:,:) = rw0r(1,:,:)


   z1 = f1->z(0)
   z9 = f1->z(nz-1)
;y   fnamg = "../dat/AOL60CGW/du_xm.uan"+expname+"a.pe.nc"
   fnamg = "../dat/AOL60CGW/mmxm.uan"+expname+"a.pg.nc"
   fg = addfile(fnamg,"r")
   zg = fg->hybrid_ht_1 / 1.e3

   gw = new((/nt,nz/),float)
;y   gw(:,:) = dim_avg( fg->du_gw      (:,{z1:},{-lat0:lat0}) )*72. + \
;y             dim_avg( fg->dudt_gwconv(:,{z1:},{-lat0:lat0}) )*86400.
fnamc = "../dat/AOL60CGW/mmxm.uan"+expname+"a.pd.nc"
fd = addfile(fnamc,"r")
gw(:,:) = dim_avg( fg->dudt_gwussp(:,{z1:z9},{-lat0:lat0}) )*86400. + \
          dim_avg( fd->dudt_gwconv(:,{z1:z9},{-lat0:lat0}) )*86400.


   fnam1 = "../dat/AOL60CGW/xm.uan"+expname(0)+"a.pj.nc"
;   fnam1r = "$DATD/clim/MERRA/"+expname(1)+".u_ypt.1979-2006.nc"

   f1 = addfile(fnam1,"r")
;   f1r = addfile(fnam1r,"r")
   zu = f1->hybrid_ht / 1.e3
;   zur = 18.6 + 6.6*log( 70./ f1r->levels )

   nzu = dimsizes(zu)
;   nzur = dimsizes(zur)

   u0 = transpose( dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) ) )
;   u0r = transpose( dim_avg( f1r->u(:ntr-1,:,{-lat0-1.5:lat0+1.5}) ) )

   ut0  = new((/nzu,nt-1/),float)
;   ut0r = new((/nzur,ntr-1/),float)
   ut0 (:,:) = (u0 (:,1:nt-1 ) - u0 (:,0:nt-2 ))/30.
;   ut0r(:,:) = (u0r(:,1:ntr-1) - u0r(:,0:ntr-2))/30.


   fnam2 = "../dat/AOL60CGW/tem/uan"+expname(0)+".tem_yzt.nc"
;   fnam2r = "/data18/kyh/dat/qbo/"+expname(1)+".tem_ypt.1979-2006.nc"
   f2 = addfile(fnam2,"r")
;   f2r = addfile(fnam2r,"r")
   za = f2->z / 1.e3
;   zar = 18.6 + 6.6*log( 70./ f2r->p )

   advz = transpose( dim_avg( f2->uadv_z(:,:,{-lat0:lat0}) ) )
;   advzr = transpose( dim_avg( f2r->uadv_z(:,:,{-lat0-1.5:lat0+1.5}) ) )


   wf = (/ rw0(:,2,:), rw0(:,3,:), rw0(:,4,:)+rw0(:,5,:), \
           rw0(:,7,:)+rw0(:,8,:), gw, \
           rw0(:,1,:)+rw0(:,8,:) /)

;   wfr = (/ rw0r(:,2,:), rw0r(:,3,:), rw0r(:,4,:)+rw0r(:,5,:), \
;            rw0r(:,7,:)+rw0r(:,8,:), rw0r(:,0,:)*0., \
;            rw0r(:,1,:)+rw0r(:,8,:) /)

   fu = new((/6,nz,nt/),float)
   do n=0, nt-1
     fu(:,:,n) = wf(:,n,:)
   end do
;   fur = new((/6,nzr,ntr/),float)
;   do n=0, ntr-1
;     fur(:,:,n) = wfr(:,n,:)
;   end do

   t = y0 + (fspan(1,nt,nt)-0.5)/12.
;   tr = 1979. + (fspan(1,ntr,ntr)-0.5)/12.

   tt = y0 + (fspan(1,nt-1,nt-1))/12.
;   ttr = 1979. + (fspan(1,ntr-1,ntr-1))/12.

;--------------------------------------- process 0
 fu = fu*30.
; fur = fur*30.
 ut0 = ut0*30.
; ut0r = ut0r*30.
 advz = advz*30.
; advzr = advzr*30.

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

;--------------------------------------- cn level
 cscale = (/4,4,1,1,1,1,4/)

 cnlev = fspan(-10,10,21)*0.5
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 7

 wi = 64 + (64./13.*2.)
 he = 15
 xf = 17 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 5 + ispan(0,7,1)*(he+5)

 mstr  = new((/ne,nl/),string)
 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:,:) = ""
 mstr(:,0) = (/"HG-C"/);,"MERRA"/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"z (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"dU/dt","ADV_Z","K","R","RG","IG","G"/)

do ie=0, ne-1
do il=0, nl-1

 if ( ie .eq. 1 .and. il .eq. 6 ) then
   continue
 end if

 if (il .lt. 2) then
   if (il .eq. 0) then
     if (ie .eq. 0) then
       contour1 = cnshadeplot(wks,tt,zu,ut0,cnlev*cscale(il),cnfil,lbopt(il))
     else
       contour1 = cnshadeplot(wks,ttr,zur,ut0r,cnlev*cscale(il),cnfil,"")
     end if
   end if
   if (il .eq. 1) then
     if (ie .eq. 0) then
       contour1 = cnshadeplot(wks,t,za,advz,cnlev*cscale(il),cnfil,lbopt(il))
     else
       contour1 = cnshadeplot(wks,tr,zar,advzr,cnlev*cscale(il),cnfil,"")
     end if
   end if
 else
   if (ie .eq. 0) then
     contour1 = cnshadeplot(wks,t,z,fu(il-2,:,:),cnlev*cscale(il),cnfil,lbopt(il))
   else
     contour1 = cnshadeplot(wks,tr,zr,fur(il-2,:,:),cnlev*cscale(il),cnfil,"")
   end if
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50+0.05 ;+ 0.55
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
 axis_range(contour1,y1,y2+1,17,39,False,False)
 axis_str(contour1,mstr(ie,il),xstr(il),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 if (ie .eq. 0) then
   contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 else
   contouru = cnlineplot(wks,tr,zur,u0r,cnlevu,0)
 end if
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(ie),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,17,39,False,False)
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
end do


 frame(wks)

end



