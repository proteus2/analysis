load "$UTIL/ncl/header.ncl"

xspace = 80   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = (/"merra","jra55","cfsr","era-int","era-int_f"/)
 datname0= (/"MERRA","JRA55","CFSR","ERA-I"  ,"ERA-I-nr"/)
 datname = (/"MERRA","JRA55","CFSR","ERA-I"  ,"ERA-I_ml"/)
 uname = (/"u","U_GRD_GDS0_ISBL","UGRD","U_GDS0_ISBL","u"/)

 ie = 4
 feb_1st = False
 y0 = 1979
 y1 = 2005 - 15 ;2003
 y2 = 2010 - 15
 lat0 = 5 ;10
 ztop = 35

 wks = cr_wks("fig/fu-tem_tz-"+y1+"-"+y2+"-ya"+lat0+"-"+expname(ie),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nw = 8
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12
 ntr = nt
 if ( feb_1st .and. ( n0 .eq. 0 ) ) then
   ntr = nt - 1
 end if

 fnam0 = "$DATD/clim/"+datname0(ie)+"/"+expname(ie)+".u_ypt.nc"
 fnam1 = "$DATD/tem/"+datname0(ie)+"/epf/"+expname(ie)+".epfsa_yz_recon0.nc"
 fnam2 = "$DATD/tem/"+datname0(ie)+"/epf/"+expname(ie)+".epf_yz_recon0.nc"
 fnam3 = "$DATD/tem/"+datname0(ie)+"/"+expname(ie)+".tem_ypt.nc"

 latrng = (/-lat0,lat0/)

   f0 = addfile(fnam0,"r")
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")

   lat = f1->lat({-lat0:lat0})
   dy05 = 0.5*(lat(1) - lat(0))
   wgt = where( lat .ge. latrng(0)+dy05 .and. lat .le. latrng(1)-dy05, 1., 0.)
   wgt = where( abs(lat-latrng(0)) .lt. dy05, 0.5*(1.+(lat-latrng(0))/dy05), wgt)
   wgt = where( abs(lat-latrng(1)) .lt. dy05, 0.5*(1.-(lat-latrng(1))/dy05), wgt)

   z = 18.6 + 6.6*log( 70./ f1->p({100:}) )
   nz = dimsizes(z)

   print(expname(ie)+"  "+abs(lat(0))+"  "+wgt(0))

   rw_i = new((/ntr,nw,nz/),float)
   rw_i(:,:,:) = dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,:nw-1,{100:},{-lat0:lat0}), wgt,1)
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                   dim_avg_wgt_Wrap( f2->epd(n0:n0+ntr-1,8:9,{100:},{-lat0:lat0}), wgt,1)
; migrating tide
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                  ( dim_avg_wgt_Wrap( f2->epd(n0:n0+ntr-1,2,{100:},{-lat0:lat0}) + \
                            f2->epd(n0:n0+ntr-1,3,{100:},{-lat0:lat0}) , wgt,1) )

   rw_i(:,3  ,:) = rw_i(:,3  ,:) + \
                 dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,13   ,{100:},{-lat0:lat0}), wgt,1)
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + \
                 dim_avg_wgt_Wrap( f1->epd(n0:n0+ntr-1,16:17,{100:},{-lat0:lat0}), wgt,1)

   rw0 = new((/nw,nz,nt/),float,1.e20)
   do iw=0, nw-1
     rw0(iw,:,nt-ntr:nt-1) = transpose(rw_i(:,iw,:))
   end do

   zname = f0->$uname(ie)$!1
   if (expname(ie) .eq. "era-int_f") then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{:25},{-lat0:lat0}) ), wgt,1) )
     pa = asciiread("/data11/data-arch/ERA-I-nr/plev/plev_ps1000_ERA-I",60,"float")
     zu = 18.6 + 6.6*log( 70./ pa(2:24) )
   end if
   if ( expname(ie) .eq. "era-int" .or. expname(ie) .eq. "jra55" .or. \
        expname(ie) .eq. "cfsr" ) then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{:100},{-lat0:lat0}) ), wgt,1) )
     zu = 18.6 + 6.6*log( 70./ f0->$zname$({:100}) )
   end if
   if (expname(ie) .eq. "merra") then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{100:},{-lat0:lat0}) ), wgt,1) )
     zu = 18.6 + 6.6*log( 70./ f0->$zname$({100:}) )
   end if

   ut   = transpose( dim_avg_wgt_Wrap( f3->u_tend (n0:n0+nt-1,{100:},{-lat0:lat0}), wgt,1) )
   advz = transpose( dim_avg_wgt_Wrap( f3->uadv_z (n0:n0+nt-1,{100:},{-lat0:lat0}), wgt,1) )
   uf   = transpose( dim_avg_wgt_Wrap( f3->u_force(n0:n0+nt-1,{100:},{-lat0:lat0}), wgt,1) )
   ue   = transpose( dim_avg_wgt_Wrap( f3->epd    (n0:n0+nt-1,{100:},{-lat0:lat0}), wgt,1) )
   advy = transpose( dim_avg_wgt_Wrap( f3->uadv_y (n0:n0+nt-1,{100:},{-lat0:lat0}) + \
                                       f3->cor    (n0:n0+nt-1,{100:},{-lat0:lat0}) , wgt,1) )

;   ureq = (/ ut - advz - advy /)

;   if (ie .eq. 0) then
     fu = new((/8,nz,nt/),float)
     t = y1 + (fspan(1,nt,nt)-0.5)/12.
;     tt = y1 + (fspan(1,nt-1,nt-1))/12.
;   end if
   fu = (/ rw0(2,:,:), rw0(4,:,:)+rw0(5,:,:), \
           rw0(6,:,:)+rw0(7,:,:), rw0(3,:,:), ue, advz, uf, ut - uf /)

; end do  ; ie

 if ( (nt - ntr) .ne. 0 ) then
   fu(0:3,:,:,0) = fu@_FillValue
;   fu(0:3,:,:,0) = fu(0:3,:,:,1)
 end if

;--------------------------------------- process 0
 fu = fu*30.

;--------------------------------------- cn level
 cscale = (/1,3,1,1,1,1,1,1,1,1/)

 cnlev = fspan(-15,15,21)
; if (datname0(ie) .eq. "ERA-I-nr") then
;   cnlev = fspan(-15,15,21)
; end if

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-50,50,11)

;--------------------------------------- contour
 nv = 8

 wi = 60 * tofloat(y2-y1+1)/8.
 he = 15
 xf = 10
 yf = 10 + ispan(0,nv-1,1)*(he+8)

 mstr  = new(nv,string,"")
 xstr  = new(nv,string,"")
 lbopt = new(nv,string)

 mstr(0) = "";datname(ie)
; xstr(nv-1) = "year"
 xstr(3) = "year"
 ystr = "height (km)"
 lbopt(:) = ""
; lbopt(nv-1) = "(m s:S:-1:N: month:S:-1:N:)"
 lbopt(3) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"Kelvin","MRG","IG","Rossby","Net_resolved","ADVz","F:B:U:N:_total","X"/) ;":F18:s:F22:U/:F18:s:F22:t - F:B:U:N:_total"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)
 do iv=0, nv-1
   if (cscale(iv) .ne. 1) then
     txtstr(iv) = txtstr(iv)+" (:F34:4:F22: "+cscale(iv)+")"
   end if
 end do

do iv=0, 3 ;nv-1

 contour1 = cnshadeplot(wks,t,z,fu(iv+4,:,:)*cscale(iv+4),cnlev,cnfil,lbopt(iv))
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

 dr_txt2(wks,txtstr(iv+4),xf+1,yf(iv)-2,22,5,0.,"CenterLeft",1)
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

