load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

it = 30
;--------------------------------------- parameters
 fdir = "$DATD/clim/ERA-I/6hrly"
; n1s = toint((fspan(4,8,5)-4)*4)

 wks = cr_wks("fig/xmean-era",600,0)
; gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
cnt = 0
do yr=1990,1990;1990, 2009
 print(yr+"")
do m=1, 1 ;1, 12
 cnt = cnt + 1

 mm = sprinti("%2.2i",m)

 fname1 = fdir+"/"+yr+"/era-int.u_ypt."+yr+mm+"_6hrly.nc"
 fname2 = fdir+"/"+yr+"/era-int.t_ypt."+yr+mm+"_6hrly.nc"
 fname4 = fdir+"/"+yr+"/era-int.gp_ypt."+yr+mm+"_6hrly.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f4 = addfile(fname4,"r")

; nt = dimsizes(n1s)

 y = tofloat( f1->g0_lat_2({17:83}) )
 ny = dimsizes(y)

 p = tofloat( f1->lv_ISBL1 )
 nz = dimsizes(p)

 u = f1->U_GDS0_ISBL(it:it,:,{-17:-83})
 t = f2->T_GDS0_ISBL(it:it,:,{-17:-83})
 z = f4->Z_GDS0_ISBL(it:it,:,{-17:-83})/9.8

 pt = t
 do k=0, nz-1
   pt(:,k,:) = t(:,k,:)*(1.e3/p(k))^0.286
 end do

 nbv = t
 lnpt = log( pt )
 do k=1, nz-2
   nbv(:,k,:) = 9.8*(lnpt(:,k+1,:) - lnpt(:,k-1,:)) / \
                    (z   (:,k+1,:) - z   (:,k-1,:))
 end do
 nbv(:,0,:) = nbv(:,1,:)
 nbv(:,nz-1,:) = nbv(:,nz-2,:)
 nbv = sqrt(nbv)

 pv = u
 re = 6371.e3
 phi = y*3.141592/180.
 sinphi = sin(phi)
 cosphi = cos(phi)
 f = 2.*7.292e-5*sinphi 
 fa = u
 tmp = u
 do j=0, ny-1
   tmp(:,:,j) = u(:,:,j)*cosphi(j)
 end do
 tmp2 = tmp
 do j=1, ny-2
   fa(:,:,j) = f(j) - (tmp(:,:,j+1) - tmp(:,:,j-1))/(sinphi(j+1) - sinphi(j-1))/re
   tmp2(:,:,j) = (pt(:,:,j+1) - pt(:,:,j-1))/(phi(j+1) - phi(j-1))/re
 end do
 tmp3 = tmp
 do k=1, nz-2
   pv(:,k,:) = -9.8*pt(:,k,:)/(p(k)*1.e2)* \
               (lnpt(:,k+1,:) - lnpt(:,k-1,:))/log(p(k+1)/p(k-1))
   tmp3(:,k,:) = 9.8/(p(k)*1.e2)* \
                 (u(:,k+1,:) - u(:,k-1,:))/log(p(k+1)/p(k-1))
 end do
 pv = (fa*pv + tmp3*tmp2)*1.e6
 pv(:,0,:) = pv(:,1,:)
 pv(:,nz-1,:) = pv(:,nz-2,:)
 pv(:,:,0) = pv(:,:,1)
 pv(:,:,ny-1) = pv(:,:,ny-2)
 
 if (cnt .eq. 1) then
   ua = u
   ta = t
   pva = pv
   nbva = nbv
 else
   ua0 = array_append_record(ua,u,0)
   ta0 = array_append_record(ta,t,0)
   pva0 = array_append_record(pva,pv,0)
   nbva0 = array_append_record(nbva,nbv,0)
   delete( ua )
   delete( ta )
   delete( pva )
   delete( nbva )
   ua = ua0
   ta = ta0
   pva = pva0
   nbva = nbva0
   delete( ua0 )
   delete( ta0 )
   delete( pva0 )
   delete( nbva0 )
 end if

 delete( u )
 delete( t )
 delete( z )
 delete( pt )
 delete( lnpt )
 delete( nbv )
 delete( pv )
 delete( fa )
 delete( tmp )
 delete( tmp2 )
 delete( tmp3 )

end do
end do

 um = dim_avg_n(ua,0)
 tm = dim_avg_n(ta,0)
 pvm = dim_avg_n(pva,0)
 nbvm = dim_avg_n(nbva,0)


;do il=0, nt-1
 il = 0

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev1 = ispan(210,305,5) ;ispan(-10,10,1)*0.25
; lbopt = "[K]" ;"[m s:S:-1:N:]"
 cnlev2 = ispan(-15,15,1)*4
 cnlev3 = ispan(0,20,1)*0.25
 cnlev4 = (/0.2,0.3,0.4,0.5,0.7,1.0,1.5,2,3,4,5,7,10,15,20,30,40,50,70,100/)
 cnlab = 1
 cnfil = ispan(2,62,3)
 cnfil(9:) = cnfil(9:) + 3
 nfil = dimsizes(cnfil)

;--------------------------------------- contour
 wi = 24
 he = 40
 xf = 20 + (/0,1,2,3,4/)*(wi+6)
 yf = 20 + (/0,1/)*(he+9)

; mstr = "U, T"
 mstr = "" ;"Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lat [deg]"
 ystr = "p [hPa]"
; if (yf(il) .ne. max(yf(:nt-1))) then
;   xstr = ""
; end if
; if (xf(il) .ne. min(xf(:nt-1))) then
;   ystr = ""
; end if
; if ( xf(il) .ne. min(xf(:nt-1)) .or. yf(il) .ne. max(yf(:nt-1)) ) then
; if ( xf(il) .ne. max(xf(:nt-1)) .or. yf(il) .ne. min(yf(:nt-1)) ) then
;   lbopt = ""
; end if

; contour1 = cnshadeplot(wks,y,p(:nz-2),tm(:nz-2,:),cnlev1,cnfil,lbopt)
 contour1 = cnlineplot(wks,y,p(:nz-2),tm(:nz-2,:),cnlev1,0)
 cnline_att(contour1,"D",0,"T",1.0,"C",2)
; cnaddline(contour1,wks,y,p(:nz-2),tm(:nz-2,:),cnlev2,cnlab,"D",0,"T",0.6,"C",2)
 cnaddline(contour1,wks,y,p(:nz-2),um(:nz-2,:),cnlev2,cnlab,"D",(/0,0,14/),"T",(/0.5,1.2,0.5/),"C",1)
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.018
 end setvalues
 set_size(contour1,xf(il),yf(0),wi,he)
 axis_range(contour1,20,80,1000,32.7,False,True)
 axis_str(contour1,mstr,xstr,ystr)
 
 cnfinalize(contour1,wks)
 
 draw(contour1)

; contour2 = cnshadeplot(wks,y,p(:nz-2),tm(:nz-2,:),cnlev3,cnfil,lbopt)
 contour2 = cnlineplot(wks,y,p(:nz-2),nbvm(:nz-2,:)*1.e2,cnlev3,0)
 cnline_att(contour2,"D",0,"T",1.0,"C",2)
; cnaddline(contour2,wks,y,p(:nz-2),nbvm(:nz-2,:)*1.e2,cnlev4,cnlab,"D",0,"T",0.6,"C",2)
 cnaddline(contour2,wks,y,p(:nz-2),pvm(:nz-2,:),cnlev4,cnlab,"D",(/0,0,14/),"T",(/0.5,1.2,0.5/),"C",1)
 setvalues contour2
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.018
 end setvalues
 set_size(contour2,xf(il),yf(1),wi,he)
 axis_range(contour2,20,80,1000,32.7,False,True)
; axis_str(contour2,mstr,xstr,ystr)
 axis_str(contour2,"",xstr,ystr)

 cnfinalize(contour2,wks)

 draw(contour2)

 delete(contour1)
 delete(contour2)

;end do

 frame(wks)

end

