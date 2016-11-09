load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "/data18/GW_jet_wrf/dat"
; n1s = toint((fspan(6.5,8,4)-4)*4)
 n1s = toint((fspan(4,8,5)-4)*4)

 wks = cr_wks("fig/xmean",600,0)
; gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/U/U__x1__08_00_00__12_00_00__6hrly.nc"
 fname2 = fdir+"/x1/T/T__x1__08_00_00__12_00_00__6hrly.nc"
 fname3 = fdir+"/x1/PRS/PRS__x1__08_00_00__12_00_00__6hrly.nc"
 fname4 = fdir+"/x1/GEOP/GEOP__x1__08_00_00__12_00_00__6hrly.nc"
 fname0 = fdir+"/x1/domain_x1.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f3 = addfile(fname3,"r")
 f4 = addfile(fname4,"r")
 f0 = addfile(fname0,"r")
 nx = 432
 ny = 648+2
 nz = 44 + 1
 iy0 = 53
 iy9 = iy0 + ny - 1
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

 nt = dimsizes(n1s)

 y = ys + ispan(0,ny-1,1)*dy

 pa = f0->ZNU(:nz-1)*990. + 10.

 u = 0.5*( f1->U(n1s,:,iy0:iy9,0:nx-1) + \
           f1->U(n1s,:,iy0:iy9,1:nx  ) )
 t = f2->T(n1s,:,iy0:iy9,:) + 300.
 lnpt = log( t )
 p = f3->P(n1s,:,iy0:iy9,:)/1.e2
 t = t*(p/1.e3)^0.286
 z = 0.5*( f4->PH(n1s,0:48,iy0:iy9,:) + \
           f4->PH(n1s,1:49,iy0:iy9,:) )/9.8

 u_p = linint1_n(log(p),u,False,log(pa),0,1)
 t_p = linint1_n(log(p),t,False,log(pa),0,1)
 z_p = linint1_n(log(p),z,False,log(pa),0,1)

 u_yzt = dim_avg(u_p)
 t_yzt = dim_avg(t_p)
 z_yzt = dim_avg(z_p)

 nbv = t_yzt
 pt_yzt = t_yzt
 do k=0, nz-1
   pt_yzt(:,k,:) = t_yzt(:,k,:)*(1.e3/pa(k))^0.286
 end do
 lnpt_yzt = log( pt_yzt )
 do k=1, nz-2
   nbv(:,k,:) = 9.8*(lnpt_yzt(:,k+1,:) - lnpt_yzt(:,k-1,:)) / \
                    (z_yzt   (:,k+1,:) - z_yzt   (:,k-1,:))
 end do
 nbv(:,0,:) = nbv(:,1,:)
 nbv(:,nz-1,:) = nbv(:,nz-2,:)
 nbv = sqrt(nbv)

 pv = u_yzt
 re = 6371.e3
 phi = y*3.141592/180.
 sinphi = sin(phi)
 cosphi = cos(phi)
 f = 2.*7.292e-5*sinphi 
 fa = u_yzt
 tmp = u_yzt
 do j=0, ny-1
   tmp(:,:,j) = u_yzt(:,:,j)*cosphi(j)
 end do
 tmp2 = tmp
 do j=1, ny-2
   fa(:,:,j) = f(j) - (tmp(:,:,j+1) - tmp(:,:,j-1))/(sinphi(j+1) - sinphi(j-1))/re
   tmp2(:,:,j) = (pt_yzt(:,:,j+1) - pt_yzt(:,:,j-1))/(phi(j+1) - phi(j-1))/re
 end do
 tmp3 = tmp
 do k=1, nz-2
   pv(:,k,:) = -9.8*pt_yzt(:,k,:)/(pa(k)*1.e2)* \
               (lnpt_yzt(:,k+1,:) - lnpt_yzt(:,k-1,:))/log(pa(k+1)/pa(k-1))
   tmp3(:,k,:) = 9.8/(pa(k)*1.e2)* \
                 (u_yzt(:,k+1,:) - u_yzt(:,k-1,:))/log(pa(k+1)/pa(k-1))
 end do
 pv = (fa*pv + tmp3*tmp2)*1.e6
 pv(:,0,:) = pv(:,1,:)
 pv(:,nz-1,:) = pv(:,nz-2,:)
 pv(:,:,0) = pv(:,:,1)
 pv(:,:,ny-1) = pv(:,:,ny-2)
 

do il=0, nt-1

;--------------------------------------- process 0
 tday  = n1s(il)
 thour = mod(tday,4)*6
 tday  = tday/4
 tmin  = 0
 tday  = 4 + tday

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
 mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "lat [deg]"
 ystr = "p [hPa]"
; if (yf(il) .ne. max(yf(:nt-1))) then
;   xstr = ""
; end if
 if (xf(il) .ne. min(xf(:nt-1))) then
   ystr = ""
 end if
; if ( xf(il) .ne. min(xf(:nt-1)) .or. yf(il) .ne. max(yf(:nt-1)) ) then
; if ( xf(il) .ne. max(xf(:nt-1)) .or. yf(il) .ne. min(yf(:nt-1)) ) then
;   lbopt = ""
; end if

; contour1 = cnshadeplot(wks,y,pa(:nz-2),t_yzt(il,:nz-2,:),cnlev1,cnfil,lbopt)
 contour1 = cnlineplot(wks,y,pa(:nz-2),t_yzt(il,:nz-2,:),cnlev1,0)
 cnline_att(contour1,"D",0,"T",1.0,"C",2)
; cnaddline(contour1,wks,y,pa(:nz-2),t_yzt(il,:nz-2,:),cnlev2,cnlab,"D",0,"T",0.6,"C",2)
 cnaddline(contour1,wks,y,pa(:nz-2),u_yzt(il,:nz-2,:),cnlev2,cnlab,"D",(/0,0,14/),"T",(/0.5,1.2,0.5/),"C",1)
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
 axis_range(contour1,-999,999,1000,-999,False,True)
 axis_str(contour1,mstr,xstr,ystr)
 
 cnfinalize(contour1,wks)
 
 draw(contour1)

; contour2 = cnshadeplot(wks,y,pa(:nz-2),t_yzt(il,:nz-2,:),cnlev3,cnfil,lbopt)
 contour2 = cnlineplot(wks,y,pa(:nz-2),nbv(il,:nz-2,:)*1.e2,cnlev3,0)
 cnline_att(contour2,"D",0,"T",1.0,"C",2)
; cnaddline(contour2,wks,y,pa(:nz-2),nbv(il,:nz-2,:)*1.e2,cnlev4,cnlab,"D",0,"T",0.6,"C",2)
 cnaddline(contour2,wks,y,pa(:nz-2),pv(il,:nz-2,:),cnlev4,cnlab,"D",(/0,0,14/),"T",(/0.5,1.2,0.5/),"C",1)
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
 axis_range(contour2,-999,999,1000,-999,False,True)
; axis_str(contour2,mstr,xstr,ystr)
 axis_str(contour2,"",xstr,ystr)

 cnfinalize(contour2,wks)

 draw(contour2)

 delete(contour1)
 delete(contour2)

end do

 frame(wks)

end

