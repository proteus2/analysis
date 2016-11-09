load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 wks = cr_wks("fig/init_yz",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,1.0,0.0,0.0)

;--------------------------------------- read data
 fdir = "/data18/kyh/dat/gwjet_wrf/x1/init"
 fdir2 = "/data18/kyh/dat/gwjet_wrf-2/x1/init"
 fname1 = fdir+"/wrfinit.nc"
 f1 = addfile(fname1,"r")

 z00= dim_avg( f1->PH(0,:,:,:) + f1->PHB(0,:,:,:) )/9.81 * 1.e-3
 lnpt0 = log( dim_avg( f1->T(0,:,:,:) ) + 300. )
 u0 = dim_avg( f1->U(0,:,:,:) )
 p0 = dim_avg( f1->P(0,:,:,:) + f1->PB(0,:,:,:) )

 y = f1->XLAT(0,:,0)
 ny = dimsizes(y)
 nz0 = 49

 z00  (1:nz0-2,:) = 0.25*(z00  (0:nz0-3,:)+2.*z00  (1:nz0-2,:)+z00  (2:nz0-1,:))
 u0   (1:nz0-2,:) = 0.25*(u0   (0:nz0-3,:)+2.*u0   (1:nz0-2,:)+u0   (2:nz0-1,:))
 lnpt0(1:nz0-2,:) = 0.25*(lnpt0(0:nz0-3,:)+2.*lnpt0(1:nz0-2,:)+lnpt0(2:nz0-1,:))

 lnrho0 = log( p0/287./exp(lnpt0)*(1.e5/p0)^0.286 )

;--------------------------------------- process 0
 z00 = exp(z00/(-7.2))
 z0 = log(0.5*(z00(0:nz0-1,:) + z00(1:nz0,:)))*(-7.2)
 
 nz1 = nz0-1

 z0h = 0.5*(z0(0:nz1-1,:) + z0(1:nz1,:))
 u0_z = (u0(1:nz1,:) - u0(0:nz1-1,:))/(z0(1:nz1,:) - z0(0:nz1-1,:))*1.e-3
 n2bv0 = 9.81*(lnpt0(1:nz1,:) - lnpt0(0:nz1-1,:))/(z0(1:nz1,:) - z0(0:nz1-1,:))*1.e-3
; term3 = (-1.)*(2.*7.292e-5*sin(45/180.*3.141592))^2/n2bv0*u0_z
 term3 = (-1.)*conform_dims((/nz0-1,ny/),(2.*7.292e-5*sin(y/180.*3.141592))^2,1)/n2bv0*u0_z
 term4 = (lnrho0(1:nz1,:) - lnrho0(0:nz1-1,:))/(z0(1:nz1,:) - z0(0:nz1-1,:))*1.e-3*term3

 nz1 = nz0-2

 tmp = 1./n2bv0
 z0hh = 0.5*(z0h(0:nz1-1,:) + z0h(1:nz1,:))
 u0_zz = (u0_z(1:nz1,:) - u0_z(0:nz1-1,:))/(z0h(1:nz1,:) - z0h(0:nz1-1,:))*1.e-3
 n2bv0_i_z = (tmp(1:nz1,:) - tmp(0:nz1-1,:))/(z0h(1:nz1,:) - z0h(0:nz1-1,:))*1.e-3
 delete(tmp)
 term3_z = (term3(1:nz1,:) - term3(0:nz1-1,:))/(z0h(1:nz1,:) - z0h(0:nz1-1,:))*1.e-3


 var0   = (/lnpt0,u0/)
 var0h  = (/u0_z,n2bv0,term4/)
 var0hh = (/u0_zz,n2bv0_i_z,term3_z/)

 z = ispan(-250,25000,250)*1.e-3
 nz = dimsizes(z)

 nd0   = dimsizes(var0  )
 nd0h  = dimsizes(var0h )
 nd0hh = dimsizes(var0hh)
 nv   = nd0  (0)
 nvh  = nd0h (0)
 nvhh = nd0hh(0)

 varc   = new((/nv  ,ny,nz0  /),float)
 varhc  = new((/nvh ,ny,nz0-1/),float)
 varhhc = new((/nvhh,ny,nz0-2/),float)
 z0c   = varc
 z0hc  = varhc
 z0hhc = varhhc
 do j=0, ny-1
   varc  (:,j,:) = var0  (:,:,j)
   varhc (:,j,:) = var0h (:,:,j)
   varhhc(:,j,:) = var0hh(:,:,j)
   z0c  (0,j,:) = z0  (:,j)
   z0hc (0,j,:) = z0h (:,j)
   z0hhc(0,j,:) = z0hh(:,j)
 end do
 do iv=1, nv-1
   z0c(iv,:,:) = z0c(0,:,:)
 end do
 do iv=1, nvh-1
   z0hc(iv,:,:) = z0hc(0,:,:)
 end do
 do iv=1, nvhh-1
   z0hhc(iv,:,:) = z0hhc(0,:,:)
 end do

 var1   = ftcurv(z0c  ,varc  ,z)
 var1h  = ftcurv(z0hc ,varhc ,z)
 var1hh = ftcurv(z0hhc,varhhc,z)

 var   = new((/nv,nz,ny/),float)
 varh  = new((/nvh,nz,ny/),float)
 varhh = new((/nvhh,nz,ny/),float)
 do j=0, ny-1
   var  (:,:,j) = var1  (:,j,:)
   varh (:,:,j) = var1h (:,j,:)
   varhh(:,:,j) = var1hh(:,j,:)
 end do

 var(0,:,:) = exp(var(0,:,:))

 ; baroclinic insta.
 beta = conform_dims((/nz,ny/),2.*7.292e-5/6370.e3*cos(y/180.*3.141592),1)
 u_yy = new((/nz,ny/),float)
 inv_dy2 = 1./((y(1) - y(0))/180.*3.141592*6370.e3)^2
 u_yy(:,1:ny-2) = (var(1,:,0:ny-3) + var(1,:,2:ny-1) - 2.*var(1,:,1:ny-2))*inv_dy2
; f02 = (2.*7.292e-5*sin(45./180.*3.141592))^2
 f02 = conform_dims((/nz,ny/),(2.*7.292e-5*sin(y/180.*3.141592))^2,1)
 n2bv = varh(1,:,:)
 u_z  = varh(0,:,:)

 binstc = varhh(2,:,:) + varh(2,:,:)

 binst = beta - u_yy + binstc

 binst@_FillValue = 1.e20
 binst(:,0   ) = 1.e20
 binst(:,ny-1) = 1.e20
 binst(0   ,:) = 1.e20
 binst(nz-1,:) = 1.e20

;--------------------------------------- cn level
 cnlev1 = ispan(200,800,10)
 cnlev2 = ispan(0,15,1)*5
 cnlev3 = ispan(-4,14,1)*1.0
 cnfil = ispan(2,26,1) ;(/22,23,24,25,26,0,11,10,9,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 50
 he = 30
 xf = 20
 yf = 20

 xstr = "latitude [deg]"
 ystr = "z [km]"

; contour1 = cnshadeplot(wks,y,z,binst*1.e11,cnlev3,cnfil,"[10:S:-11:N:]")
; cnaddline(contour1,wks,y,z,var(0,:,:),cnlev1,0, \
;           "D",0,"T",0.6,"C",2)
 contour1 = cnlineplot(wks,y,z,var(0,:,:),cnlev1,0)
 cnline_att(contour1,"D",0,"T",0.8,"C",colind)

 cnaddline(contour1,wks,y,z,var(1,:,:),cnlev2,1, \
           "D",0,"T",1.0,"C",1)
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
 set_size(contour1,xf,yf,wi,he)
 axis_range(contour1,0,90,0,15,False,False)
 axis_str(contour1,"",xstr,ystr)
 setvalues contour1
  "tmXBMode"          : "Manual"
  "tmXBTickSpacingF"  : 15
  "tmXBMinorPerMajor" : 2
 end setvalues
 
 cnfinalize(contour1,wks)
 
 draw(contour1)

 frame(wks)

end

