load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

ta = fspan(4,8,2*96+1)  ; 2 times per hour
do nnn=0, dimsizes(ta)-1
;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 ds = 300
 tims = ta(nnn)
 n1s = toint(288*(tims-4))
 ipt0 = (/216,367/) + (/0,0/)
 iitv = (/1,0/)
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 tims = round(tims*1.e4,0)*1.e-4
 ctims = tims+""
 if (tims(0) .eq. tofloat(toint(tims(0)))) then
   ctims = tims+".0"
 end if
 wks = cr_wks("fig/prt"+ds+"_hz54-"+var1+"-"+ctims,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nx = 432
 ny = 648+2
 nz = 101
 dz = 250
 dx = 1./10.8
 dy = 1./10.8
 xs = 20. + 0.5*dx  ; yn = 80. + 0.5*dy
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

nt = dimsizes(n1s)

do il=0, nt-1

 tday = toint(tims(il)+1.e-4)
 tmp = (tims(il)+1.e-4 - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)

 v1 = new((/nz,ny,nx/),float)
 v2 = new((/nz,ny,nx/),float)
 v3 = new((/nz,ny,nx/),float)
 sq_rho = new((/nz,ny,nx/),float)

 z = ispan(0,(nz-1)*dz,dz)
 do k=0, nz-1
   fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d"+ds+"km/"+var1name+"__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fname2 = fdir+"/x6_5min/U_z/prt_d"+ds+"km/mean_U__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fname3 = fdir+"/x6_5min/V_z/prt_d"+ds+"km/mean_V__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   f1 = addfile(fname1,"r")
   f2 = addfile(fname2,"r")
   f3 = addfile(fname3,"r")

;   v1(k,:,:) = f1->$var1name$(n1s(il),:,:)
;   v2(k,:,:) = f2->mean_U(n1s(il),:,:)
;   v3(k,:,:) = f3->mean_V(n1s(il),:,:)
   v1(k,0,:) = f1->$var1name$(n1s(il),ipt0(1),:)
   v2(k,0,:) = f2->mean_U(n1s(il),ipt0(1),:)
   v3(k,0,:) = f3->mean_V(n1s(il),ipt0(1),:)

   fnamep = fdir+"/x6_5min/PRS_z/prt_d"+ds+"km/mean_PRS__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   fnamet = fdir+"/x6_5min/T_z/prt_d"+ds+"km/mean_T__z"+sprinti("%5.5i",z(k))+"__x6__08_00_00__12_00_00__5min.nc"
   ppz = addfile(fnamep,"r")
   ttz = addfile(fnamet,"r")
;   vp  = ppz->mean_P(n1s(il),:,:)
;   vpt = ttz->mean_T(n1s(il),:,:) + 300.
;   sq_rho(k,:,:) = sqrt( vp/287./(vpt*(vp/1.e5)^0.286) )
   vp  = ppz->mean_P(n1s(il),ipt0(1),:)
   vpt = ttz->mean_T(n1s(il),ipt0(1),:) + 300.
   sq_rho(k,0,:) = sqrt( vp/287./(vpt*(vp/1.e5)^0.286) )
 end do

;--------------------------------------- process 0
 x = xs + ispan(0,nx-1,1)*dx - 40
 y = ys + ispan(0,ny-1,1)*dy

 nxy = 430 ;101
 xx = x(ipt0(0)-nxy/2*iitv(0):ipt0(0)+nxy/2*iitv(0):abs(iitv(0)))
; yy = y(ipt0(1)-nxy/2*iitv(1):ipt0(1)+nxy/2*iitv(1):abs(iitv(1)))
 yy = y(ipt0(1))

 v1i = new((/nz,nxy/),float)
 v2i = new((/nz,nxy/),float)
 v3i = new((/nz,nxy/),float)
 sqr = new((/nz,nxy/),float)
 do ij=0, nxy-1
;   v1i(:,ij) = v1(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
;   v2i(:,ij) = v2(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
;   v3i(:,ij) = v3(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
;   sqr(:,ij) = sq_rho(:,ipt0(1)+(ij-nxy/2)*iitv(1),ipt0(0)+(ij-nxy/2)*iitv(0))
   v1i(:,ij) = v1(:,0,ipt0(0)+(ij-nxy/2)*iitv(0))
   v2i(:,ij) = v2(:,0,ipt0(0)+(ij-nxy/2)*iitv(0))
   v3i(:,ij) = v3(:,0,ipt0(0)+(ij-nxy/2)*iitv(0))
   sqr(:,ij) = sq_rho(:,0,ipt0(0)+(ij-nxy/2)*iitv(0))
 end do
 v1i = v1i*sqr
; v23i = (v2i*iitv(0) + v3i*iitv(1))/sqrt(iitv(0)^2+iitv(1)^2)
 
;--------------------------------------- cn level
 if (var1 .eq. "W") then
   scl = 1.e-3
   cnlev1 = ispan(-10,10,1)*1.
   lbopt = "[mm s:S:-1:N:]"
 end if
 if (var1 .eq. "U") then
   scl = 1.
   cnlev1 = ispan(-10,10,1)*0.25
   lbopt = "[m s:S:-1:N:]"
 end if
 if (var1 .eq. "V") then
   scl = 1.
   cnlev1 = ispan(-10,10,1)*0.25
   lbopt = "[m s:S:-1:N:]"
 end if
 cnlev2 = ispan(-10,10,1)*4
 cnlab = 1
 cnfil = ispan(2,65,3)
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 20 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 20 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+9)

 mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr = "longitude [deg]"
 xstr2 = "latitude [deg]"
 ystr = "z [km]"
 if (yf(il) .ne. max(yf(:nt-1))) then
   xstr = ""
 end if
 if (xf(il) .ne. min(xf(:nt-1))) then
   ystr = ""
 end if
; if ( xf(il) .ne. min(xf(:nt-1)) .or. yf(il) .ne. max(yf(:nt-1)) ) then
 if ( xf(il) .ne. max(xf(:nt-1)) .or. yf(il) .ne. min(yf(:nt-1)) ) then
   lbopt = ""
 end if

 contour1 = cnshadeplot(wks,xx,z/1.e3,v1i/scl,cnlev1,cnfil,lbopt)
; cnaddline(contour1,wks,xx,z/1.e3,v23i,cnlev2,cnlab,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.25
;  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues
; labt = ispan(toint(yy(0))/2*2,toint(yy(nxy-1)),2)
; minort = fspan(min(labt)-2,max(labt)+2,(max(labt)-min(labt)+4)*2+1)
; setvalues contour1
;  "tmXUseBottom" : False
;  "tmXTLabelsOn" : True
;  "tmXTMode"   : "Explicit"
;  "tmXTValues" : linint1(yy,xx,False,labt,0)
;  "tmXTLabels" : labt
;  "tmXTMinorValues" : linint1(yy,xx,False,minort,0)
; end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-20,20,-999,999,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)
 
 draw(contour1)

 delete(contour1)

; dr_txt(wks,xstr2,xf(0)+wi/2.,yf(0)-6.5,22,6.0,0.)

end do

 frame(wks)

end do  ; ta

end

