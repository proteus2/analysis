load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf/x6_5min"
; fdname = "domain_x6.nc"
 ftail = "__x6__08_00_00__12_00_00__5min.nc"
 vars = (/"PRS","T","U"/)
; tims = fspan(4,8,96+1) ;(/6.5/)  ; [day]
; tims = fspan(4,8,9) ;(/6.5/)  ; [day]
 tims = ispan(4,8,1)*1.
 levs = (/500,1500,3000,5000,8000/)
 p_app = (/930,830,670,500,350/)
 ni_shift = 1
 ymin = 25
 ymax = 75
 
; wks = cr_wks("fig/multilayer_z.ps",600,0)
 wks = cr_wks("fig/mult-b",600,0)
; gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"cmp_b2r")
; gsn_define_colormap(wks,"cmp_flux")
 gsn_define_colormap(wks,"precip3_16lev")
 colind = NhlNewColor(wks,1.,0.,0.)

;--------------------------------------- main body
 nv = dimsizes(vars)
 nt = dimsizes(tims)
 nz = dimsizes(levs)

 fnamep = fdir+"/PRS_z/prt_d300km/mean_PRS__z"+sprinti("%5.5i",levs)+ftail
 fnamet = fdir+"/T_z/prt_d300km/mean_T__z"+sprinti("%5.5i",levs)+ftail
 fnameu = fdir+"/U_z/prt_d300km/mean_U__z"+sprinti("%5.5i",levs)+ftail
 fnamev = fdir+"/V_z/prt_d300km/mean_V__z"+sprinti("%5.5i",levs)+ftail
 fnames = new((/2,nv,nz/),string)
 do iv=0, nv-1
   fnames(0,iv,:) = fdir+"/"+vars(iv)+"_z/prt_d300km/"+ \
                    "mean_"+vars(iv)+"__z"+sprinti("%5.5i",levs)+ftail
   fnames(1,iv,:) = fdir+"/"+vars(iv)+"_z/prt_d300km/"+ \
                    "prt_"+vars(iv)+"__z"+sprinti("%5.5i",levs)+ftail
 end do
 
 fd = addfile(fnames(0,0,0),"r")
 x = fd->lon
 y = fd->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 dx = x(1) - x(0)
 dy = y(1) - y(0)
 nt0 = dimsizes( fd->t )

 tims = round(tims*1.e4,0)*1.e-4
 time = round(fspan(4,8,nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(time,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 nxp = nx*2 + (ny-1)*ni_shift

do n=0, nt-1  ; global loop
do k=0, nz-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   fm = addfile(fnames(0,iv,k),"r")
   fp = addfile(fnames(1,iv,k),"r")
   if ( vars(iv) .eq. "U" .or. vars(iv) .eq. "V" ) then
     uuz = addfile(fnameu(k),"r")
     vvz = addfile(fnamev(k),"r")
     vs (iv,:,:) = uuz->mean_U(it(n),:,:)
     vs2 = vvz->mean_V(it(n),:,:)
     vs(iv,:,:) = (/ sqrt(vs(iv,:,:)^2 + vs2(:,:)^2) /)
   end if
   if (vars(iv) .eq. "W") then
     ppz = addfile(fnamep(k),"r")
     ttz = addfile(fnamet(k),"r")
     vs2 = ppz->mean_P(it(n),:,:)
     vs3 = ttz->mean_T(it(n),:,:) + 300.
     sq_rho = sqrt( vs2/287./(vs3*(vs2/1.e5)^0.286) )
     vs(iv,:,:) = fp->prt_W(it(n),:,:) * (sq_rho*1.e3)
   end if 
   if (vars(iv) .eq. "PRS") then
     vs(iv,:,:) = fm->mean_P(it(n),:,:) / 100.
   end if
   if (vars(iv) .eq. "T") then
     vs(iv,:,:) = fm->mean_T(it(n),:,:) + 300.
   end if
 end do

; process --------------------------------------------------
 xpl = new(nxp,float,1.e20)
 xpl(:nx-1) = x(:)
 xpl(nx:) = xpl(nx-1) + (x(1)-x(0))*ispan(1,nxp-nx,1)

 vpl = new((/nv,ny,nxp/),float,1.e20)
 do j=0, ny-1
   vpl(:,j,j*ni_shift:j*ni_shift+nx-1) = vs(:,j,:)
   vpl(:,j,j*ni_shift+nx:j*ni_shift+nx+nx-1) = vs(:,j,:)
 end do

; plot setting ---------------------------------------------
 mstr = ""
 xstr = ""
 if (k .eq. nz-1) then
   tday = toint(tims(n)+1.e-4)
   tmp = (tims(n)+1.e-4 - tday)*24.
   thour = toint(tmp)
   tmin = toint((tmp - thour)*60.)
   mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 end if
; if (k .eq. 0) then
;   xstr = "lon [deg]"
; end if
 ystr = "" ;"lat [deg]"

 scl = (/1.,1.,1./)
 cnlab = (/1,1,1/)
 lbopt = (/"","[m s:S:-1:N:]","[mm s:S:-1:N:]"/)

 cnlev1 = ispan(-10,10,1)*10 + 270.
 cnlev2 = ispan(0,15,1)*3.
 if (levs(k) .lt. 5.e3) then
   cnlev0 = p_app(k) + ispan(-20,20,1)*8.;60.
 else
   cnlev0 = p_app(k) + ispan(-20,20,1)*8.
 end if
 if ( k .ne. nz-1 .or. n .ne. nt-1 ) then
   lbopt = ""
 end if
 cnfil = ispan(2,18,1)
 nfil = dimsizes(cnfil)
; cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; plot -----------------------------------------------------
 pscl = 0.4
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xpl(nxp-1) - xpl(0))*rlonlat*pscl
 he = (ymax - ymin)*sqrt( 1. - (ni_shift*rlonlat*dx/dy)^2 )*pscl
 he = he*(ymax - ymin)/(y(ny-1) - y(0))
 xf = 10 + (/0,1,2,3,4,5/)*(wi-10)
 yf = 20 + ispan(nz-1,0,1)*(he-0.3)

 contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:)/scl(2),cnlev2,cnfil,lbopt(2))
 cnaddline(contour1,wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,cnlab(0), \
           "D",0,"T",0.6,"C",1)
 if (k .eq. 0) then
   cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnlab(1), \
             "D",0,"T",0.6,"C",colind)
 end if
 setvalues contour1
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tmYRBorderOn" : False
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
  "lbOrientation" : "Vertical"
  "pmLabelBarOrthogonalPosF" : -1.20
  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.1
  "pmLabelBarHeightF"        : 1.1
  "lbLabelFontHeightF" : 0.012
  "lbTitleOffsetF"     : 0.006
  "lbTitleFontHeightF" : 0.015
 end setvalues
 set_size(contour1,xf(n),yf(k),wi,he)
 axis_range(contour1,-999,999,ymin,ymax,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 xpl1 = x(0)+(ymin-y(0))/dy*dx*ni_shift
 xpl4 = x(0)+(ymax-y(0))/dy*dx*ni_shift
 xpl2 = xpl1 + 2.*(x(nx-1) - x(0))
 xpl3 = xpl4 + 2.*(x(nx-1) - x(0))
 boxres = True
 boxres@gsLineThicknessF = 1.5
 dumm = gsn_add_polyline(wks,contour1,(/xpl1,xpl2,xpl3,xpl4,xpl1/), \
                                      (/ymin,ymin,ymax,ymax,ymin/),boxres)

 draw(contour1)

 dr_txt(wks,(toint(levs(k)/1.e2)/10.)+" km",xf(0)-4,yf(k)+he-1.1,22,5.5,0.)

 delete(contour1)

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop
;frame(wks)
end do  ; global loop
frame(wks)

end

