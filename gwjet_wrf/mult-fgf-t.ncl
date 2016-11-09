load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9
sz_mstr = p2n(5)

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf/x6_5min"
; fdname = "domain_x6.nc"
 ftail = "__x6__08_00_00__12_00_00__5min.nc"
 vars = (/"PRS","FGF","W","T0"/)
 lev = (/250,1500,3000,5000,250/)
 lev0 = 250
 p_app = (/950,830,670,500,950/)
 tims = (/6.5,7.25/) ;fspan(4.25,8.0,16)
; tims = (/7.+3.5/24./) ;fspan(4.25,8.0,16)
 ni_shift = 0
 ymin = 30
 ymax = 70
 xmin = (/10,20,20/)
 xmax = (/50,60,60/)
 
; wks = cr_wks("fig/multilayer_z.ps",600,0)
 wks = cr_wks("fig/mult-fgf-t-2",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
; gsn_define_colormap(wks,"cmp_flux")
; gsn_define_colormap(wks,"precip3_16lev")
 colind = NhlNewColor(wks,0.,1.,0.)
 colind2 = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- main body
 nv = dimsizes(vars)
 nz = dimsizes(lev)
 nt = dimsizes(tims)

 fnamep = fdir+"/PRS_z/prt_d300km/mean_PRS__z"+sprinti("%5.5i",lev)+ftail
 fnamet = fdir+"/T_z/prt_d300km/mean_T__z"+sprinti("%5.5i",lev)+ftail
 fnamet0 = fdir+"/T_z/prt_d300km/mean_T__z"+sprinti("%5.5i",lev0)+ftail
 fnamet0p = fdir+"/T_z/prt_d300km/prt_T__z"+sprinti("%5.5i",lev0)+ftail
 fnameu = fdir+"/U_z/prt_d300km/mean_U__z"+sprinti("%5.5i",lev)+ftail
 fnamev = fdir+"/V_z/prt_d300km/mean_V__z"+sprinti("%5.5i",lev)+ftail
 fnamen1 = fdir+"/nbe_z/mean_d300km/d_nbe-maj__z"+sprinti("%5.5i",lev)+ftail
 fnamen2 = fdir+"/nbe_z/mean_d300km/d_nbe-curv__z"+sprinti("%5.5i",lev)+ftail
 fnamen3 = fdir+"/nbe_z/mean_d300km/d_nbe-nl_min__z"+sprinti("%5.5i",lev)+ftail
 fnamen4 = fdir+"/nbe_z/mean_d300km/d_nbe-ageo__z"+sprinti("%5.5i",lev)+ftail
 fnamen5 = fdir+"/nbe_z/mean_d300km/div__z"+sprinti("%5.5i",lev)+ftail
 fnamen6 = fdir+"/nbe_z/mean_d300km/d_nbe-fvor__z"+sprinti("%5.5i",lev)+ftail
 fnamef1 = fdir+"/fgf_z/fgf-t__z"+sprinti("%5.5i",lev)+ftail
 fnames = new((/2,nz,nv/),string)
 do k=0, nz-1
 do iv=0, nv-1
   fnames(0,k,iv) = fdir+"/"+vars(iv)+"_z/prt_d300km/"+ \
                    "mean_"+vars(iv)+"__z"+sprinti("%5.5i",lev(k))+ftail
   fnames(1,k,iv) = fdir+"/"+vars(iv)+"_z/prt_d300km/"+ \
                    "prt_"+vars(iv)+"__z"+sprinti("%5.5i",lev(k))+ftail
 end do
 end do
 
 fd = addfile(fnames(0,0,0),"r")
 x = fd->lon - 40
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
do k=0, nz-2 ;nz-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   if ( vars(iv) .eq. "UV" ) then
     uuz = addfile(fnameu(iv),"r")
     vvz = addfile(fnamev(iv),"r")
     vs (iv,:,:) = uuz->mean_U(it(n),:,:)
     vs2 = vvz->mean_V(it(n),:,:)
     vs(iv,:,:) = (/ sqrt(vs(iv,:,:)^2 + vs2(:,:)^2) /)
   end if
   if (vars(iv) .eq. "W") then
     fp = addfile(fnames(1,k,iv),"r")
     vs(iv,:,:) = fp->prt_W(it(n),:,:) * 1.e3
     vs(iv,:,:) = vs(iv,:,:)*exp(-lev(k)/(2.*7.e3))
   end if
   if (vars(iv) .eq. "V") then
     fp = addfile(fnames(1,k,iv),"r")
     vs(iv,:,:) = fp->prt_V(it(n),:,:) * 10.
     vs(iv,:,:) = vs(iv,:,:)*exp(-lev(k)/(2.*7.e3))
   end if
   if (vars(iv) .eq. "PRS") then
     fm = addfile(fnames(0,k,iv),"r")
     vs(iv,:,:) = fm->mean_P(it(n),:,:) / 100.
   end if
   if (vars(iv) .eq. "T") then
     fm = addfile(fnames(0,k,iv),"r")
     vs(iv,:,:) = fm->mean_T(it(n),:,:) + 300.
   end if
   if (vars(iv) .eq. "T0") then
     fm = addfile(fnamet0,"r")
     fp = addfile(fnamet0p,"r")
     vs(iv,:,:) = fm->mean_T(it(n),:,:) + 300.
     vs(iv,:,:) = vs(iv,:,:) + fp->prt_T(it(n),:,:)
   end if
   if (vars(iv) .eq. "D_NBE") then
     vs(iv,:,:) = 0.
     nnz = addfile(fnamen1(k),"r")
     vs(iv,:,:) = vs(iv,:,:) + nnz->d_nbe_maj(it(n),:,:)
;     nnz = addfile(fnamen2(k),"r")
;     vs(iv,:,:) = vs(iv,:,:) + nnz->curv(it(n),:,:)
;     nnz = addfile(fnamen3(k),"r")
;     vs(iv,:,:) = vs(iv,:,:) + nnz->nl_min(it(n),:,:)
;     nnz = addfile(fnamen4(k),"r")
;     vs(iv,:,:) = vs(iv,:,:) - nnz->ageo(it(n),:,:)
     vs(iv,:,:) = vs(iv,:,:)*1.e9
   end if
   if (vars(iv) .eq. "FGF") then
     ffz = addfile(fnamef1(k),"r")
     vs(iv,:,:) = ffz->fgf(it(n),:,:)
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
 ystr = ""
 tday = toint(tims(n)+1.e-4)
 tmp = (tims(n)+1.e-4 - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)
 txtstr = "Day "+tday+",  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 if (n .eq. 0) then
   mstr = "z = "+(lev(k)/1.e3)+" km"
   if (lev(k) .lt. 1000) then
     mstr = "z = "+lev(k)+" m"
   end if
 end if
 if (n .eq. nt-1) then
   xstr = "longitude [deg]"
 end if
 if (k .eq. 0) then
   ystr = "latitude [deg]"
 end if

 scl = (/1.,1.,1.,1./)
 cnlab = (/0,0,0,0/)
 lbopt = (/"","","",""/)

 cnlev3 = ispan(-20,20,1)*8 + 300
; cnlev1 = ispan(1,10,1)*5.
 cnlev1 = (/1,10,100,1000/)
; cnlev2 = ispan(-10,10,1)*1.5
 cnlev2 = ispan(-5,5,1)*4.
; cntmp1 = 2^(ispan(1,10,1)*0.5)
; cntmp0 = 2^(ispan(10,1,1)*0.5)*(-1.)
; cnlev2(:9 ) = cntmp0
; cnlev2(11:) = cntmp1
 if (lev(k) .lt. 5.e3) then
   cnlev0 = p_app(k) + ispan(-20,20,1)*8.
 else
   cnlev0 = p_app(k) + ispan(-20,20,1)*8.
 end if
 if ( n .eq. nt-1 .and. k .eq. nz-4 ) then
; if ( k .eq. nz-4 ) then
   lbopt = (/"","","[mm s:S:-1:N:]","[K]"/)
 end if
; cnfil = ispan(2,65,3)
; nfil = dimsizes(cnfil)
; cnfil(nfil/2-2:nfil/2+1) = 0
 cnfil = array_append_record( ispan(2,32,6), ispan(35,65,6), 0 )
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0

; plot -----------------------------------------------------
 pscl = 0.6 ;0.45
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xmax(0) - xmin(0))*rlonlat*pscl
 he = (ymax - ymin)*sqrt( 1. - (ni_shift*rlonlat*dx/dy)^2 )*pscl
 xf = 10 + (/0,1,2,3,4/)*(wi+5) + (/0,0,0,0,1/)*3
 yf = 20 + (/0,1,2,3,4/)*(he+5)

if (k .ne. nz-1) then
 contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:)/scl(2),cnlev2,cnfil,lbopt(2))
end if
if (k .eq. nz-1) then
; cnaddline(contour1,wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,cnlab(0), \
;           "D",0,"T",1.0,"C",colind2)
 contour1 = cnlineplot(wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,0)
 cnline_att(contour1,"D",0,"T",1.0,"C",11)
end if
 cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnlab(1), \
           "D",0,"T",0.4,"C",(/1,-1,1/))
; cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1*0.1,cnlab(1), \
;           "D",0,"T",0.5,"C",(/1,-1,1/)*colind2)
if (k .ne. nz-1) then
 cnaddline(contour1,wks,xpl,y,vpl(3,:,:)/scl(3),cnlev3,cnlab(3), \
           "D",0,"T",0.7,"C",colind)
end if
 setvalues contour1
  "tmXBLabelsOn" : False
;  "tmYLLabelsOn" : False
  "pmLabelBarOrthogonalPosF" : 0.10
;  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.1
  "lbTitlePosition"  : "Right"
  "lbTitleDirection" : "Across"
  "lbLabelFontHeightF" : 0.04
  "lbTitleOffsetF"     : 0.04
 end setvalues
; if (k .eq. 0) then
   setvalues contour1
    "tmXBLabelsOn" : True
   end setvalues
; end if
; if (k .eq. 0) then
;   setvalues contour1
;    "tmYLLabelsOn" : True
;   end setvalues
; end if
 set_size(contour1,xf(k),yf(n),wi,he)
 axis_range(contour1,xmin(n),xmax(n),ymin,ymax,False,False)
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

 delete(contour1)

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop
 dr_txt(wks,txtstr,xf(0)+wi/2-2,yf(n)+he-1.2,22,4.,0.)
end do  ; global loop

; dr_txt(wks,"(a)",xf(0   )-3.0,yf(0)-3.5,25,6.,0.)
; dr_txt(wks,"(b)",xf(nz-1)-3.0,yf(0)-3.5,25,6.,0.)

frame(wks)

end

