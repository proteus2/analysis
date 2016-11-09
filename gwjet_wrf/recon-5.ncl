load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8
;sz_mstr = p2n(5)

 xy0p = (/ (/42.1,50/), \
           (/29,59/), \
           (/54,59/), \
           (/30,41.0/), \
           (/28.9,45.1/) /)
 iitv0 = (/ (/1,0/), \
            (/1,0 /), \
            (/1,1 /), \
            (/1,0 /), \
            (/1,-1/) /)
 xyl0 = (/6,6,6,24,4.2/)

begin

; wks = cr_wks("fig/multilayer_z.ps",600,0)
 wks = cr_wks("fig/recon-5-8",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
; gsn_define_colormap(wks,"cmp_flux")
; gsn_define_colormap(wks,"precip3_16lev")
 colind = NhlNewColor(wks,0.,1.,0.)

do iw=0, 4  ; global loop
;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf/x6_5min"
; fdname = "domain_x6.nc"
 if (iw .lt. 3) then
   ftail = "__x6__10_00_00__11_12_00__1hrly-w"+(iw+1)+".nc"
   t_09 = (/6,7.5/)
   tims = 7.0
   ymin = 48
   ymax = 70
 else
   ftail = "__x6__11_00_00__12_00_00__1hrly-w"+(iw+1)+".nc"
   t_09 = (/7,8/)
   tims = 7.5
   ymin = 30
   ymax = 46
 end if
 fmtail = "__x6__08_00_00__12_00_00__5min.nc"
 vars = (/"PRS","T","W"/)
 lev = (/250,250,8000/)
 p_app = (/950,950,350/)
; vars = (/"W","W","W"/)
; lev = (/8000,8000,8000/)
; p_app = (/350,350,350/)
 ni_shift = 0
 xmin = 20
 xmax = 60
 
;--------------------------------------- main body
 nv = dimsizes(vars)
 nt = dimsizes(tims)

 fnamep = fdir+"/PRS_z/prt_d300km/mean_PRS__z"+sprinti("%5.5i",lev)+fmtail
 fnamet = fdir+"/T_z/prt_d300km/mean_T__z"+sprinti("%5.5i",lev)+fmtail
 fnameu = fdir+"/U_z/prt_d300km/mean_U__z"+sprinti("%5.5i",lev)+fmtail
 fnamev = fdir+"/V_z/prt_d300km/mean_V__z"+sprinti("%5.5i",lev)+fmtail
 fnames = new((/2,nv/),string)
 do iv=0, nv-1
   fnames(0,iv) = fdir+"/"+vars(iv)+"_z/prt_d300km/"+ \
                    "mean_"+vars(iv)+"__z"+sprinti("%5.5i",lev(iv))+fmtail
   fnames(1,iv) = fdir+"/"+vars(iv)+"_z/prt_d300km/recon_xy/"+ \
                    "prt_"+vars(iv)+"__z"+sprinti("%5.5i",lev(iv))+ftail
 end do
 
 fd = addfile(fnames(1,2),"r")
 fmd = addfile(fnames(0,2),"r")
 x = fd->lon - 40
 y = fd->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 dx = x(1) - x(0)
 dy = y(1) - y(0)
 nt0 = dimsizes( fd->t )
 nt0m = dimsizes( fmd->t )

 tims = round(tims*1.e4,0)*1.e-4

 time = round(fspan(t_09(0),t_09(1),nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(time,tims)
 timem = round(fspan(4,8,nt0m)*1.e4,0)*1.e-4  ; [day]
 itm = get1Dindex(timem,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 nxp = nx*2 + (ny-1)*ni_shift

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   if ( vars(iv) .eq. "U" .or. vars(iv) .eq. "V" ) then
     uuz = addfile(fnameu(iv),"r")
     vvz = addfile(fnamev(iv),"r")
     vs (iv,:,:) = uuz->mean_U(itm,{y},:)
     vs2 = vvz->mean_V(itm,{y},:)
     vs(iv,:,:) = (/ sqrt(vs(iv,:,:)^2 + vs2(:,:)^2) /)
   end if
   if (vars(iv) .eq. "W") then
     fp = addfile(fnames(1,iv),"r")
     vs(iv,:,:) = fp->prt_W(it,{y},:) * 1.e3
   end if 
   if (vars(iv) .eq. "PRS") then
     fm = addfile(fnames(0,iv),"r")
     vs(iv,:,:) = fm->mean_P(itm,{y},:) / 100.
   end if
   if (vars(iv) .eq. "T") then
     fm = addfile(fnames(0,iv),"r")
     vs(iv,:,:) = fm->mean_T(itm,{y},:) + 300.
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
 tday = toint(tims+1.e-4)
 tmp = (tims+1.e-4 - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)
 mstr = "Day "+tday+",  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 if (iw .eq. 0) then
   mstr1 = mstr
 end if
 if (iw .eq. 3) then
   mstr2 = mstr
 end if
 xstr = "longitude [deg]"
 if ( iw .eq. 0 .or. iw .eq. 3 ) then
   ystr = "latitude [deg]"
 end if

 scl = (/1.,1.,1./)
 cnlab = (/0,0,0/)
 lbopt = (/"","",""/)

 cnlev1 = ispan(-10,10,1)*10 + 270.
; cnlev2 = ispan(-10,10,1)*1.5
 cnlev2 = ispan(-5,5,1)*3.
; cntmp1 = 2^(ispan(1,10,1)*0.5)
; cntmp0 = 2^(ispan(10,1,1)*0.5)*(-1.)
; cnlev2(:9 ) = cntmp0
; cnlev2(11:) = cntmp1
 cnlev0 = p_app(0) + ispan(-20,20,1)*8.
 if ( iw .eq. 4 ) then
   lbopt = (/"","[K]","[mm s:S:-1:N:]"/)
 end if
; cnfil = ispan(2,65,3)
; nfil = dimsizes(cnfil)
; cnfil(nfil/2-2:nfil/2+1) = 0
 cnfil = array_append_record( ispan(2,32,6), ispan(35,65,6), 0 )
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0

; plot -----------------------------------------------------
 pscl = 0.8 ;0.45
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xmax - xmin)*rlonlat*pscl
 he = (ymax - ymin)*sqrt( 1. - (ni_shift*rlonlat*dx/dy)^2 )*pscl
 xf = 10 + (/0,1,2,0,1,2/)*(wi+5)
 yf = 20 + (/0,0,0,1,1,1/)*(he+18)

 contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:)/scl(2),cnlev2,cnfil,lbopt(2))
 cnaddline(contour1,wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,cnlab(0), \
           "D",0,"T",0.4,"C",1)
; cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnlab(1), \
;           "D",0,"T",0.5,"C",colind)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.20
;  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.2
  "pmLabelBarHeightF"        : 0.08
  "lbTitlePosition"  : "Right"
  "lbTitleDirection" : "Across"
  "lbLabelFontHeightF" : 0.03
  "lbTitleOffsetF"     : 0.04
 end setvalues

 xylp = iitv0(iw,:)*xyl0(iw)
 plres = True
 plres@gsLineColor      = colind
 plres@gsLineThicknessF = 1.5
 dumm = gsn_add_polyline(wks,contour1,xy0p(iw,0)+(/0,xylp(0)/), \
                                      xy0p(iw,1)+(/0,xylp(1)/),plres)

 set_size(contour1,xf(iw),yf(iw),wi,he)
 axis_range(contour1,xmin,xmax,ymin,ymax,False,False)
 axis_str(contour1,"",xstr,ystr)
 setvalues contour1
  "tmYLMode" : "Manual"
  "tmYLTickStartF" : 20
  "tmYLTickSpacingF"  : 10
  "tmYLMinorPerMajor" : 4
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

 dr_txt(wks,(iw+1),xf(iw)+1.5,yf(iw)+1.5,22,4.5,0.)

 delete(time)

end do  ; global loop

 dr_txt(wks,"(a)",xf(0)-6,yf(0)-3.1,25,6,0.)
 dr_txt(wks,"(b)",xf(0)-6,yf(3)-3.1,25,6,0.)
 dr_txt(wks,mstr1,xf(0)+6,yf(0)-2.8,22,5,0.)
 dr_txt(wks,mstr2,xf(0)+6,yf(3)-2.8,22,5,0.)

frame(wks)

end

