load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 fdir = "/data18/GW_jet_wrf/dat/x6_5min"
 fdname = "domain_x6.nc"
 ftail = "__x6__08_00_00__12_00_00__5min.nc"
 vars = (/"U","GEOP"/)
 tims = fspan(4,8,96+1) ;(/6.5/)  ; [day]
 levs = (/1,9,17,22,28/)
 ni_shift = 1
 ymin = 25
 ymax = 75
 
 wks = cr_wks("fig/multilayer.ps",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- main body
 nv = dimsizes(vars)
 nt = dimsizes(tims)
 nz = dimsizes(levs)

 fnamed = fdir+"/"+fdname
 fnamep = fdir+"/PRS/PRS__ml"+sprinti("%3.3i",levs)+ftail
 fnamez0 = fdir+"/GEOP/GEOP__mls"+sprinti("%3.3i",levs-1)+ftail
 fnamez1 = fdir+"/GEOP/GEOP__mls"+sprinti("%3.3i",levs  )+ftail
 fnames = new((/nv,nz/),string)
 fnamess = new((/nv,nz/),string)
 do iv=0, nv-1
   if ( vars(iv) .eq. "W" .or. vars(iv) .eq. "GEOP" ) then
     fnames(iv,:) = fdir+"/"+vars(iv)+"/"+ \
                    vars(iv)+"__mls"+sprinti("%3.3i",levs)+ftail
     fnamess(iv,:) = fdir+"/"+vars(iv)+"/"+ \
                    vars(iv)+"__mls"+sprinti("%3.3i",levs-1)+ftail
     if ( vars(iv) .eq. "GEOP" .and. levs(0) .eq. 1 ) then
       fnamesp = fdir+"/PRS/PRS__ml001"+ftail
     end if
   else
     fnames(iv,:) = fdir+"/"+vars(iv)+"/"+ \
                    vars(iv)+"__ml"+sprinti("%3.3i",levs)+ftail
   end if
 end do
 
 fd = addfile(fnamed,"r")
 x = fd->XLONG(0,:)  ; uniform lat-lon grid
 y = fd->XLAT (:,0)  ; uniform lat-lon grid
 nx = dimsizes(x)
 ny = dimsizes(y)

 tims = round(tims*1.e4,0)*1.e-4
 f1 = addfile(fnames(0,0),"r")
 nt0 = dimsizes( f1->$vars(0)$(:,0,0))
 t = round(fspan(4,8,nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(t,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 z = new(nz,float)
 p = new(nz,float)
 do k=0, nz-1
   fp = addfile(fnamep(k),"r")
   fz0 = addfile(fnamez0(k),"r")
   fz1 = addfile(fnamez1(k),"r")
   p(k) = avg( fp->P(0,:,:) )/100.
   p2 = fp->P(it,:,:)/100.
   z(k) = avg( 0.5*(fz0->PH(0,:,:) + fz1->PH(0,:,:)) )/9.806/1.e3
   z2 = 0.5*(fz0->PH(it,:,:) + fz1->PH(it,:,:))/9.806/1.e3
   print("Z:  "+z(k)+"  "+(exp((z(k)-max(z2))/7.)-1.)+"  "+(exp((z(k)-min(z2))/7.)-1.))
   print("P:  "+p(k)+"  "+(min(p2)/p(k)-1.)+"  "+(max(p2)/p(k)-1.))
 end do

 nxp = nx*2 + (ny-1)*ni_shift

do n=0, nt-1  ; global loop
do k=0, nz-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   f1 = addfile(fnames(iv,k),"r")
   if (vars(iv) .eq. "U") then
     vs(iv,:,:) = 0.5*( f1->$vars(iv)$(it(n),:,0:nx-1) + \
                        f1->$vars(iv)$(it(n),:,1:nx  ) )
   end if
   if (vars(iv) .eq. "V") then
     vs(iv,:,:) = 0.5*( f1->$vars(iv)$(it(n),0:ny-1,:) + \
                        f1->$vars(iv)$(it(n),1:ny  ,:) )
   end if
   if (vars(iv) .eq. "PRS") then
     vs(iv,:,:) = f1->P(it(n),:,:) / 100.
   end if
   if (vars(iv) .eq. "GEOP") then
     if (levs(k) .ne. 1) then
       f1s = addfile(fnamess(iv,k),"r")
       vs(iv,:,:) = 0.5*( f1->PH (it(n),:,:) + \
                          f1s->PH(it(n),:,:) ) / 9.806
     else
       f1sp = addfile(fnamesp,"r")
       vs(iv,:,:) = f1sp->P(it(n),:,:) / 100.
     end if
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

; if (vars(0) .eq. "W") then
;   scl = 1.e-3
;   cnlev1 = ispan(-10,10,1)*4
;   cnlev2 = ispan(-10,10,1)*4
;   cnlab = 0
;   lbopt = "[mm s:S:-1:N:]"
; end if
 if (vars(0) .eq. "U") then
   scl = 1.
;   cnlev1 = ispan(-10,10,1)*0.25
;   cnlev2 = ispan(1,20,1)*3
   cnlab = 1
   lbopt = "[m s:S:-1:N:]"
 end if
; if (vars(0) .eq. "V") then
;   scl = 1.
;   cnlev1 = ispan(-10,10,1)*0.25
;   cnlev2 = ispan(-10,10,1)*2
;   cnlab = 0
;   lbopt = "[m s:S:-1:N:]"
; end if
 cnlev1 = ispan(-10,10,1)*4
 if (levs(k) .ne. 1) then
   cnlev2 = z(k)*1.e3 + ispan(-10,10,1)*60.
 else
   cnlev2 = 1000. + ispan(-10,10,1)*5.
 end if
 if (k .ne. nz-1) then
   lbopt = ""
 end if
 cnfil = ispan(2,65,3)
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; plot -----------------------------------------------------
 pscl = 0.7
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xpl(nxp-1) - xpl(0))*rlonlat*pscl
 he = (ymax - ymin)*sqrt( 1. - (ni_shift*rlonlat*dx/dy)^2 )*pscl
 he = he*(ymax - ymin)/(y(ny-1) - y(0))
 xf = 20 + (/0,1,2/)*(wi+6)
 yf = 20 + ispan(nz-1,0,1)*(he-1.5)

 contour1 = cnshadeplot(wks,xpl,y,vpl(0,:,:)/scl,cnlev1,cnfil,lbopt)
 cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl,cnlev2,cnlab, \
           "D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
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
 set_size(contour1,xf(0),yf(k),wi,he)
 axis_range(contour1,-999,999,ymin,ymax,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 if (levs(k) .eq. 1) then
   dr_txt(wks,"sfc",xf(0)-1,yf(k)+he-1.1,22,5.5,0.)
 else
   dr_txt(wks,(round(p(k)/10.,3)*10)+" hPa",xf(0)-4,yf(k)+he-3.8,22,5.5,0.)
   dr_txt(wks,"(~"+round(z(k),3)+" km)",xf(0)-4,yf(k)+he-1.1,22,5.5,0.)
 end if

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop
frame(wks)
end do  ; global loop

end

; dumm = new(5,graphic)
; dumm(0) = gsn_add_polyline(wks,contour1,(/31,35/),(/52,48/),False)
; dumm(1) = gsn_add_polyline(wks,contour1,(/16,24/),(/61.4,61.4/),False)
; dumm(2) = gsn_add_polyline(wks,contour1,(/4,8/),(/58.8,62.8/),False)
; dumm(3) = gsn_add_polyline(wks,contour1,(/13,17/),(/40.5,44.5/),False)
; dumm(4) = gsn_add_polyline(wks,contour1,(/10,11.5/),(/44,42.5/),False)

; delete(contour1)

; dr_txt(wks,mstr,xf(0)+wi/2.,yf(0)-2.5,22,5.5,0.)
; dr_txt(wks,"W1",xf(0)+wi-2,yf(0)+he/2+3,21,4,0.)
; dr_txt(wks,"W2",xf(0)+wi/2+1,yf(0)+he/5,21,4,0.)
; dr_txt(wks,"W3",xf(0)+wi/4,yf(0)+he/5,21,4,0.)
; dr_txt(wks,"W4",xf(0)+wi/3,yf(0)+he/5*4-2,21,4,0.)
; dr_txt(wks,"W5",xf(0)+wi/3-3.5,yf(0)+he/5*4-4.5,21,4,0.)
