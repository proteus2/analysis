load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf/x6_1hrly/prt_d300km"
; fdname = "domain_x6.nc"
 ftail = "__x6__08_00_00__12_00_00__1hrly.nc"
 vars = (/"GEOP","U"/)
; tims = fspan(4,8,96+1) ;(/6.5/)  ; [day]
; tims = fspan(4,8,9) ;(/6.5/)  ; [day]
 tims = (/8./) ;fspan(4,8,9) ;(/6.5/)  ; [day]
 levs = (/1,5,9,13,17,22/) ;,28/)
 ni_shift = 1
 ymin = 25
 ymax = 75
 
; wks = cr_wks("fig/multilayer0b.ps",600,0)
 wks = cr_wks("fig/multilayer0b-"+tims(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"cmp_b2r")
 gsn_define_colormap(wks,"cmp_flux")
 colind = NhlNewColor(wks,0.,1.,0.)

;--------------------------------------- main body
 nv = dimsizes(vars)
 nt = dimsizes(tims)
 nz = dimsizes(levs)

; fnamed = fdir+"/"+fdname
; fnamep = fdir+"/PRS/prt_PRS__ml"+sprinti("%3.3i",levs)+ftail
; fnamet = fdir+"/T/prt_T__ml"+sprinti("%3.3i",levs)+ftail
; fnamez0 = fdir+"/GEOP/prt_GEOP__mls"+sprinti("%3.3i",levs-1)+ftail
; fnamez1 = fdir+"/GEOP/prt_GEOP__mls"+sprinti("%3.3i",levs  )+ftail
 fnamep = fdir+"/prt_PRS__ml"+sprinti("%3.3i",levs)+ftail
 fnamet = fdir+"/prt_T__ml"+sprinti("%3.3i",levs)+ftail
 fnamez0 = fdir+"/prt_GEOP__mls"+sprinti("%3.3i",levs-1)+ftail
 fnamez1 = fdir+"/prt_GEOP__mls"+sprinti("%3.3i",levs  )+ftail
 fnames = new((/nv,nz/),string)
 fnamess = new((/nv,nz/),string)
 do iv=0, nv-1
   if ( vars(iv) .eq. "W" .or. vars(iv) .eq. "GEOP" ) then
;     fnames(iv,:) = fdir+"/"+vars(iv)+"/"+ \
     fnames(iv,:) = fdir+"/"+ \
                    "prt_"+vars(iv)+"__mls"+sprinti("%3.3i",levs)+ftail
;     fnamess(iv,:) = fdir+"/"+vars(iv)+"/"+ \
     fnamess(iv,:) = fdir+"/"+ \
                    "prt_"+vars(iv)+"__mls"+sprinti("%3.3i",levs-1)+ftail
     if ( vars(iv) .eq. "GEOP" .and. levs(0) .eq. 1 ) then
;       fnamesp = fdir+"/PRS/prt_PRS__ml001"+ftail
       fnamesp = fdir+"/prt_PRS__ml001"+ftail
     end if
   else
;     fnames(iv,:) = fdir+"/"+vars(iv)+"/"+ \
     fnames(iv,:) = fdir+"/"+ \
                    "prt_"+vars(iv)+"__ml"+sprinti("%3.3i",levs)+ftail
   end if
 end do
 
 fd = addfile(fnamep(0),"r")
 x = fd->lon
 y = fd->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 dx = x(1) - x(0)
 dy = y(1) - y(0)

 tims = round(tims*1.e4,0)*1.e-4
 fp = addfile(fnamep(0),"r")
 nt0 = dimsizes( fp->mean_P(:,0,0) )
 time = round(fspan(4,8,nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(time,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 z = new(nz,float)
 p = new(nz,float)
 t = new(nz,float)
 do k=0, nz-1
   fp = addfile(fnamep(k),"r")
   ft = addfile(fnamet(k),"r")
   fz0 = addfile(fnamez0(k),"r")
   fz1 = addfile(fnamez1(k),"r")
   p(k) = avg( fp->mean_P(0,:,:) )/100.
   p2 = fp->mean_P(it,:,:)/100.
   z(k) = avg( 0.5*(fz0->mean_PH(0,:,:) + fz1->mean_PH(0,:,:)) )/9.806/1.e3
   z2 = 0.5*(fz0->mean_PH(it,:,:) + fz1->mean_PH(it,:,:))/9.806/1.e3
   t(k) = avg( ft->mean_T(0,:,:) )+300. * (p(k)/1.e3)^0.286
   print("Z:  "+z(k)+"  "+(exp((z(k)-max(z2))/7.)-1.)+"  "+(exp((z(k)-min(z2))/7.)-1.))
   print("P:  "+p(k)+"  "+(min(p2)/p(k)-1.)+"  "+(max(p2)/p(k)-1.))
 end do
 sq_rho = sqrt( p*100./287./t )

 nxp = nx*2 + (ny-1)*ni_shift

do n=0, nt-1  ; global loop
do k=0, nz-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   f1 = addfile(fnames(iv,k),"r")
   if ( vars(iv) .eq. "U" .or. vars(iv) .eq. "V" ) then
     vs(iv,:,:) = f1->$("mean_"+vars(iv))$(it(n),:,:)
   end if
   if (vars(iv) .eq. "W") then
     f1s = addfile(fnamess(iv,k),"r")
     vs(iv,:,:) = 0.5*( f1->$("prt_"+vars(iv))$(it(n),:,:) + \
                        f1s->$("prt_"+vars(iv))$(it(n),:,:) )*sq_rho(k)*1.e3
   end if 
   if (vars(iv) .eq. "PRS") then
     vs(iv,:,:) = f1->P(it(n),:,:) / 100.
   end if
   if (vars(iv) .eq. "GEOP") then
     if (levs(k) .ne. 1) then
       f1s = addfile(fnamess(iv,k),"r")
       vs(iv,:,:) = 0.5*( f1->mean_PH (it(n),:,:) + \
                          f1s->mean_PH(it(n),:,:) ) / 9.806
     else
       f1sp = addfile(fnamesp,"r")
       vs(iv,:,:) = f1sp->mean_P(it(n),:,:) / 100.
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

 scl = (/1.,1.,1./)
 cnlab = (/1,1,1/)
 lbopt = (/"","[m s:S:-1:N:]","[mm s:S:-1:N:]"/)

 cnlev1 = ispan(-10,10,1)*5
; cnlev1 = (/20,30,40,50/)
 cnlev2 = ispan(-10,10,1)*2.
 if (levs(k) .ne. 1) then
   if (p(k) .gt. 600.) then
     cnlev0 = z(k)*1.e3 + ispan(-20,20,1)*40.
   else
     cnlev0 = z(k)*1.e3 + ispan(-20,20,1)*80.
   end if
 else
   cnlev0 = 1000. + ispan(-20,20,1)*5.
 end if
 if (k .ne. nz-1) then
   lbopt = ""
 end if
 cnfil = ispan(2,23,1)
 nfil = dimsizes(cnfil)
; cnfil(nfil/2-1:nfil/2) = 0
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; plot -----------------------------------------------------
 pscl = 0.7
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xpl(nxp-1) - xpl(0))*rlonlat*pscl
 he = (ymax - ymin)*sqrt( 1. - (ni_shift*rlonlat*dx/dy)^2 )*pscl
 he = he*(ymax - ymin)/(y(ny-1) - y(0))
 xf = 20 + (/0,1,2/)*(wi+6)
 yf = 20 + ispan(nz-1,0,1)*(he-1.5)

; contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:)/scl(2),cnlev2,cnfil,lbopt(2))
 contour1 = cnshadeplot(wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnfil,lbopt(1))
 cnaddline(contour1,wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,cnlab(0), \
           "D",0,"T",0.25,"C",1)
; cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnlab(1), \
;           "D",0,"T",1.2,"C",colind)
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

 xpl1 = x(0)+(ymin-y(0))/dy*dx*ni_shift
 xpl4 = x(0)+(ymax-y(0))/dy*dx*ni_shift
 xpl2 = xpl1 + 2.*(x(nx-1) - x(0))
 xpl3 = xpl4 + 2.*(x(nx-1) - x(0))
 boxres = True
 boxres@gsLineThicknessF = 1.5
 dumm = gsn_add_polyline(wks,contour1,(/xpl1,xpl2,xpl3,xpl4,xpl1/), \
                                      (/ymin,ymin,ymax,ymax,ymin/),boxres)

 draw(contour1)

 if (levs(k) .eq. 1) then
   dr_txt(wks,"sfc",xf(0)-1,yf(k)+he-1.1,22,5.5,0.)
 else
   dr_txt(wks,(round(p(k)/10.,3)*10)+" hPa",xf(0)-4,yf(k)+he-3.8,22,5.5,0.)
   if (z(k) .gt. 10.) then
     dr_txt(wks,"("+round(z(k),3)+" km)",xf(0)-4,yf(k)+he-1.1,22,5.5,0.)
   else
     dr_txt(wks,"("+(round(z(k)*10.,3)/10.)+" km)",xf(0)-4,yf(k)+he-1.1,22,5.5,0.)
   end if
 end if

 delete(contour1)

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop
frame(wks)
end do  ; global loop

end

