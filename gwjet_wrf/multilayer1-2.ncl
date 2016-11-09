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
 vars = (/"GEOP","U","W"/)
; tims = fspan(4,8,96+1) ;(/6.5/)  ; [day]
; tims = fspan(4,8,9) ;(/6.5/)  ; [day]
 tims = (/5.0/) ;fspan(4,8,9) ;(/6.5/)  ; [day]
 levs = (/1,5,9,13,17,22/) ;,28/)
 lev2 = (/1000,920,850,700,500,350/)*1.
 ni_shift = 1
 ymin = 25
 ymax = 75
 
; wks = cr_wks("fig/multilayer0b.ps",600,0)
 wks = cr_wks("fig/multilayer1-2-"+tims(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
; gsn_define_colormap(wks,"cmp_flux")
; gsn_define_colormap(wks,"precip3_16lev")
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
 fnamez = fdir+"/prt_GEOP__mls"+sprinti("%3.3i",ispan(0,49,1))+ftail
 fnamew = fdir+"/prt_W__mls"+sprinti("%3.3i",ispan(0,49,1))+ftail
 fnameu = fdir+"/prt_U__ml"+sprinti("%3.3i",ispan(1,49,1))+ftail
 fnamev = fdir+"/prt_V__ml"+sprinti("%3.3i",ispan(1,49,1))+ftail
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

 nay = 26
 nax = toint(nay/cos(y*3.141592/180.))
 nay = nay*2+1
 nax = nax*2+1

 tims = round(tims*1.e4,0)*1.e-4
 fp = addfile(fnamep(0),"r")
 nt0 = dimsizes( fp->mean_P(:,0,0) )
 time = round(fspan(4,8,nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(time,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 z = new(nz,float)
 t = new(nz,float)
 do k=0, nz-1
   ft = addfile(fnamet(k),"r")
   t(k) = avg( ft->mean_T(0,:,:) )+300. * (lev2(k)/1.e3)^0.286
 end do
 sq_rho = sqrt( lev2*100./287./t )

 nxp = nx*2 + (ny-1)*ni_shift

lnlev2 = log(lev2)
ez = (/1, 0.993, 0.983, 0.97, 0.954, 0.934, 0.909, 0.88, 0.8514746, 0.8229493, \
    0.7944239, 0.7658985, 0.7135208, 0.6640677, 0.6174032, 0.5733964, \
    0.5319219, 0.4928589, 0.4560919, 0.4215095, 0.3890053, 0.3584769, \
    0.3298262, 0.3029593, 0.2777859, 0.2542197, 0.2321779, 0.2115813, \
    0.1923539, 0.1744233, 0.1577198, 0.1421769, 0.1277312, 0.1143217, \
    0.1018904, 0.09038174, 0.0797426, 0.06992226, 0.06087225, 0.05254621, \
    0.04489981, 0.03789053, 0.03147761, 0.0256217, 0.02028462, 0.01542884, \
    0.01101657, 0.007007994, 0.003357153, 0/)
ep = (/0.9965, 0.988, 0.9765, 0.962, 0.944, 0.9215, 0.8945, 0.8657373, 0.837212, \
    0.8086866, 0.7801613, 0.7397097, 0.6887943, 0.6407355, 0.5953999, \
    0.5526592, 0.5123904, 0.4744754, 0.4388007, 0.4052574, 0.3737411, \
    0.3441516, 0.3163928, 0.2903726, 0.2660028, 0.2431988, 0.2218796, \
    0.2019676, 0.1833886, 0.1660715, 0.1499483, 0.1349541, 0.1210264, \
    0.1081061, 0.09613608, 0.08506217, 0.07483242, 0.06539726, 0.05670923, \
    0.04872301, 0.04139517, 0.03468407, 0.02854965, 0.02295316, 0.01785673, \
    0.0132227, 0.00901228, 0.005182574, 0.001678576/)
coefz1 = new(50,float)
do k=1, 48
  coefz1(k) = (ez(k) - ep(k))/(ep(k-1) - ep(k))
end do
coefz2 = 1. - coefz1

do n=0, nt-1  ; global loop
 pa = new((/49,ny,nx/),float)
 do k=0, 48
  fpname = fdir+"/prt_PRS__ml"+sprinti("%3.3i",k+1)+ftail
  ffp = addfile(fpname,"r")
  pa(k,:,:) = ffp->mean_P(it(n),:,:) + ffp->prt_P(it(n),:,:)
 end do
 pa = pa/100.
 lnpa = log(pa)
 lnpaz = new((/50,ny,nx/),float)
 do k=1, 48
;  paz(k,:,:) = (/wrf_interp_3d_z(pa,conform_dims((/49,ny,nx/),ep,0),ez(k)), \
  lnpaz(k,:,:) = coefz1(k)*pa(k-1,:,:) + coefz2(k)*pa(k,:,:)
 end do
 lnpaz = log(lnpaz)
 zz = new((/50,ny,nx/),float)
 ww = new((/50,ny,nx/),float)
 uu = new((/49,ny,nx/),float)
 vv = new((/49,ny,nx/),float)
 do k=0, 49
   ffz = addfile(fnamez(k),"r")
   ffw = addfile(fnamew(k),"r")
;   zz(k,:,:) = ffz->mean_PH(it(n),:,:)/9.806
   zz(k,:,:) = (ffz->mean_PH(it(n),:,:)+ffz->prt_PH(it(n),:,:))/9.806
   ww(k,:,:) = (ffw->mean_W (it(n),:,:)+ffw->prt_W (it(n),:,:))
 end do
 do k=0, 48
   uuz = addfile(fnameu(k),"r")
   vvz = addfile(fnamev(k),"r")
;   zz(k,:,:) = ffz->mean_PH(it(n),:,:)/9.806
   uu(k,:,:) = (uuz->mean_U(it(n),:,:)+uuz->prt_U(it(n),:,:))
   vv(k,:,:) = (vvz->mean_V(it(n),:,:)+vvz->prt_V(it(n),:,:))
 end do 
do k=0, nz-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 vs = new((/nv,ny,nx/),float)
 vs2 = new((/nv,ny,nx/),float)
 do iv=0, nv-1
   f1 = addfile(fnames(iv,k),"r")
   if ( vars(iv) .eq. "U" .or. vars(iv) .eq. "V" ) then
     if (lev2(k) .ne. 1000.) then
       vs(iv,:,:) = (/ wrf_interp_3d_z(uu,lnpa,lnlev2(k)) /)
       vs2(iv,:,:) = (/ wrf_interp_3d_z(vv,lnpa,lnlev2(k)) /)
       mmm = new((/ny,nx/),float)
       do j=0, ny-1
         mmm(j,:) = runave(vs(iv,j,:),nax(j),-1)
       end do
       mmm = runave_n(mmm,nay,0,0)
       vs(iv,:,:) = mmm
       do j=0, ny-1
         mmm(j,:) = runave(vs2(iv,j,:),nax(j),-1)
       end do
       mmm = runave_n(mmm,nay,0,0)
       vs2(iv,:,:) = mmm
       delete(mmm)
       vs(iv,:,:) = (/ sqrt(vs(iv,:,:)^2 + vs2(iv,:,:)^2) /)
     else
       vs(iv,:,:) = (/ f1->$("mean_"+vars(iv))$(it(n),:,:) /)
       fnamesv = fdir+"/prt_V__ml"+sprinti("%3.3i",levs(k))+ftail
       f2 = addfile(fnamesv,"r")
       vs(iv,:,:) = (/ sqrt(vs(iv,:,:)^2 + f2->mean_V(it(n),:,:)^2) /)
     end if
   end if
   if (vars(iv) .eq. "W") then
     if (lev2(k) .ne. 1000.) then
       vs(iv,:,:) = (/ wrf_interp_3d_z(ww(1:48,:,:),lnpaz(1:48,:,:),lnlev2(k))*sq_rho(k)*1.e3 /)
       mmm = new((/ny,nx/),float)
       do j=0, ny-1
         mmm(j,:) = runave(vs(iv,j,:),nax(j),-1)
       end do
       mmm = runave_n(mmm,nay,0,0)
       vs(iv,:,:) = vs(iv,:,:) - mmm
       delete(mmm)
     else
       f1s = addfile(fnamess(iv,k),"r")
       vs(iv,:,:) = 0.5*( f1->$("prt_"+vars(iv))$(it(n),:,:) + \
                          f1s->$("prt_"+vars(iv))$(it(n),:,:) )*sq_rho(k)*1.e3
     end if
   end if 
   if (vars(iv) .eq. "PRS") then
     vs(iv,:,:) = f1->P(it(n),:,:) / 100.
   end if
   if (vars(iv) .eq. "GEOP") then
     if (lev2(k) .ne. 1000.) then
;       f1s = addfile(fnamess(iv,k),"r")
       vs(iv,:,:) = (/ wrf_interp_3d_z(zz(1:48,:,:),lnpaz(1:48,:,:),lnlev2(k)) /)
       mmm = new((/ny,nx/),float)
       do j=0, ny-1
         mmm(j,:) = runave(vs(iv,j,:),nax(j),-1)
       end do
       mmm = runave_n(mmm,nay,0,0)
       vs(iv,:,:) = mmm
       delete(mmm)
       z(k) = avg(vs(iv,:,:))/1.e3
print(z(k)+"")
;       vs(iv,:,:) = 0.5*( f1->mean_PH (it(n),:,:) + \
;                          f1s->mean_PH(it(n),:,:) ) / 9.806
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

; cnlev1 = ispan(0,15,1)*3.
 cnlev1 = (/20,30,40,50/)
 cnlev2 = ispan(-10,10,1)*2.
 if (lev2(k) .ne. 1000.) then
   if (lev2(k) .gt. 600.) then
     cnlev0 = z(k)*1.e3 + ispan(-20,20,1)*80.;60.
   else
     cnlev0 = z(k)*1.e3 + ispan(-20,20,1)*80.
   end if
 else
   cnlev0 = 1000. + ispan(-20,20,1)*8.
 end if
 if (k .ne. nz-1) then
   lbopt = ""
 end if
 cnfil = ispan(2,65,3)
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

 contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:)/scl(2),cnlev2,cnfil,lbopt(2))
 cnaddline(contour1,wks,xpl,y,vpl(0,:,:)/scl(0),cnlev0,cnlab(0), \
           "D",0,"T",0.25,"C",1)
 cnaddline(contour1,wks,xpl,y,vpl(1,:,:)/scl(1),cnlev1,cnlab(1), \
           "D",0,"T",0.6,"C",colind)
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

 if (lev2(k) .eq. 1000.) then
   dr_txt(wks,"sfc",xf(0)-1,yf(k)+he-1.1,22,5.5,0.)
 else
   dr_txt(wks,toint(lev2(k))+" hPa",xf(0)-4,yf(k)+he-1.1,22,5.5,0.)
 end if

 delete(contour1)

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop
frame(wks)
end do  ; global loop

end

