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
 lev = 8000
 lons = (/ 250, 150, 400 ,400, 280 /)
 ts   = (/ (/-2.,2./), (/-2.,2./), (/-4.,4./), (/-12,6/), (/-6,6/) /)
 ymin = 30
 ymax = 70
 xmin = 20
 xmax = 60

 lon2 = (/44,38,54,37,31/)
 lat2 = (/52,60,62,42,44/)
 
; wks = cr_wks("fig/multilayer_z.ps",600,0)
 wks = cr_wks("fig/zoomin-4",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
; gsn_define_colormap(wks,"cmp_flux")
; gsn_define_colormap(wks,"precip3_16lev")
 colind = NhlNewColor(wks,0.,1.,0.)

;--------------------------------------- main body
 nys = dimsizes(lons)

 fname = fdir+"/W_z/prt_d300km/"+ \
          "prt_W__z"+sprinti("%5.5i",lev)+ftail
 fnamep = fdir+"/PRS_z/prt_d300km/"+ \
          "mean_PRS__z"+sprinti("%5.5i",250)+ftail
 fnameu = fdir+"/U_z/prt_d300km/"+ \
          "mean_U__z"+sprinti("%5.5i",lev)+ftail
 fnamev = fdir+"/V_z/prt_d300km/"+ \
          "mean_V__z"+sprinti("%5.5i",lev)+ftail
 fnamet1 = fdir+"/T_z/prt_d300km/"+ \
           "mean_T__z"+sprinti("%5.5i",lev-250)+ftail
 fnamet2 = fdir+"/T_z/prt_d300km/"+ \
           "mean_T__z"+sprinti("%5.5i",lev+250)+ftail
 
 fd = addfile(fname,"r")
 x = fd->lon - 40
 y = fd->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 dx = x(1) - x(0)
 dy = y(1) - y(0)
 nt0 = dimsizes( fd->t )/2+1
 t = fspan(6,8,nt0)

 time = round(fspan(6,8,nt0)*1.e4,0)*1.e-4  ; [day]

 it = ispan(nt0-1,nt0*2-2,1)

print("No. of lon, lat :  "+nx+", "+ny)

 nxp = nx*2

; read -----------------------------------------------------
 fw = addfile(fname,"r")
 fu = addfile(fnameu,"r")
 fv = addfile(fnamev,"r")
 fp = addfile(fnamep,"r")
 ft1 = addfile(fnamet1,"r")
 ft2 = addfile(fnamet2,"r")
 vsa = fw->prt_W(it,:,:) * 1.e3
 pa = fp->mean_P(it,:,:) * 1.e3

; process --------------------------------------------------
 xpl = new(nxp,float,1.e20)
 xpl(:nx-1) = x(:)
 xpl(nx:) = xpl(nx-1) + (x(1)-x(0))*ispan(1,nxp-nx,1)
 xpl!0 = "x"
 xpl&x = xpl

 vpl = new((/nys,nt0,ny,nxp/),float,1.e20)
 ppl = new((/nys,nt0,ny,nxp/),float,1.e20)

; vsx = vs
; vsx = 0.
; do n=0, nt0-1
; do j=0, nys-1
;   vsx(n,j,:) = fspan(vs(n,j,0),vs(n,j,nx-1),nx) - (vs(n,j,0)+vs(n,j,nx-1))*0.5
; end do
; end do
;vs = vs - vsx

 do j=0, nys-1
   vpl(j,:,:,:nx-1) = vsa(:,:,:)
   vpl(j,:,:,nx:nx+nx-1) = vsa(:,:,:)
   ppl(j,:,:,1:nx-2) = pa(:,:,2:nx-1) - pa(:,:,0:nx-3)
   ppl(j,:,:,0) = pa(:,:,1) - pa(:,:,nx-1)
   ppl(j,:,:,nx-1) = pa(:,:,0) - pa(:,:,nx-2)
   ppl(j,:,:,nx:nx+nx-1) = ppl(j,:,:,:nx-1)
 end do
 vpl!1 = "t"
 vpl!2 = "y"
 vpl!3 = "x"
 vpl&t = t
 vpl&y = y
 vpl&x = xpl

 vpl(0,:,{:48},:) = 0.
 vpl(0,:,{58:},:) = 0.
 vpl(1,:,{:54},:) = 0.
 vpl(1,:,{64:},:) = 0.
 vpl(2,:,{:55},:) = 0.
 vpl(2,:,{65:},:) = 0.
 vpl(3,:,{:36},:) = 0.
 vpl(3,:,{46:},:) = 0.
 vpl(4,:,{:38},:) = 0.
 vpl(4,:,{48:},:) = 0.

 vpl(0,:,:,:) = where(ppl(0,:,:,:) .lt. 0., 0., vpl(0,:,:,:))
 vpl(1,:,:,:) = where(ppl(1,:,:,:) .gt. 0., 0., vpl(1,:,:,:))
 vpl(2,:,:,:) = where(ppl(2,:,:,:) .lt. 0., 0., vpl(2,:,:,:))
 vpl(3,:,:,:) = where(ppl(3,:,:,:) .lt. 0., 0., vpl(3,:,:,:))
 vpl(4,:,:,:) = where(ppl(4,:,:,:) .gt. 0., 0., vpl(4,:,:,:))
; vpl = runave_Wrap(vpl,3,-1)
; vpl = runave_n_Wrap(vpl,3,1,1)
 vpl(2,:,:,:nx+60) = 0.
 vpl(2,:,:,nx+nx-1-60:) = 0.

 vpl0 = new((/nys,ny,nxp/),float,1.e20)
 time0 = new(nys,float)
 maxlon = new(nys,float)
 maxlat = new(nys,float)
 maxlat  = new(nys,float)
 dl = new((/nys,2/),float)

 t0 = (/6.75,6.75,6.75,7.25,7.25/)
 dt = 0.75

 do j=0, nys-1

   it0 = closest_val(t0(j),t)
   av1d = ndtooned(abs(vpl(j,{t0(j):t0(j)+dt},:,nx:nx+nx-1)))
   dsizes = dimsizes(vpl(j,{t0(j):t0(j)+dt},:,nx:nx+nx-1))
   ij = ind_resolve(maxind(av1d),dsizes)
   time0(j) = t(ij(0,0)+it0)
   maxlat(j) = y(ij(0,1))
   maxlon(j) = xpl(ij(0,2)+nx)
   print(max(abs(vpl(j,{t0(j):t0(j)+dt},:,nx:nx+nx-1)))+" at "+maxlon(j)+" (t = "+time0(j)+")")

   vs0 = fw->prt_W(nt0-1+ij(0,0)+it0,:,:) * 1.e3
   do jj=0, ny-1
     vpl0(j,jj,:nx-1) = vs0(jj,:)
     vpl0(j,jj,nx:nx+nx-1) = vs0(jj,:)
   end do

   lhs = fspan(0.,40.,432)
   lhc = lhs(closest_val(4.,lhs))
   vsh = fw->prt_W(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc:maxlon(j)+lhc}) * 1.e3
   ushu = fu->mean_U(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc:maxlon(j)+lhc})
   ushv = fv->mean_V(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc:maxlon(j)+lhc})
   nsh = sqrt( log( ( ft2->mean_T(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc:maxlon(j)+lhc}) + 300. )/  \
                    ( ft1->mean_T(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc:maxlon(j)+lhc}) + 300. )  \
                  )/500. * 9.806 )

   nd = dimsizes(vsh)

   if (j .eq. 0) then
     nh = nd(1)
     vplh = new((/nys,nt0,nh/),float,1.e20)
     uplh = new((/nys,nt0,nh/),float,1.e20)
     nplh = new((/nys,nt0,nh/),float,1.e20)
     xs = (xpl({maxlon(j)-lhc:maxlon(j)+lhc}) - maxlon(j))*111.
     hh = new((/nys,nh/),float)
   end if
   if (nh .ne. nd(2)) then
     print("smaller xdim : "+nd(2)+" < "+nh+" for W"+(j+1))
     if (maxlon(j)+lhc .gt. 60.) then
       vsh2 = fw->prt_W(it,{maxlat(j)-lhc:maxlat(j)+lhc},{20.:maxlon(j)+lhc-40.}) * 1.e3
;       ushu2 = fu->mean_U(it,{maxlat(j)-lhc:maxlat(j)+lhc},{20.:maxlon(j)+lhc-40.})
;       ushv2 = fv->mean_V(it,{maxlat(j)-lhc:maxlat(j)+lhc},{20.:maxlon(j)+lhc-40.})
     else
       vsh2 = fw->prt_W(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc+40.:60.}) * 1.e3
;       ushu2 = fu->mean_U(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc+40.:60.})
;       ushv2 = fv->mean_V(it,{maxlat(j)-lhc:maxlat(j)+lhc},{maxlon(j)-lhc+40.:60.})
     end if
   end if
 
   if (j .eq. 0) then
     if ( (nh/2)/2*2 .eq. (nh/2) ) then
       i9 = nh/2
     else
       i9 = nh/2-1
     end if
     do i=-i9, i9-1, 2
       vplh(j,:,i+nh/2) = vsh(:,nh/2-i/2,i+nh/2)
       uplh(j,:,i+nh/2) = ushu(:,nh/2-i/2,i+nh/2)*(1.0 /sqrt(1.^2+0.5^2)) +  \
                          ushv(:,nh/2-i/2,i+nh/2)*(-0.5/sqrt(1.^2+0.5^2))
       nplh(j,:,i+nh/2) = nsh(:,nh/2-i/2,i+nh/2)
     end do
     do i=-i9, i9-1, 2
       vplh(j,:,i+nh/2+1) = 0.5*(vsh(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 vsh(:,nh/2-(i/2+1),i+nh/2+1))
       uplh(j,:,i+nh/2+1) = 0.5*(ushu(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 ushu(:,nh/2-(i/2+1),i+nh/2+1))* \
                                (1.0 /sqrt(1.^2+0.5^2)) +  \
                            0.5*(ushv(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 ushv(:,nh/2-(i/2+1),i+nh/2+1))* \
                                (-0.5/sqrt(1.^2+0.5^2))
       nplh(j,:,i+nh/2+1) = 0.5*(nsh(:,nh/2-i/2    ,i+nh/2+1)+ \
                                  nsh(:,nh/2-(i/2+1),i+nh/2+1))
     end do
     hha = sqrt((0.5*xs)^2 + (xs*cos(maxlat(j)*3.141592/180.))^2)
     dl(j,:) = (/1,-0.5/)
   end if
   if (j .eq. 1) then
     do i=0, nh-1
       vplh(j,:,i) = vsh(:,nh/2,i)
       uplh(j,:,i) = ushu(:,nh/2,i)
       nplh(j,:,i) = nsh(:,nh/2,i)
     end do
     hha = abs(xs*cos(maxlat(j)*3.141592/180.))
     dl(j,:) = (/1,0/)
   end if
   if (j .eq. 2) then
;     if (maxlon(j)+lhc .gt. 60.) then
;       do i=0, nd(2)-1
;         vplh(j,:,i) = vsh(:,i,i)
;       end do
;       do i=nd(2), nh-1
;         vplh(j,:,i) = vsh2(:,i,i-nd(2))
;       end do
;     else
;       print("correct this.")
;;       do i=0, nd(2)-1
       do i=0, nh-1
         vplh(j,:,i) = vsh(:,i,i)
         uplh(j,:,i) = ushu(:,i,i)*(1.0/sqrt(1.^2+1.^2)) +  \
                       ushv(:,i,i)*(1.0/sqrt(1.^2+1.^2))
         nplh(j,:,i) = nsh(:,i,i)
       end do
;       do i=nd(2), nh-1
;         vplh(j,:,i) = vsh2(:,i,i-nd(2))
;       end do
;     end if
     hha = sqrt(xs^2 + (xs*cos(maxlat(j)*3.141592/180.))^2)
     dl(j,:) = (/1,1/)
   end if
   if (j .eq. 3) then
     do i=0, nh-1
       vplh(j,:,i) = vsh(:,i,i)
       uplh(j,:,i) = ushu(:,i,i)*(1.0/sqrt(1.^2+1.^2)) +  \
                     ushv(:,i,i)*(1.0/sqrt(1.^2+1.^2))
       nplh(j,:,i) = nsh(:,i,i)
     end do
     hha = sqrt(xs^2 + (xs*cos(maxlat(j)*3.141592/180.))^2)
     dl(j,:) = (/1,1/)
   end if
   if (j .eq. 4) then
     if ( (nh/2)/2*2 .eq. (nh/2) ) then
       i9 = nh/2
     else
       i9 = nh/2-1
     end if
     do i=-i9, i9-1, 2
       vplh(j,:,i+nh/2) = vsh(:,nh/2-i/2,i+nh/2)
       uplh(j,:,i+nh/2) = ushu(:,nh/2-i/2,i+nh/2)*(1.0 /sqrt(1.^2+0.5^2)) +  \
                          ushv(:,nh/2-i/2,i+nh/2)*(-0.5/sqrt(1.^2+0.5^2))
       nplh(j,:,i+nh/2) = nsh(:,nh/2-i/2,i+nh/2)
     end do
     do i=-i9, i9-1, 2
       vplh(j,:,i+nh/2+1) = 0.5*(vsh(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 vsh(:,nh/2-(i/2+1),i+nh/2+1))
       uplh(j,:,i+nh/2+1) = 0.5*(ushu(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 ushu(:,nh/2-(i/2+1),i+nh/2+1))* \
                                (1.0 /sqrt(1.^2+0.5^2)) +  \
                            0.5*(ushv(:,nh/2-i/2    ,i+nh/2+1)+ \
                                 ushv(:,nh/2-(i/2+1),i+nh/2+1))* \
                                (-0.5/sqrt(1.^2+0.5^2))
       nplh(j,:,i+nh/2+1) = 0.5*(nsh(:,nh/2-i/2    ,i+nh/2+1)+ \
                                  nsh(:,nh/2-(i/2+1),i+nh/2+1))
     end do
     hha = sqrt((0.5*xs)^2 + (xs*cos(maxlat(j)*3.141592/180.))^2)
     dl(j,:) = (/1,-0.5/)
   end if
 
   hh(j,:) = where(xs .ge. 0., hha, -hha)

   dl(j,0) = dl(j,0)*tofloat(lhc)*lons(j)/hh(j,nh-1)
   dl(j,1) = dl(j,1)*tofloat(lhc)*lons(j)/hh(j,nh-1)

   delete( vsh )

 end do

 vplh!1 = "t"
 vplh&t = t
 uplh!1 = "t"
 uplh&t = t
 nplh!1 = "t"
 nplh&t = t
 
; plot setting ---------------------------------------------
 tday = toint(time0+1.e-4)
 tmp = (time0+1.e-4 - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)
 mstr = "Day "+tday+",  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 xstr1 = (/"","","","","longitude [deg]"/)
 ystr1 = "latitude [deg]"
 xstr2 = (/"","","","","distance [km]"/)
 ystr2 = "time [hr]"

 lbopt1 = (/"","","","","[mm s:S:-1:N:]"/)
 lbopt2 = (/"","","",""," "/)
; cnfil = ispan(2,65,3)
; nfil = dimsizes(cnfil)
; cnfil(nfil/2-2:nfil/2+1) = 0
 cnfil = array_append_record( ispan(2,32,6), ispan(35,65,6), 0 )
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0

 cnfil2 = array_append_record( ispan(2,26,6), ispan(35,65,6), 0 )
 nfil2 = dimsizes(cnfil2)
 cnfil2(nfil2/2) = 0

; plot -----------------------------------------------------
 pscl = 0.6 ;0.45
 rlonlat = 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi0 = (xmax - xmin)*rlonlat*pscl
 he0 = (ymax - ymin)*pscl
 xf0 = 10
 yf0 = 20 + ispan(0,nys-1,1)*(he0+5)

 scl = (/1,1,1,5,1/)
do j=0, nys-1

 cnlev = ispan(-5,5,1)*10.
; if (j .eq. 0) then
;   cnlev = (/-100,-50,-20,-10,-5,0,5,10,20,50,100/)
;   cnlev = ispan(-5,5,1)*15.
; end if
 cnlev2 = (/-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9/)
 cnlev2b = ispan(-40,40,1)
 cnlev2n = ispan(120,180,5)*0.01

 contour1 = cnshadeplot(wks,xpl,y,vpl0(j,:,:),cnlev,cnfil,lbopt1(j))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.6
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.8
  "pmLabelBarHeightF"        : 0.1*0.5
  "lbTitlePosition"  : "right"
  "lbTitleDirection" : "across"
  "lbLabelFontHeightF" : 0.04
  "lbTitleOffsetF"     : 0.04
 end setvalues
 set_size(contour1,xf0,yf0(j),wi0,he0)
 axis_range(contour1,xmin,xmax,ymin,ymax,False,False)
 axis_str(contour1,"",xstr1(j),ystr1)

 if (j .eq. 0) then
   dr_txt(wks,"(a)",xf0-4,yf0(j)-3.0,25,6,0.)
 end if

 ly = 4

 xf = xf0 + (wi0+12)

 he = he0
 wi = he
 yf = yf0(j)
 lx = tofloat(wi)/he0*ly/rlonlat

 contour3 = cnshadeplot(wks,xpl,y,vpl0(j,:,:),cnlev,cnfil,"")
 set_size(contour3,xf,yf,wi,he)
 axis_range(contour3,lon2(j)-lx,lon2(j)+lx,lat2(j)-ly,lat2(j)+ly,False,False)
 axis_str(contour3,"",xstr1(j),ystr1)

 if (j .eq. 0) then
   dr_txt(wks,"(b)",xf-4,yf-3.0,25,6,0.)
 end if

 xf = xf + (wi+12)

 wi = wi0*2.5
 he = he0
 yf = yf0(j)

 contour2 = cnshadeplot(wks,hh(j,:),(t-time0(j))*24,vplh(j,:,:)/max(abs(vplh(j,{ts(j,0)/24.+time0(j):ts(j,1)/24.+time0(j)},:)))*scl(j),cnlev2,cnfil2,lbopt2(j))
 cnaddline(contour2,wks,hh(j,:),(t-time0(j))*24,uplh(j,:,:),cnlev2b,1,"D",(/0,0,14/),"T",(/0.5,1.0,0.5/),"C",1)
 cnaddline(contour2,wks,hh(j,:),(t-time0(j))*24,nplh(j,:,:)*1.e2,cnlev2n,1,"D",14,"T",0.5,"C",1)

 setvalues contour2
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 0.6
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
  "pmLabelBarHeightF"        : 0.1*0.8
  "lbTitlePosition"   : "right"
  "lbTitleDirection"  : "across"
  "lbLabelAutoStride" : False
  "lbLabelStride"     : 1
  "lbLabelFontHeightF" : 0.04*0.4
  "lbTitleOffsetF"     : 0.04
 end setvalues
 set_size(contour2,xf,yf,wi,he)
 axis_range(contour2,-lons(j),lons(j),ts(j,0),ts(j,1),False,False)
 axis_str(contour2,"",xstr2(j),ystr2)

 if (j .eq. 0) then
   dr_txt(wks,"(c)",xf-4,yf-3.0,25,6,0.)
 end if

 plres = True
 plres@gsLineColor      = colind
 plres@gsLineThicknessF = 1.5
 dumm = gsn_add_polyline(wks,contour1,  \
         (/lon2(j)-lx,lon2(j)-lx,lon2(j)+lx,lon2(j)+lx,lon2(j)-lx/),  \
         (/lat2(j)-ly,lat2(j)+ly,lat2(j)+ly,lat2(j)-ly,lat2(j)-ly/),  \
                         plres)

 plres3 = True
 plres3@gsLineColor      = colind
 plres3@gsLineThicknessF = 1.8
 dumm3 = gsn_add_polyline(wks,contour3,  \
         (/maxlon(j)-dl(j,0),maxlon(j)+dl(j,0)/),  \
         (/maxlat(j)-dl(j,1),maxlat(j)+dl(j,1)/),plres3)

 plres3t = True
 plres3t@txFont        = 21
 plres3t@txFontHeightF = 0.008
 plres3t@txFontColor   = 1
 dumm3t1 = gsn_add_text(wks,contour3,"A",  \
           maxlon(j)-dl(j,0)-0.7,maxlat(j)-dl(j,1)+0.5,plres3t)
 dumm3t2 = gsn_add_text(wks,contour3,"B",  \
           maxlon(j)+dl(j,0)+0.7,maxlat(j)+dl(j,1)-0.5,plres3t)
 if (j .eq. 3) then
   dumm3t1b = gsn_add_text(wks,contour3,"A",  \
              maxlon(j)-dl(j,0)+0.7,maxlat(j)-dl(j,1)-0.3,plres3t)
 end if

 plres3d = True
 plres3d@gsMarkerIndex = 1
 plres3d@gsMarkerColor = colind
 plres3d@gsMarkerSizeF = 0.025
 dumm3d = gsn_add_polymarker(wks,contour3,  \
         maxlon(j),maxlat(j),plres3d)
 plres3d2 = True
 plres3d2@gsMarkerIndex = 1
 plres3d2@gsMarkerColor = 0
 plres3d2@gsMarkerSizeF = 0.010
 dumm3d2 = gsn_add_polymarker(wks,contour3,  \
         maxlon(j),maxlat(j),plres3d2)

 cnfinalize(contour1,wks)
 draw(contour1)
 delete(contour1)

 cnfinalize(contour3,wks)
 draw(contour3)
 delete(contour3)

 cnfinalize(contour2,wks)
 draw(contour2)
 delete(contour2)

 delete(dumm)
 delete(dumm3)

 dr_txt(wks,mstr(j),xf0+0.5*wi0-1,yf0(j)+he0-1.5,22,4.5,0.)

; xf = xf + (wi+10)

end do  ; j

 dr_txt(wks,"A",xf,yf+he+3,21,6.0,0.)
 dr_txt(wks,"B",xf+wi,yf+he+3,21,6.0,0.)

frame(wks)

end

