load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 im = 4+ispan(0, 143, 1)
 im3= im
 lat0 = 15 ;20
 kz = 68

 nmon = dimsizes(im)
 nmon3= dimsizes(im3)

 pi = acos(-1.)
 kl = 21
 mstar = 2.*pi/4.3e3
 a0 = 1.023659
 cl0 = 3.4191e-9 *1.5
 c_eq = 0.5
 beta = 0.1047
 m_min = 2.*pi/20.e3
 beq_ms = 2.3e-11/mstar

 nc = 30
 cmax = 60.
 dc0 = cmax/nc
 dc = dc0*sqrt(2.)
 norm = 1./sqrt(2.)  ; yh. for old output

 r2_earth = 6371229.^2
 omega = 7.292116e-5
 g = 9.80665
 h_scale = 287.05*239.145/g/1.e3

 cc_head = beta/pi*sqrt(2.)
 pm1 = 5./3. - 1.
 c2mp = 2. - 5./3.

 wks = cr_wks("z_fc-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/tr_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   lat = f1->latitude({-lat0:lat0})
   z   = f1->hybrid_ht  /1.e3
   zr  = f1->hybrid_ht_1/1.e3
   pt1   = f1->theta      (im,:,{-lat0:lat0},0)
   rhor1 = f1->unspecified(im,:,{-lat0:lat0},0)/r2_earth

   ny = dimsizes(lat)
   nz = dimsizes(z)

   fnam1 = "../dat/u_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   ur1 = f1->u(im,:,{-lat0:lat0},0)

   fnam0 = "../dat/spec_"+expname(0)+".nc"

   f0  = addfile(fnam0,"r")
   mfs0e = f0->unspecified  (im,:,{-lat0:lat0},:)/dc*norm
   mfs0w = f0->unspecified_1(im,:,{-lat0:lat0},:)/dc*norm*(-1.)

   fnam3 = "../dat/tr_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   pt3   = f3->theta      (im3,:,{-lat0:lat0},0)
   rhor3 = f3->unspecified(im3,:,{-lat0:lat0},0)/r2_earth

   fnam3 = "../dat/u_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   ur3 = f3->u(im3,:,{-lat0:lat0},0)

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfeun2 = dim_avg_n( dim_avg_n( f3->field420(im3,:,{-lat0:lat0},0), 2 ), 0 )
   mfwun2 = dim_avg_n( dim_avg_n( f3->field422(im3,:,{-lat0:lat0},0), 2 ), 0 )

;--------------------------------------- process 0
 nt = nmon
 nt3 = nmon3

 dzr = new(nz-1,float)
 az0 = new(nz-1,float)
 az1 = new(nz-1,float)
 do k=0, nz-2
   dzr(k) = zr(k+1) - zr(k)
   az0(k) = z(k) - zr(k)
   az1(k) = zr(k+1) - z(k)
 end do
 az0 = az0/dzr
 az1 = az1/dzr
 dzl = z(nz-1) - z(nz-2)

 dz2m = new(nz-1,float)
 do k=1, nz-2
   dz2m(k) = (z(k+1) - z(k-1))*1.e3
 end do

; ------ u, N, rho for 1 -----
 u1 = new((/nt,nz,ny/),float)
 do k=0, nz-2
   u1(:,k,:) = ur1(:,k,:)*az1(k)+ur1(:,k+1,:)*az0(k)
 end do
 u1(:,nz-1,:) = u1(:,nz-2,:)

 rho1 = new((/nt,nz,ny/),float)
 do k=0, nz-2
   rho1(:,k,:) = rhor1(:,k,:)*az1(k)+rhor1(:,k+1,:)*az0(k)
 end do
 rho1(:,nz-1,:) = rho1(:,nz-2,:)*exp(-dzl/h_scale)

 nbv1 = new((/nt,nz,ny/),float)
 do k=1, nz-2
   nbv1(:,k,:) = g/pt1(:,k,:)*(pt1(:,k+1,:)-pt1(:,k-1,:))/dz2m(k)
 end do
 nbv1(:,0,:) = nbv1(:,1,:)
 nbv1(:,nz-1,:) = nbv1(:,nz-2,:)
 do n=0, nt-1
 do k=0, nz-1
 do j=0, ny-1
   nbv1(n,k,j) = sqrt( max( (/1.e-4,nbv1(n,k,j)/) ) )
 end do
 end do
 end do

; ------ u, N, rho for 3 -----
 u3 = new((/nt3,nz,ny/),float)
 do k=0, nz-2
   u3(:,k,:) = ur3(:,k,:)*az1(k)+ur3(:,k+1,:)*az0(k)
 end do
 u3(:,nz-1,:) = u3(:,nz-2,:)

 rho3 = new((/nt3,nz,ny/),float)
 do k=0, nz-2
   rho3(:,k,:) = rhor3(:,k,:)*az1(k)+rhor3(:,k+1,:)*az0(k)
 end do
 rho3(:,nz-1,:) = rho3(:,nz-2,:)*exp(-dzl/h_scale)

 nbv3 = new((/nt3,nz,ny/),float)
 do k=1, nz-2
   nbv3(:,k,:) = g/pt3(:,k,:)*(pt3(:,k+1,:)-pt3(:,k-1,:))/dz2m(k)
 end do
 nbv3(:,0,:) = nbv3(:,1,:)
 nbv3(:,nz-1,:) = nbv3(:,nz-2,:)
 do n=0, nt3-1
 do k=0, nz-1
 do j=0, ny-1
   nbv3(n,k,j) = sqrt( max( (/1.e-4,nbv3(n,k,j)/) ) )
 end do
 end do
 end do


; ------ mfs for 1 -----
 fcor2 = new(ny,float)
 fcor2(:) = (2.*omega*sin(lat(:)/180.*pi))^2

 cl0y = ( (tanh((lat-60.)/10.)-tanh((lat+60.)/10.))* \
          (1.0-c_eq)*0.5 + 1. )*cl0

 c = fspan(-nc*dc,nc*dc,nc*2+1)
 mfs1_l = new((/2,nt,ny,nc*2+1/),float)
 mfs1_l = 0.
 do n=0, nt-1
 do j=0, ny-1
 do i=0, nc*2
   c_int = c(i) - u1(n,kl,j)
   if (c_int .ne. 0.) then
     mwn = nbv1(n,kl,j)/abs(c_int)
     tmp = rho1(n,kl,j)*a0*cl0y(j)/nbv1(n,kl,j)
     if (c_int .gt. 0.) then
       iv = 0
     else
       iv = 1
     end if
     if ( mwn .ge. m_min .and. mwn .lt. mstar ) then
       mfs1_l(iv,n,j,i) = tmp*mwn^3/mstar^4
     else
       if ( mwn .ge. mstar ) then
         mfs1_l(iv,n,j,i) = tmp/mwn
       end if
     end if
   end if
 end do
 end do
 end do

 mfs1 = new((/2,nt,ny,nz,nc*2+1/),float)
 mfs1 = 0.

 mfs1(:,:,:,kl,:) = mfs1_l(:,:,:,:)
 do n=0, nt-1
 do j=0, ny-1
   omin = sqrt( max( (/fcor2(j), nbv1(n,kl,j)*beq_ms/) ) )
   do k=kl+1, nz-1
   do i=0, nc*2
     c_int = c(i) - u1(n,k,j)
     if ( ( mfs1(0,n,j,k-1,i) .ne. 0. .or. mfs1(1,n,j,k-1,i) .ne. 0. ) \
          .and. c_int .ne. 0. ) then
       mwn = nbv1(n,k,j)/abs(c_int)
;       omin = sqrt( max( (/fcor2(j), nbv1(n,k,j)*beq_ms/) ) )
       omin_nbv = omin/nbv1(n,k,j)
       cc = cc_head*nbv1(n,k,j)^2*omin_nbv^pm1*(1.-omin_nbv^c2mp) \
            /(1.-omin_nbv^pm1)
       mfs1s = rho1(n,k,j)*a0*cc/nbv1(n,k,j)/mwn
       if (c_int .gt. 0.) then
         mfs1(0,n,j,k,i) = min( (/mfs1s, mfs1(0,n,j,k-1,i)/) )
       else
         mfs1(1,n,j,k,i) = min( (/mfs1s, mfs1(1,n,j,k-1,i)/) )
       end if
     end if
   end do
   end do
 end do
 end do

; ------ mfs for 3 -----
 c = fspan(-nc*dc,nc*dc,nc*2+1)
 mfs3_l = new((/2,nt3,ny,nc*2+1/),float)
 mfs3_l = 0.
 do n=0, nt3-1
 do j=0, ny-1
 do i=0, nc*2
   c_int = c(i) - u3(n,kl,j)
   if (c_int .ne. 0.) then
     mwn = nbv3(n,kl,j)/abs(c_int)
     tmp = rho3(n,kl,j)*a0*cl0/nbv3(n,kl,j)
     if (c_int .gt. 0.) then
       iv = 0
     else
       iv = 1
     end if
     if ( mwn .ge. m_min .and. mwn .lt. mstar ) then
       mfs3_l(iv,n,j,i) = tmp*mwn^3/mstar^4
     else
       if ( mwn .ge. mstar ) then
         mfs3_l(iv,n,j,i) = tmp/mwn
       end if
     end if
   end if
 end do
 end do
 end do

 mfs3 = new((/2,nt3,ny,nz,nc*2+1/),float)
 mfs3 = 0.

 mfs3(:,:,:,kl,:) = mfs3_l(:,:,:,:)
 do n=0, nt3-1
 do j=0, ny-1
   omin = sqrt( max( (/fcor2(j), nbv3(n,kl,j)*beq_ms/) ) )
   do k=kl+1, nz-1
   do i=0, nc*2
     c_int = c(i) - u3(n,k,j)
     if ( ( mfs3(0,n,j,k-1,i) .ne. 0. .or. mfs3(1,n,j,k-1,i) .ne. 0. ) \
          .and. c_int .ne. 0. ) then
       mwn = nbv3(n,k,j)/abs(c_int)
;       omin = sqrt( max( (/fcor2(j), nbv3(n,k,j)*beq_ms/) ) )
       omin_nbv = omin/nbv3(n,k,j)
       cc = cc_head*nbv3(n,k,j)^2*omin_nbv^pm1*(1.-omin_nbv^c2mp) \
            /(1.-omin_nbv^pm1)
       mfs3s = rho3(n,k,j)*a0*cc/nbv3(n,k,j)/mwn
       if (c_int .gt. 0.) then
         mfs3(0,n,j,k,i) = min( (/mfs3s, mfs3(0,n,j,k-1,i)/) )
       else
         mfs3(1,n,j,k,i) = min( (/mfs3s, mfs3(1,n,j,k-1,i)/) )
       end if
     end if
   end do
   end do
 end do
 end do


 mfs_avg = new((/2,4,nt,nc*2+1/),float)
 mfs_avg(0,0,:,:) = dim_avg_n(mfs0e(:,kz,:,:),1)*1.e3
 mfs_avg(1,0,:,:) = dim_avg_n(mfs0w(:,kz,:,:),1)*1.e3
 mfs_avg(:,1,:,:) = dim_avg_n(mfs1(:,:,:,kz,:),2)*1.e3
 mfs_avg(:,2,:,:) = mfs_avg(:,0,:,:) + mfs_avg(:,1,:,:)
 mfs_avg(:,3,:,:) = dim_avg_n(mfs3(:,:,:,kz,:),2)*1.e3

;print( (dim_sum(mfs_avg(0,3,21:68,:))*dc/1.e3/mfeun2(21:68))+"   "+ \
;        dim_avg_n(dim_avg(nbv3(:,21:68,:)),0) )
; check mfs_avg(0,1,:,:) also

 dzm = new(nz,float)
 do k=1, nz-1
   dzm(k) = (z(k) - z(k-1))*1.e3
 end do

 wds0e = new(dimsizes(mfs0e),float)
 wds0w = new(dimsizes(mfs0w),float)
 wds1  = new(dimsizes(mfs1 ),float)
 wds3  = new(dimsizes(mfs3 ),float)
 wds0e = 0.
 wds0w = 0.
 wds1  = 0.
 wds3  = 0.
 do k=3, nz-1
 do i=0, nc*2
   wds0e(:,k,:,i) = (mfs0e(:,k-1,:,i)-mfs0e(:,k,:,i))/dzm(k)/rhor1(:,k,:)
   wds0w(:,k,:,i) = (mfs0w(:,k-1,:,i)-mfs0w(:,k,:,i))/dzm(k)/rhor1(:,k,:)
 end do
 end do
 do iv=0, 1
 do k=1, nz-1
 do i=0, nc*2
   wds1(iv,:,:,k,i) = (mfs1(iv,:,:,k-1,i)-mfs1(iv,:,:,k,i))/dzm(k)/rhor1(:,k,:)
   wds3(iv,:,:,k,i) = (mfs3(iv,:,:,k-1,i)-mfs3(iv,:,:,k,i))/dzm(k)/rhor3(:,k,:)
 end do
 end do
 end do

 wds_avg = new((/2,4,nt,nc*2+1/),float)
 wds_avg(0,0,:,:) = dim_avg_n(wds0e(:,kz,:,:),1)*86400.
 wds_avg(1,0,:,:) = dim_avg_n(wds0w(:,kz,:,:),1)*86400.
 wds_avg(:,1,:,:) = dim_avg_n(wds1(:,:,:,kz,:),2)*86400.
 wds_avg(:,2,:,:) = wds_avg(:,0,:,:) + wds_avg(:,1,:,:)
 wds_avg(:,3,:,:) = dim_avg_n(wds3(:,:,:,kz,:),2)*86400.
 wds_avg(1,:,:,:) = wds_avg(1,:,:,:)*(-1.)

 wds = new((/4,nc*2+1,nt/),float)
 do i=0, nc*2
   wds(:,i,:) = wds_avg(0,:,:,i) + wds_avg(1,:,:,i)
 end do

 ss = new((/4,nc*2+1,12/),float)
 do n=0, 11
   ss(:,:,n) = dim_avg(wds(:,:,ispan(n,nt-1,12)))
 end do

 do n=0, nt-1
   wds(:,:,n) = wds(:,:,n) - ss(:,:,mod(n,12))
 end do

   nf  = nt/2
   fr  = fspan(1,nf,nf)/(int2flt(nt)/12.)

   fftu  = ezfftf(wds)

   pws  = ( fftu(0,:,:,:)^2 + fftu(1,:,:,:)^2 )*0.5*nt/12. 

 mfs_avg = mfs_avg - 1.e-6

;--------------------------------------- cn level
 cnlev2 = (/1,2,5,10,20,50,100,200/)*0.01
 cnfil2 = (/0,11,9,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 32
 he = 42
 xf = 12 + (/0,1,2,3/)*40
 yf = 10 + (/0,1,2/)*56 + (/0,0,20/)

 mstr = (/ (/"CGW MF","NGW MF","Sum","NGW MF (CTL)"/), (/"","","",""/), \
           (/"CGW forcing","NGW forcing","Sum","NGW forcing (CTL)"/) /)
 lbopt = (/ (/"","","",""/), (/"","(10:S:-3:N: mPa / m s:S:-1:N:)","",""/), \
            (/"","(10:S:-3:N: (m s:S:-1:N: day:S:-1:N:):S:2:N: / m s:S:-1:N:)","",""/) /)
 ystr = (/"c:B:x:N: (m s:S:-1:N:)","","",""/)

 mstr2 = (/"CGW forcing","NGW forcing","Sum","NGW forcing (CTL)"/)
 lbopt2 = (/"","(10:S:-3:N: m s:S:-1:N: day:S:-1:N: / m s:S:-1:N:)","",""/)

do iv=0, 3

 contour2 = cnshadeplot(wks,fr,c,pws(iv,:,:)*1.e3,cnlev2,cnfil2,lbopt(2,iv))
 setvalues contour2
   "pmLabelBarOrthogonalPosF" : 0.06
   "pmLabelBarParallelPosF"   : 0.7
   "pmLabelBarKeepAspect"     : False
   "pmLabelBarWidthF"         : 2.2
   "pmLabelBarHeightF"        : 0.12
   "lbLabelAutoStride" : False
   "lbLabelFontHeightF" : 0.03
   "lbTitlePosition"    : "Right"
   "lbTitleDirection"   : "Across"
   "lbTitleOffsetF"     : 0.03
   "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour2,xf(iv),yf(2),wi,he)
 axis_range(contour2,-999,2.5,-999,999,False,False)
 axis_str(contour2,mstr(2,iv),"frequency (cycle yr:S:-1:N:)",ystr(iv))

 cnfinalize(contour2,wks)


 draw(contour2)

end do


 frame(wks)

end



