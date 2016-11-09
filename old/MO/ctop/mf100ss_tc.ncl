load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+12*12
 lat1 = -5
 lat2 = 5
 kz = 47

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

 wks = cr_wks("mf100ss_tc-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   nz = kz+2

   fnam1 = "$PORTAL/MO/dat/tr_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   lat = f1->latitude({lat1:lat2})
   z   = f1->hybrid_ht  (0:nz-1)/1.e3
   zr  = f1->hybrid_ht_1(0:nz-1)/1.e3
   pt1   = f1->theta      (4:nmon-1,:nz-1,{lat1:lat2},0)
   rhor1 = f1->unspecified(4:nmon-1,:nz-1,{lat1:lat2},0)/r2_earth

   ny = dimsizes(lat)
;   nz = dimsizes(z)

   fnam1 = "$PORTAL/MO/dat/u_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   ur1 = f1->u(4:nmon-1,:nz-1,{lat1:lat2},0)

   fnam1 = "$PORTAL/MO/dat/ng_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   mfeun1 = dim_avg_n( f1->field420(4:nmon-1,kz,{lat1:lat2},0), 1 )
   mfwun1 = dim_avg_n( f1->field422(4:nmon-1,kz,{lat1:lat2},0), 1 )

   fnam0 = "$PORTAL/MO/dat/spec_"+expname(0)+".nc"

   f0  = addfile(fnam0,"r")
   mfs0e = f0->unspecified  (4:nmon-1,kz,{lat1:lat2},:)/dc*norm
   mfs0w = f0->unspecified_1(4:nmon-1,kz,{lat1:lat2},:)/dc*norm*(-1.)

   fnam3 = "$PORTAL/MO/dat/tr_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   pt3   = f3->theta      (4:nmon-1,:nz-1,{lat1:lat2},0)
   rhor3 = f3->unspecified(4:nmon-1,:nz-1,{lat1:lat2},0)/r2_earth

   fnam3 = "$PORTAL/MO/dat/u_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   ur3 = f3->u(4:nmon-1,:nz-1,{lat1:lat2},0)

   fnam3 = "$PORTAL/MO/dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfeun2 = dim_avg_n( f3->field420(4:nmon-1,kz,{lat1:lat2},0), 1 )
   mfwun2 = dim_avg_n( f3->field422(4:nmon-1,kz,{lat1:lat2},0), 1 )

;--------------------------------------- process 0
 nt = nmon-4

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
 u3 = new((/nt,nz,ny/),float)
 do k=0, nz-2
   u3(:,k,:) = ur3(:,k,:)*az1(k)+ur3(:,k+1,:)*az0(k)
 end do
 u3(:,nz-1,:) = u3(:,nz-2,:)

 rho3 = new((/nt,nz,ny/),float)
 do k=0, nz-2
   rho3(:,k,:) = rhor3(:,k,:)*az1(k)+rhor3(:,k+1,:)*az0(k)
 end do
 rho3(:,nz-1,:) = rho3(:,nz-2,:)*exp(-dzl/h_scale)

 nbv3 = new((/nt,nz,ny/),float)
 do k=1, nz-2
   nbv3(:,k,:) = g/pt3(:,k,:)*(pt3(:,k+1,:)-pt3(:,k-1,:))/dz2m(k)
 end do
 nbv3(:,0,:) = nbv3(:,1,:)
 nbv3(:,nz-1,:) = nbv3(:,nz-2,:)
 do n=0, nt-1
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
 mfs3_l = new((/2,nt,ny,nc*2+1/),float)
 mfs3_l = 0.
 do n=0, nt-1
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

 mfs3 = new((/2,nt,ny,nz,nc*2+1/),float)
 mfs3 = 0.

 mfs3(:,:,:,kl,:) = mfs3_l(:,:,:,:)
 do n=0, nt-1
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

 t = fspan(0,13,14)

 mfs_avg = new((/3,2,nc*2+1,14/),float)
 do n=1, 12
   mfs_avg(0,0,:,n) = dim_avg_n(mfs0e(ispan(n-1,nt-1,12),:,:),(/0,1/))
   mfs_avg(0,1,:,n) = dim_avg_n(mfs0w(ispan(n-1,nt-1,12),:,:),(/0,1/))
   mfs_avg(1,:,:,n) = dim_avg_n(mfs1(:,ispan(n-1,nt-1,12),:,kz,:),(/1,2/))
   mfs_avg(2,:,:,n) = dim_avg_n(mfs3(:,ispan(n-1,nt-1,12),:,kz,:),(/1,2/))
 end do
 mfs_avg = mfs_avg*1.e3

  mfun = new((/2,2,12/),float)
  do n=0, 11
    mfun(0,0,n) = avg(mfeun1(ispan(n,nt-1,12)))*1.e3
    mfun(0,1,n) = avg(mfwun1(ispan(n,nt-1,12)))*1.e3
    mfun(1,0,n) = avg(mfeun2(ispan(n,nt-1,12)))*1.e3
    mfun(1,1,n) = avg(mfwun2(ispan(n,nt-1,12)))*1.e3
  end do
  cor1 = dim_sum_n(mfs_avg(1,:,:,1:12),1)*dc/mfun(0,:,0:11)
  cor2 = dim_sum_n(mfs_avg(2,:,:,1:12),1)*dc/mfun(1,:,0:11)
  print(cor2)

 do i=0, nc*2
   mfs_avg(1,:,i,1:12) = mfs_avg(1,:,i,1:12)/cor1(:,0:11)
   mfs_avg(2,:,i,1:12) = mfs_avg(2,:,i,1:12)/cor2(:,0:11)
 end do

 mfs_avg(1,:,:,:) = mfs_avg(0,:,:,:) + mfs_avg(1,:,:,:)

 mfs_avg(:,:,:,0 ) = mfs_avg(:,:,:,12)
 mfs_avg(:,:,:,13) = mfs_avg(:,:,:,1 )


 mfs_avg = mfs_avg - 1.e-6

;--------------------------------------- cn level
 cnlev = (/0,5,10,20,40,60,80,100,200,400,600,800,1000,2000/)*0.01
 cnfil = (/0,11,9 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22  ,21  ,19/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1/)*38
 yf = 10 + (/0,1,2/)*30

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-MF","Westward-MF"/)
 xstr = (/"","","month"/)
 ystr = (/"c:B:x:N: (m s:S:-1:N:)",""/)
 txtstr = (/"CGW", "CGW + NGW", "NGW (CTL)"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(10:S:-2:N: mPa / m s:S:-1:N:)",""/)

do iv=0, 2
do idr=0, 1

 contour1 = cnshadeplot(wks,t,c,mfs_avg(iv,idr,:,:)*1.e2,cnlev,cnfil,lbopt(iv,idr))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.05
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.3
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,-999,999,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
;  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 0
;  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

  dr_txt2(wks,txtstr(iv),xf(0)+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do
end do

  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,8,0.)
  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end



