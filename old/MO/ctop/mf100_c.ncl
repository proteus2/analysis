load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+12*12
 lat1 = -5
 lat2 = 5
 kz = 56  ; 47

 pi = acos(-1.)
 kl = 21
 mstar = 2.*pi/4.3e3
 a0 = 1.023659
 cl0 = 3.4191e-9 *1.5
 c_eq = 0.5
 beta = 0.1047
 m_min = 2.*pi/20.e3
 beq_ms = 2.3e-11/mstar
 mstar4 = mstar^4

 nc = 30
 cmax = 60.
 dc0 = cmax/nc
 dc = dc0*sqrt(2.)
 norm = 1./sqrt(2.)  ; yh. for old output

 r_earth = 6371.229
 omega = 7.292116e-5
 g = 9.80665
 h_scale = 287.05*239.145/g/1.e3

 cc_head = beta/pi*sqrt(2.)
 pm1 = 5./3. - 1.
 c2mp = 2. - 5./3.

 wks = cr_wks("fig/mf100_c-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/tr_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   lat = f1->latitude({lat1:lat2})
   z   = f1->hybrid_ht  (kl:kz+1)/1.e3
   zr  = f1->hybrid_ht_1(kl:kz+1)/1.e3
   pt1   = f1->theta      (4:nmon-1,kl-1:kz+1,{lat1:lat2},0)
   rhor1 = f1->unspecified(4:nmon-1,kl  :kz+1,{lat1:lat2},0)

   ny = dimsizes(lat)

   fnam1 = "../dat/u_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   ur1 = f1->u(4:nmon-1,kl:kz+1,{lat1:lat2},0)

   fnam1 = "../dat/ng_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   mfeun1 = f1->field420(4:nmon-1,kl:kz,{lat1:lat2},0)
   mfwun1 = f1->field422(4:nmon-1,kl:kz,{lat1:lat2},0)

   fnam0 = "../dat/spec_"+expname(0)+".nc"

   f0  = addfile(fnam0,"r")
   mfs0e = f0->unspecified  (4:nmon-1,kl:kz,{lat1:lat2},:)/dc*norm
   mfs0w = f0->unspecified_1(4:nmon-1,kl:kz,{lat1:lat2},:)/dc*norm*(-1.)
   mfs0e(:,0:1,:,:) = f0->unspecified  (4:nmon-1,0:1,{lat1:lat2},:)/dc*norm
   mfs0w(:,0:1,:,:) = f0->unspecified_1(4:nmon-1,0:1,{lat1:lat2},:)/dc*norm*(-1.)

   fnam3 = "../dat/tr_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   pt3   = f3->theta      (4:nmon-1,kl-1:kz+1,{lat1:lat2},0)
   rhor3 = f3->unspecified(4:nmon-1,kl  :kz+1,{lat1:lat2},0)

   fnam3 = "../dat/u_"+expname(1)+".nc"
   f3 = addfile(fnam3,"r")
   ur3 = f3->u(4:nmon-1,kl:kz+1,{lat1:lat2},0)

   fnam3 = "../dat/ng_"+expname(1)+".nc"

   f3 = addfile(fnam3,"r")
   mfeun2 = f3->field420(4:nmon-1,kl:kz,{lat1:lat2},0)
   mfwun2 = f3->field422(4:nmon-1,kl:kz,{lat1:lat2},0)

   kz = kz - kl
   nz = kz + 1

   do k=0, nz
     r2 = (r_earth+zr(k))^2*1.e6
     rhor1(:,k,:) = rhor1(:,k,:)/r2
     rhor3(:,k,:) = rhor3(:,k,:)/r2
   end do

;--------------------------------------- process 0
 nt = nmon-4

 dzr = new(nz,float)
 az0 = new(nz,float)
 az1 = new(nz,float)
 do k=0, nz-1
   dzr(k) = zr(k+1) - zr(k)
   az0(k) = z(k) - zr(k)
   az1(k) = zr(k+1) - z(k)
 end do
 az0 = az0/dzr
 az1 = az1/dzr
 dzl = z(nz) - z(nz-1)

 dz2m = new(nz,float)
 do k=1, nz-1
   dz2m(k) = (z(k+1) - z(k-1))*1.e3
 end do

; ------ u, N, rho for 1 -----
 u1   = new((/nt,nz,ny/),float)
 rho1 = new((/nt,nz,ny/),float)
 nbv1 = new((/nt,nz,ny/),float)
 do k=0, nz-1
   u1  (:,k,:) = ur1  (:,k,:)*az1(k)+ur1  (:,k+1,:)*az0(k)
   rho1(:,k,:) = rhor1(:,k,:)*az1(k)+rhor1(:,k+1,:)*az0(k)
   nbv1(:,k,:) = g/pt1(:,k+1,:)*(pt1(:,k+2,:)-pt1(:,k,:))/dz2m(k)
 end do
 do n=0, nt-1
 do k=0, nz-1
 do j=0, ny-1
   nbv1(n,k,j) = sqrt( max( (/1.e-4,nbv1(n,k,j)/) ) )
 end do
 end do
 end do

; ------ u, N, rho for 3 -----
 u3   = new((/nt,nz,ny/),float)
 rho3 = new((/nt,nz,ny/),float)
 nbv3 = new((/nt,nz,ny/),float)
 do k=0, nz-1
   u3  (:,k,:) = ur3  (:,k,:)*az1(k)+ur3  (:,k+1,:)*az0(k)
   rho3(:,k,:) = rhor3(:,k,:)*az1(k)+rhor3(:,k+1,:)*az0(k)
   nbv3(:,k,:) = g/pt3(:,k+1,:)*(pt3(:,k+2,:)-pt3(:,k,:))/dz2m(k)
 end do
 do n=0, nt-1
 do k=0, nz-1
 do j=0, ny-1
   nbv3(n,k,j) = sqrt( max( (/1.e-4,nbv3(n,k,j)/) ) )
 end do
 end do
 end do


 fcor2 = new(ny,float)
 fcor2(:) = (2.*omega*sin(lat(:)/180.*pi))^2

 cl0_1 = ( (tanh((lat-60.)/10.)-tanh((lat+60.)/10.))* \
          (1.0-c_eq)*0.5 + 1. )*cl0
 cl0_3 = cl0_1*0. + cl0
 cl0_3d_1 = new((/nt,nz,ny/),float)
 cl0_3d_3 = new((/nt,nz,ny/),float)
 do j=0, ny-1
   cl0_3d_1(:,:,j) = cl0_1(j)
   cl0_3d_3(:,:,j) = cl0_3(j)
 end do

 c = fspan(-nc*dc,nc*dc,nc*2+1)

 sig = new(nz,float)
 sig = 3.
 coef = new((/nz,31/),float)
 do k=0, nz-1
   coef(k,:) = exp(-((fspan(-15*dc,15*dc,31)/sig(k))^2))
   coef(k,:) = coef(k,:)/sum(coef(k,:))
 end do

; ------ mfs for 1 -----
 mfs1_l = new((/2,nt,ny,nc*2+1/),float)
 mfs1_l = 0.

 c_int = new((/nt,nz,ny,nc*2+1/),float)
 mwn   = new((/nt,nz,ny,nc*2+1/),float)
 mwn0  = new((/nt,nz,ny,nc*2+1,31/),float)
 do i=0, nc*2
   c_int(:,:,:,i) = c(i) - u1(:,:,:)
   mwn  (:,:,:,i) = nbv1(:,:,:)/abs(c_int(:,:,:,i))
   do ii=0, 30
     mwn0(:,:,:,i,ii) = nbv1(:,:,:)/abs(c_int(:,:,:,i)+(ii-15)*dc)
   end do
 end do
 tmp = rho1(:,:,:)*a0*cl0_3d_1(:,:,:)/nbv1(:,:,:)

 do n=0, nt-1
 do j=0, ny-1
 do i=0, nc*2
   if (c_int(n,0,j,i) .gt. 0.) then
     iv = 0
   else
     iv = 1
   end if
   if ( mwn(n,0,j,i) .ge. m_min .and. mwn(n,0,j,i) .lt. mstar ) then
     mfs1_l(iv,n,j,i) = tmp(n,0,j)*mwn(n,0,j,i)^3/mstar4
   else
     if ( mwn(n,0,j,i) .ge. mstar ) then
       mfs1_l(iv,n,j,i) = tmp(n,0,j)/mwn(n,0,j,i)
     end if
   end if
 end do
 end do
 end do

 mfs1 = new((/2,nt,ny,nz,nc*2+1/),float)
 mfs1 = 0.

 omin = new((/nt,ny/),float)
 do n=0, nt-1
 do j=0, ny-1
   omin(n,j) = sqrt( max( (/fcor2(j), nbv1(n,0,j)*beq_ms/) ) )
;   omin = sqrt( max( (/fcor2(j), nbv1(n,k,j)*beq_ms/) ) )
 end do
 end do
 omin_nbv = new((/nt,nz,ny/),float) 
 do k=0, nz-1
   omin_nbv(:,k,:) = omin(:,:)/nbv1(:,k,:)
 end do
 cc = cc_head*nbv1(:,:,:)^2*omin_nbv(:,:,:)^pm1*(1.-omin_nbv(:,:,:)^c2mp) \
        /(1.-omin_nbv(:,:,:)^pm1)
 tmp = rho1(:,:,:)*a0*cc(:,:,:)/nbv1(:,:,:)
 mfs1s = new((/nt,nz,ny,nc*2+1/),float)
 mfs1s0 = new((/nt,nz,ny,nc*2+1,31/),float)
 do i=0, nc*2
   mfs1s(:,:,:,i) = tmp(:,:,:)/mwn(:,:,:,i)
   do ii=0, 30
     mfs1s0(:,:,:,i,ii) = tmp(:,:,:)/mwn0(:,:,:,i,ii)
   end do
 end do

 mfs1(:,:,:,0,:) = (/ mfs1_l(:,:,:,:) /)
 do k=1, nz-1
 do ii=0, 30
   lcri = ( ((c_int(:,k,:,:)+(ii-15)*dc)*c_int(:,k-1,:,:)) / \
         abs((c_int(:,k,:,:)+(ii-15)*dc)*c_int(:,k-1,:,:)) + 1.)*0.5
   do iv=0, 1
     mfs1(iv,:,:,k,:) = mfs1(iv,:,:,k,:) + coef(k,ii)* \
                      dim_min_n( (/mfs1s0(:,k,:,:,ii), mfs1(iv,:,:,k-1,:)/), 0 )* \
                      lcri(:,:,:)
   end do
 end do
 end do

; ------ mfs for 3 -----
 mfs3_l = new((/2,nt,ny,nc*2+1/),float)
 mfs3_l = 0.

 c_int = new((/nt,nz,ny,nc*2+1/),float)
 mwn   = new((/nt,nz,ny,nc*2+1/),float)
 mwn0  = new((/nt,nz,ny,nc*2+1,31/),float)
 do i=0, nc*2
   c_int(:,:,:,i) = c(i) - u3(:,:,:)
   mwn  (:,:,:,i) = nbv3(:,:,:)/abs(c_int(:,:,:,i))
   do ii=0, 30
     mwn0(:,:,:,i,ii) = nbv3(:,:,:)/abs(c_int(:,:,:,i)+(ii-15)*dc)
   end do
 end do
 tmp = rho3(:,:,:)*a0*cl0_3d_3(:,:,:)/nbv3(:,:,:)

 do n=0, nt-1
 do j=0, ny-1
 do i=0, nc*2
   if (c_int(n,0,j,i) .gt. 0.) then
     iv = 0
   else
     iv = 1
   end if
   if ( mwn(n,0,j,i) .ge. m_min .and. mwn(n,0,j,i) .lt. mstar ) then
     mfs3_l(iv,n,j,i) = tmp(n,0,j)*mwn(n,0,j,i)^3/mstar4
   else
     if ( mwn(n,0,j,i) .ge. mstar ) then
       mfs3_l(iv,n,j,i) = tmp(n,0,j)/mwn(n,0,j,i)
     end if
   end if
 end do
 end do
 end do

 mfs3 = new((/2,nt,ny,nz,nc*2+1/),float)
 mfs3 = 0.

 omin = new((/nt,ny/),float)
 do n=0, nt-1
 do j=0, ny-1
   omin(n,j) = sqrt( max( (/fcor2(j), nbv3(n,0,j)*beq_ms/) ) )
;   omin = sqrt( max( (/fcor2(j), nbv3(n,k,j)*beq_ms/) ) )
 end do
 end do
 omin_nbv = new((/nt,nz,ny/),float)
 do k=0, nz-1
   omin_nbv(:,k,:) = omin(:,:)/nbv3(:,k,:)
 end do
 cc = cc_head*nbv3(:,:,:)^2*omin_nbv(:,:,:)^pm1*(1.-omin_nbv(:,:,:)^c2mp) \
        /(1.-omin_nbv(:,:,:)^pm1)
 tmp = rho3(:,:,:)*a0*cc(:,:,:)/nbv3(:,:,:)
 mfs3s = new((/nt,nz,ny,nc*2+1/),float)
 mfs3s0 = new((/nt,nz,ny,nc*2+1,31/),float)
 do i=0, nc*2
   mfs3s(:,:,:,i) = tmp(:,:,:)/mwn(:,:,:,i)
   do ii=0, 30
     mfs3s0(:,:,:,i,ii) = tmp(:,:,:)/mwn0(:,:,:,i,ii)
   end do
 end do

 mfs3(:,:,:,0,:) = (/ mfs3_l(:,:,:,:) /)
 do k=1, nz-1
 do ii=0, 30
   lcri = ( ((c_int(:,k,:,:)+(ii-15)*dc)*c_int(:,k-1,:,:)) / \
         abs((c_int(:,k,:,:)+(ii-15)*dc)*c_int(:,k-1,:,:)) + 1.)*0.5
   do iv=0, 1
     mfs3(iv,:,:,k,:) = mfs3(iv,:,:,k,:) + coef(k,ii)* \
                      dim_min_n( (/mfs3s0(:,k,:,:,ii), mfs3(iv,:,:,k-1,:)/), 0 )* \
                      lcri(:,:,:)
   end do
 end do
 end do


 mfs_avg_m = new((/3,2,nz,ny,nc*2+1,12/),float)
 do n=0, 11
   mfs_avg_m(0,0,:,:,:,n) = dim_avg_n(mfs0e(ispan(n,nt-1,12),:kz,:,:),0)
   mfs_avg_m(0,1,:,:,:,n) = dim_avg_n(mfs0w(ispan(n,nt-1,12),:kz,:,:),0)
   do k=0, kz
     mfs_avg_m(1,:,k,:,:,n) = dim_avg_n(mfs1(:,ispan(n,nt-1,12),:,k,:),1)
     mfs_avg_m(2,:,k,:,:,n) = dim_avg_n(mfs3(:,ispan(n,nt-1,12),:,k,:),1)
   end do
 end do
 mfs_avg_m = mfs_avg_m*1.e3

  mfun = new((/2,2,nz,ny,12/),float)
  do n=0, 11
    mfun(0,0,:,:,n) = dim_avg_n(mfeun1(ispan(n,nt-1,12),:kz,:),0)*1.e3
    mfun(0,1,:,:,n) = dim_avg_n(mfwun1(ispan(n,nt-1,12),:kz,:),0)*1.e3
    mfun(1,0,:,:,n) = dim_avg_n(mfeun2(ispan(n,nt-1,12),:kz,:),0)*1.e3
    mfun(1,1,:,:,n) = dim_avg_n(mfwun2(ispan(n,nt-1,12),:kz,:),0)*1.e3
  end do
  cor1 = dim_sum_n(mfs_avg_m(1,:,:,:,:,:),3)*dc/mfun(0,:,:,:,:)
  cor2 = dim_sum_n(mfs_avg_m(2,:,:,:,:,:),3)*dc/mfun(1,:,:,:,:)
;print(dim_sum_n(dim_avg_n(dim_sum_n(mfs_avg_m(0:1,:,1,:,:,:),3)*dc,(/2,3/)),0))
;print(dim_avg_n(dim_sum_n(mfs_avg_m(2,:,0,:,:,:),2)*dc,(/1,2/)))
print(dim_avg_n(cor2(:,:,:,:),(/2,3/)))

; cor1 = cor1 + 0.000001
; cor2 = cor2 + 0.000001
; do i=0, nc*2
;   mfs_avg_m(1,:,:,:,i,:) = mfs_avg_m(1,:,:,:,i,:)/cor1(:,:,:,:)
;   mfs_avg_m(2,:,:,:,i,:) = mfs_avg_m(2,:,:,:,i,:)/cor2(:,:,:,:)
; end do

 mfs_avg = dim_avg_n(dim_avg(mfs_avg_m),3)

 mfs_avg(1,:,:,:) = mfs_avg(0,:,:,:) + mfs_avg(1,:,:,:)

;--------------------------------------- cn level
 cnlev = (/0,5,10,20,40,60,80,100,200,400,600,800,1000,2000/)*0.01
 cnfil = (/0,11,9 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22  ,21  ,19/)

;--------------------------------------- contour
 wi = 32
 he = 24
 xf = 12 + (/0,1/)*38
 yf = 10 + (/0,1,2/)*30

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"100 hPa","20 hPa"/)
 xstr = "c:B:x:N: (m s:S:-1:N:)"
 ystr = (/"MF spec. (mPa / m s:S:-1:N:)",""/)
 xmax = (/60,40/)
 xmin = (/-60,-40/)

; kk = (/47-kl,kz/)
 kk = (/0,kz/)
 kp1 = 1

 mfs_avg = mfs_avg + 1.e-18

 plot = xylineplot(wks,c,mfs_avg(2,0:1,kk(0),:),"D",(/0,16/),"T",0.7,"C",1)
 xyaddline(plot,c,mfs_avg(1,0:1,kk(0)+kp1,:),"D",(/0,16/),"T",0.7,"C",2)
 xyaddline(plot,c,mfs_avg(2,0:1,kk(1),:),"D",(/0,16/),"T",1.6,"C",1)
 xyaddline(plot,c,mfs_avg(1,0:1,kk(1),:),"D",(/0,16/),"T",1.6,"C",2)

 set_size(plot,xf(0),yf(0),wi,he)
 axis_range(plot,xmin(0),xmax(0),1.e-3,3.e-1,False,True)
 axis_str(plot,"",xstr,ystr(0))

 draw(plot)

 frame(wks)

 mfs_avg = mfs_avg - 1.e-18


 wks2 = cr_wks("fig/mf100_c2-"+expname(0),600,0)
; gsn_define_colormap(wks2,"StepSeq25")
 colind = NhlNewColor(wks2,0.4,0.4,0.4)

 plot2 = xylineplot(wks2,c,mfs_avg(2,0:1,kk(0),:),"D",(/0,16/),"T",0.7,"C",1)
 xyaddline(plot2,c,mfs_avg(1,0:1,kk(0)+kp1,:),"D",(/0,16/),"T",0.7,"C",2)
 xyaddline(plot2,c,mfs_avg(2,0:1,kk(1),:),"D",(/0,16/),"T",1.6,"C",1)
 xyaddline(plot2,c,mfs_avg(1,0:1,kk(1),:),"D",(/0,16/),"T",1.6,"C",2)

 set_size(plot2,xf(0),yf(0),wi,he)
 axis_range(plot2,xmin(1),xmax(1),0,999,False,False)
 axis_str(plot2,"",xstr,ystr(0))

 draw(plot2)

 frame(wks2)


end



