program cond_eof_analysis

  use eofpack
  use netcdfio

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxp0 = 45, nyp0 = 18, nxa0 = 80, nya0 = 108
  integer, parameter :: nxa = 27, nya = 36
  integer, parameter :: nxo = 320, nyo = 108
  integer, parameter :: nt = 1200, mlag = 240
  real,    parameter :: dt = 1./12.
  real,    parameter :: missvd = 1.e32, missvo = 1.e32
!--------------------------------------------------------------
  integer :: i,j,n, ncid, ii, jj, tag, missn, nn, k
  real    :: pi

  real, dimension(nxa,nya,nt) :: toa, toare

  real    :: lona(nxa), lata(nya)
  real    :: lono(nxo), lato(nyo), t1(nt), inino(nt), ininoo(nt)
  real    :: lona0(nxa0), lata0(nya0)
  real    :: normo(nyo), normoa(nxa,nya), corr(nxa,nya,mlag*2+1)
  real    :: lag(mlag*2+1)
  real    :: eof(nxa,nya,10)
  real    :: pc(nt,10), corrk(mlag*2+1,10), lagk(10), cork(10)

  real    :: avg1(nt), avg2(nt), tavg1, tavg2
  real    :: lag2(nxa,nya), cor2(nxa,nya)
  real, dimension(:,:,:), allocatable :: toa1, toa2, top0, toa0

  integer :: tmpi, lagg
  real    :: tmp1(nxa0), tmp2(nya0), tmp3(nxo), tmp4(nyo)
  real*8  :: tmp5(2,nya0), tmp6(2,nyo)
  real    :: tmpr1, tmpr2, tmpr3, summ, summ1, summ2

  character*128 :: ofname1, ofname2
  character*128 :: fdir(2)
  character*128 :: fn_to

  real :: reof(nxa,nya,10), pcsub(nt,10)


  pi = acos(-1.)

! read toarea ------------------------------------------------------------
  fdir(1) = '/export30/kyh/MIROC/'
  write(fn_to,'(a)') trim(fdir(1))//'tosa.nc'

! axis
  call opennc(fn_to,ncid)
  call get1d(ncid,'lon',nxo,tmp3)
  call get1d(ncid,'lat',nyo,tmp4)
  lono = real(tmp3)
  lato = real(tmp4)
  lona0(:61) = lono(260:) - 360.
  lona0(62:) = lono(1:19)
  lata0(:) = lato(:)

  do n=1, nt
    t1(n) = real(n-1)/12.
  enddo

! variables
  allocate(top0(nxp0,nyp0,nt))                        ; top0 = 0.
  allocate(toa0(nxa0,nya0,nt))                        ; toa0 = 0.
  allocate(toa1(61,nya0,nt))                          ; toa1 = 0.
  allocate(toa2(19,65,nt))                            ; toa2 = 0.
  call geta3d(ncid,'TOSA',170,nxp0,46,nyp0,1,nt,top0)
  call geta3d(ncid,'TOSA',260,61,1,nya0,1,nt,toa1)
  call geta3d(ncid,'TOSA',1,19,1,65,1,nt,toa2)
  call get1d(ncid,'norm',nyo,normo)
  call closenc(ncid)
  toa0 = 1.e32
  toa0(:61,:,:) = toa1(:,:,:)
  toa0(62:,:65,:) = toa2(:,:,:)
  deallocate(toa1)  ;  deallocate(toa2)

  missn = 0
  jj = 0 
  do j=2, nya0, 3
    jj = jj + 1 
    ii = 0
    do i=2, nxa0, 3 
      ii = ii + 1
      toa(ii,jj,:) = toa0(i,j,:)
      tag = 0
!      do n=1, nt 
        if (toa0(i,j,1) .eq. missvd)  tag = 1
!      enddo
      if (tag .eq. 1) then
        missn = missn + 1
!        toa(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(toa0)

  ii = 0
  do i=2, nxa0, 3
    ii = ii + 1
    lona(ii) = lona0(i)
  enddo
  jj = 0 
  do j=2, nya0, 3
    jj = jj + 1
    lata(jj) = lata0(j)
    normoa(:,jj) = normo(j)
  enddo

! (1) Nino3.4 index & SSTA over Atlantic ocean -------------------------
  ininoo = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyp0
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxp0
        if (top0(i,j,n) .ne. missvd) then
          tmpi = tmpi + 1 
          tmpr1 = tmpr1 + top0(i,j,n)
        end if
      enddo
      ininoo(n) = ininoo(n) + tmpr1 * normo(45+j)
      tmpr2 = tmpr2 + tmpi * normo(45+j)
    enddo
    ininoo(n) = ininoo(n) / tmpr2
  enddo

  call detrend_season(nt,ininoo)

  do n=2, nt-1
    inino(n) = (ininoo(n-1)+ininoo(n)+ininoo(n+1)) / 3.
  enddo
  inino(1) = ininoo(1)
  inino(nt)= ininoo(nt)


! (1) EOF analysis -----------------------------------------------------
if (0 .eq. 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lagg = 7 !>0
  do j=1, nya 
  do i=1, nxa 
  if (toa(i,j,1) .ne. missvd) then 
      tmpr1= 0. 
      tmpr2= 0. 
      do n=1, nt-lagg
        tmpr1 = tmpr1 + inino(n+lagg) 
        tmpr2 = tmpr2 + toa(i,j,n) 
      enddo 
      tmpr1 = tmpr1 / (nt-lagg) 
      tmpr2 = tmpr2 / (nt-lagg) 
      summ = 0. 
      summ1= 0. 
      summ2= 0. 
      do n=1, nt-lagg
        summ = summ + (inino(n+lagg)-tmpr1) * (toa(i,j,n)-tmpr2)
        summ1= summ1+ (inino(n+lagg)-tmpr1)**2
        summ2= summ2+ (toa(i,j,n)-tmpr2)**2
      enddo 
    do n=1, nt-lagg
      toa(i,j,n) = toa(i,j,n) - inino(n+lagg)*summ/summ1
    enddo 
 
  end if 
  enddo 
  enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

  ofname1 = '../eofatl'
  call eofnc(ofname1,nxa,lona,nya,lata,nt,t1,toa,missn,missvo,&
                  30,1,normoa)

  call opennc('../eofatl.nc',ncid)
  call geta3d(ncid,'EOF',1,nxa,1,nya,1,10,eof)
  call closenc(ncid)
  ofname1 = '../reof'
  call reofnc(ofname1,nxa,lona,nya,lata,nt,t1,toa,missn,missvo, &
              10,eof,50,1,normoa)

! (2) Conditional EOF analysis -----------------------------------------

! EOF vs NINO3.4 lag corr.
  call opennc('../reof.nc',ncid)
  call geta2d(ncid,'PC',1,nt,1,10,pc)
  call get3d(ncid,'REOF',nxa,nya,10,reof)
  call closenc(ncid)

do jj=1, 10

  do k=1, 10
    do nn=-mlag, 0
      tmpr1= 0.
      tmpr2= 0.
      do n=1, nt+nn
        tmpr1 = tmpr1 + inino(n)
        tmpr2 = tmpr2 + pc(n-nn,k)
      enddo
      tmpr1 = tmpr1 / (nt+nn)
      tmpr2 = tmpr2 / (nt+nn)

      summ = 0.
      summ1= 0.
      summ2= 0.
      do n=1, nt+nn
        summ = summ + (inino(n)-tmpr1) * (pc(n-nn,k)-tmpr2)
        summ1= summ1+ (inino(n)-tmpr1)**2
        summ2= summ2+ (pc(n-nn,k)-tmpr2)**2
      enddo
      corrk(mlag+1+nn,k) = summ / sqrt(summ1*summ2)
    enddo
    do nn=1, mlag
      tmpr1= 0.
      tmpr2= 0.
      do n=1, nt-nn
        tmpr1 = tmpr1 + inino(n+nn)
        tmpr2 = tmpr2 + pc(n,k)
      enddo
      tmpr1 = tmpr1 / (nt-nn)
      tmpr2 = tmpr2 / (nt-nn)

      summ = 0.
      summ1= 0.
      summ2= 0.
      do n=1, nt-nn
        summ = summ + (inino(n+nn)-tmpr1) * (pc(n,k)-tmpr2)
        summ1= summ1+ (inino(n+nn)-tmpr1)**2
        summ2= summ2+ (pc(n,k)-tmpr2)**2
      enddo
      corrk(mlag+1+nn,k) = summ / sqrt(summ1*summ2)
    enddo
    lagk(k) = real(maxloc(abs(corrk(mlag+1:,k)),1))-1
    cork(k) = corrk(int(lagk(k))+1+mlag,k)
  enddo
print*, int(lagk)
print*, cork

  do n=1, mlag*2+1
    lag(n) = real(n-mlag-1)
  enddo

  if (jj .eq. 1) then
  call out2d('../coreof.nc',1,(/'Corr'/),corrk(mlag+1:,:5),&
             'lag',mlag+1,lag(mlag+1:),'com',5,(/1.,2.,3.,4.,5./),'correlation')
  end if

if (1 .eq. 1) then
! PC sub. !!!!!!!!!!!!!!!!!!!!!!!!!
  do k=1, 10
    lagg = int(lagk(k))
      tmpr1= 0.
      tmpr2= 0.
      do n=1, nt-lagg
        tmpr1 = tmpr1 + inino(n+lagg)
        tmpr2 = tmpr2 + pc(n,k)
      enddo
      tmpr1 = tmpr1 / (nt-lagg)
      tmpr2 = tmpr2 / (nt-lagg)
      summ = 0.
      summ1= 0.
      summ2= 0.
      do n=1, nt-lagg
        summ = summ + (inino(n+lagg)-tmpr1) * (pc(n,k)-tmpr2)
        summ1= summ1+ (inino(n+lagg)-tmpr1)**2
        summ2= summ2+ (pc(n,k)-tmpr2)**2
      enddo

    if (abs(summ/sqrt(summ1*summ2)) .ge. 0.1) then
    do n=1, nt-lagg
      pc(n,k) = pc(n,k) - inino(n+lagg)*summ/summ1
    enddo
    end if
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

enddo   ! jj

! reconstruct
toare=0.
do n=1, nt
do j=1, nya
do i=1, nxa
do k=1, 10
if (toa(i,j,1) .ne. missvd) then
toare(i,j,n) = toare(i,j,n) + pc(n,k)*reof(i,j,k)
else
toare(i,j,n) = missvd
endif
enddo
enddo
enddo
enddo
call out3d2d('../re.nc',1,(/'re'/),toare,'lon',nxa,lona,'lat',nya,lata,&
             'time',nt,t1,1,(/'PC'/),pc,'time_pc',nt,t1,&
             'com',10,(/1.,2.,3.,4.,5.,6.,7.,8.,9.,10./),'reconstructed')

  ofname1 = '../reofre'
  call reofnc(ofname1,nxa,lona,nya,lata,nt,t1,toare,missn,missvo, &
              10,reof,50,1,normoa)


! toa vs NINO3.4 lag corr.
  do j=1, nya
  do i=1, nxa

  if (toa(i,j,1) .ne. missvd) then
    do nn=-mlag, 0
      tmpr1= 0.
      tmpr2= 0.
      do n=1, nt+nn
        tmpr1 = tmpr1 + inino(n)
        tmpr2 = tmpr2 + toa(i,j,n-nn)
      enddo
      tmpr1 = tmpr1 / (nt+nn)
      tmpr2 = tmpr2 / (nt+nn)

      summ = 0.
      summ1= 0.
      summ2= 0.
      do n=1, nt+nn
        summ = summ + (inino(n)-tmpr1) * (toa(i,j,n-nn)-tmpr2)
        summ1= summ1+ (inino(n)-tmpr1)**2
        summ2= summ2+ (toa(i,j,n-nn)-tmpr2)**2
      enddo
      corr(i,j,mlag+1+nn) = summ / sqrt(summ1*summ2)
!      corr(i,j,mlag+1+nn) = summ / (nt+nn)
    enddo

    do nn=1, mlag
      tmpr1= 0.
      tmpr2= 0.
      do n=1, nt-nn
        tmpr1 = tmpr1 + inino(n+nn)
        tmpr2 = tmpr2 + toa(i,j,n)
      enddo 
      tmpr1 = tmpr1 / (nt-nn)
      tmpr2 = tmpr2 / (nt-nn)

      summ = 0. 
      summ1= 0.
      summ2= 0.
      do n=1, nt-nn 
        summ = summ + (inino(n+nn)-tmpr1) * (toa(i,j,n)-tmpr2)
        summ1= summ1+ (inino(n+nn)-tmpr1)**2
        summ2= summ2+ (toa(i,j,n)-tmpr2)**2
      enddo
      corr(i,j,mlag+1+nn) = summ / sqrt(summ1*summ2)
!      corr(i,j,mlag+1+nn) = summ / (nt-nn)
    enddo
    lag2(i,j) = real(maxloc(abs(corr(i,j,mlag+1:)),1))-1
    cor2(i,j) = corr(i,j,int(lag2(i,j))+1+mlag)
  else
    corr(i,j,:) = missvd
    lag2(i,j) = missvd
    cor2(i,j) = missvd
  end if

  enddo
  enddo

  call out3d('../corlag.nc',1,(/'Corr'/),corr(:,:,mlag+1:),'lon',nxa,lona,&
             'lat',nya,lata,'lag',mlag+1,lag(mlag+1:),'correlation')

  call out2d('../cormax.nc',2,(/'Icorr','Corr'/),(/lag2,cor2/),&
             'lon',nxa,lona,'lat',nya,lata,'correlation')



!  ofname1 = '../eofatl_rem' 
!  call cond_eofnc(ofname1,nxa,lona,nya,lata,nt,t1,toa,missn,missvo,inino,&
!                  30,1,normoa)

!  call opennc('../eofatl_rem.nc',ncid)
!  call geta3d(ncid,'EOF',1,nxa,1,nya,1,10,eof)
!  call closenc(ncid)
!  ofname1 = '../reof_rem'
!  call reofnc(ofname1,nxa,lona,nya,lata,nt,t1,toa,missn,missvo, &
!              10,eof,50,1,normoa)


  stop

end program


subroutine detrend_season(nt,x)

  integer, intent(in)    :: nt
  real,    intent(inout) :: x(nt)

  integer :: i,j,n, nn(12)
  real    :: ave(12)

  n   = 0
  nn  = 0
  ave = 0.
  do i=1, nt
    n = n + 1
    ave(n) = ave(n) + x(i)
    nn(n) = nn(n) + 1
    if (n .eq. 12)  n = 0
  enddo
  ave(:) = ave(:) / nn(:)

  n = 0
  do i=1, nt
    n = n + 1
    x(i) = x(i) - ave(n)
    if (n .eq. 12)  n = 0
  enddo

end subroutine detrend_season

