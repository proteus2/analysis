MODULE ssapack

  use netcdfio

  public :: ssanc

  contains

subroutine ssanc(ofname,nx,x,ny,y,ndata,t,mdata,missn,missv,maxlag,nout, &
                 normopt,norm)

  implicit none

  integer, intent(in) :: nx, ny, ndata, missn, maxlag, nout, normopt
  real,    intent(in) :: x(nx), y(ny), t(ndata), mdata(nx,ny,ndata), missv
  real,    intent(in) :: norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n,ii,jj
  integer :: nch, iiend, tag, nc, nc2, mlch, ier

  real, dimension(ndata,nx*ny-missn) :: datum
  real, dimension(nx*ny-missn) :: avg, sd
  real :: eival(maxlag*(nx*ny-missn))
  real :: eivec(maxlag*(nx*ny-missn),maxlag*(nx*ny-missn))
  real :: eivec4(nx,ny,maxlag,nout)
  real :: stgsym(maxlag*(nx*ny-missn)*(maxlag*(nx*ny-missn)+1)/2)
  real :: rc(ndata,nout,nx*ny-missn), pc(ndata,nout), ts(ndata,nx*ny-missn)
  real :: rc4(nx,ny,ndata,nout), ts3(nx,ny,ndata)
  real :: com(nout), xlag(maxlag), dt
  real :: cov, sum_eival
  real :: tmp1d(maxlag*(nx*ny-missn)*(maxlag*(nx*ny-missn)+1))

  character*3 :: covopt
  data covopt/'cov'/


  nch  = nx * ny - missn
  mlch = maxlag * nch

  nc = 0
  do j=1, ny
  do i=1, nx
    tag = 0
    do n=1, ndata
      if (mdata(i,j,n) .eq. missv)  tag = 1
    enddo
    if (tag .eq. 0) then
      nc = nc + 1
      datum(:,nc) = mdata(i,j,:)
      if (normopt .eq. 1)  datum(:,nc) = mdata(i,j,:) * norm(i,j)
    end if
  enddo
  enddo

  ! consider discrete missing data
!  do j=1, nch
!  do i=2, ndata-1
!    if (datum(i,j) .eq. missv) then
!      if (datum(i-1,j).ne.missv .and. datum(i+1,j).ne.missv) then
!        datum(i,j) = (datum(i-1,j) + datum(i+1,j)) / 2.
!      endif
!    endif
!  enddo
!  enddo

  ! zero mean, variance and cross-variance
  call mean(datum,ndata,nch,missv,covopt,sd,avg)

  ! generate cov-matrix
  nc = 0
  do j=1, nch
  do i=1, maxlag
  do jj=1, j
    if (jj .eq. j)  iiend = i
    if (jj .lt. j)  iiend = maxlag
    do ii=1, iiend
      nc = nc + 1
      call covf(j,i,jj,ii,datum,ndata,cov,maxlag,nch,missv)
      stgsym(nc) = cov / float(mlch)
    enddo
  enddo
  enddo
  enddo

  ! calculate
  call symtrx(stgsym,mlch,eival,eivec,mlch,tmp1d,ier)
  call pcf(datum,ndata,eival,eivec,maxlag,mlch,nch,sd,avg,covopt, &
           missv,nout,nout,pc,rc,ts)

  sum_eival = 0.
  do ii=1, mlch
    sum_eival = sum_eival + eival(ii)
  enddo
  print*, 'sum of the eigen-values :', sum_eival
  eival(:) = eival(:) /sum_eival*100.

  nc = 0  ;  nc2 = 0
  do j=1, ny
  do i=1, nx
    tag = 0
    do n=1, ndata
      if (mdata(i,j,n) .eq. missv)  tag = 1
    enddo
    if (tag .eq. 0) then
      nc = nc + 1
      rc4(i,j,:,:) = rc(:,:nout,nc)
      ts3(i,j,:) = ts(:,nc)
      do ii=1, maxlag
        nc2 = nc2 + 1
        eivec4(i,j,ii,:) = eivec(nc2,:nout)
      enddo
    else
      rc4(i,j,:,:) = missv
      ts3(i,j,:)   = missv
      eivec4(i,j,:,:) = missv
    end if
  enddo
  enddo

  if (normopt .eq. 1) then
    do j=1, ny
    do i=1, nx
      if (ts3(i,j,1) .ne. missv) then
        rc4(i,j,:,:)    = rc4(i,j,:,:)    / norm(i,j)
        ts3(i,j,:)      = ts3(i,j,:)      / norm(i,j)
        eivec4(i,j,:,:) = eivec4(i,j,:,:) / norm(i,j)
      end if
    enddo
    enddo
  end if

  dt = t(2) - t(1)
  do i=1, maxlag
    xlag(i) = real(i*dt)
  enddo
  do j=1, nout
    com(j) = real(j)
  enddo

  ! output
  if (nch .eq. 1) then
    call out2d1d(trim(ofname)//'_compo.nc',2,(/'PC','RC'/),        &
                 (/pc,rc4(1,1,:,:)/),'t',ndata,t,'com',nout,com,   &
                 2,(/'recon','data'/),(/ts3(1,1,:),mdata(1,1,:)/), &
                 't2',ndata,t,                                     &
         'principal and reconstructed compo., recon. and original data')
    call out2d1d(trim(ofname)//'_eigen.nc',1,(/'ei_vec'/),         &
                 eivec4(1,1,:,:),'lag',maxlag,xlag,'com',nout,com, &
                 1,(/'ei_val_100'/),eival,'com2',nout,com,         &
                 'eigen-values and vectors')
  else if (nx .eq. 1) then
    call out3d2d2d(trim(ofname)//'_compo.nc',1,(/'RC'/),rc4(1,:,:,:), &
                   'y',ny,y,'t',ndata,t,'com',nout,com,               &
                   1,(/'PC'/),pc,'t2',ndata,t,'com2',nout,com,        &
                   2,(/'recon','data'/),(/ts3(1,:,:),mdata/),         &
                   'y2',ny,y,'t3',ndata,t,                            &
         'principal and reconstructed compo., recon. and original data')
    call out3d1d(trim(ofname)//'_eigen.nc',1,(/'ei_vec'/),eivec4(1,:,:,:),&
                 'y',ny,y,'lag',maxlag,xlag,'com',nout,com,           &
                 1,(/'ei_val_100'/),eival,'com2',nout,com,            & 
                 'eigen-values and vectors')
  else if (ny .eq. 1) then
    call out3d2d2d(trim(ofname)//'_compo.nc',1,(/'RC'/),rc4(:,1,:,:), &
                   'x',nx,x,'t',ndata,t,'com',nout,com,               &
                   1,(/'PC'/),pc,'t2',ndata,t,'com2',nout,com,        &
                   2,(/'recon','data'/),(/ts3(:,1,:),mdata/),         &
                   'x2',nx,x,'t3',ndata,t,                            &
         'principal and reconstructed compo., recon. and original data')
    call out3d1d(trim(ofname)//'_eigen.nc',1,(/'ei_vec'/),eivec4(:,1,:,:),&
                 'x',nx,x,'lag',maxlag,xlag,'com',nout,com,           & 
                 1,(/'ei_val_100'/),eival,'com2',nout,com,            &
                 'eigen-values and vectors')
  else
    call out4d3d2d(trim(ofname)//'_compo.nc',1,(/'RC'/),rc4,     &
                   'x',nx,x,'y',ny,y,'t',ndata,t,'com',nout,com, &
                   2,(/'recon','data'/),(/ts3,mdata/),           &
                   'x2',nx,x,'y2',ny,y,'t3',ndata,t,             &
                   1,(/'PC'/),pc,'t2',ndata,t,'com2',nout,com,   &
         'principal and reconstructed compo., recon. and original data')
    call out4d1d(trim(ofname)//'_eigen.nc',1,(/'ei_vec'/),eivec4,    &
                 'x',nx,x,'y',ny,y,'lag',maxlag,xlag,'com',nout,com, &
                 1,(/'ei_val_100'/),eival,'com2',nout,com,           &
                 'eigen-values and vectors')
  end if

  return

end subroutine ssanc


END module ssapack
