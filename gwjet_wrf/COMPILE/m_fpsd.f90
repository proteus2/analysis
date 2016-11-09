MODULE fpsd

  use fft

  implicit none

  double precision, parameter, private ::  two_pi = 2.d0*                &
                                                  3.14159265358979323846


CONTAINS


SUBROUTINE f_psd1d(nx,var,psd,dx,k)

  integer,               intent(in)  ::  nx
  real, dimension(nx),   intent(in)  ::  var
  real, dimension(nx/2), intent(out) ::  psd
  real,                  intent(in),  optional ::  dx
  real, dimension(nx/2), intent(out), optional ::  k

  integer                    ::  nk
  real,    dimension(0:nx-1) ::  k2
  complex, dimension(0:nx-1) ::  fc

  nk = nx/2

  if ( present(dx) .and. present(k) ) then
    call fft1d_f(nx,var,fc,dx,k2)
    k(:) = k2(1:nk)
  else
    call fft1d_f(nx,var,fc)
  end if

  psd(:) = ( real (fc(1:nk))*real (fc(1:nk)) + &
             aimag(fc(1:nk))*aimag(fc(1:nk)) )*2.

  if ( present(dx) )  psd(:) = psd(:)/real(two_pi/(nx*dx))

END subroutine f_psd1d

SUBROUTINE f_psd2d(nx,ny,var,psd,dx,dy,k,l)

  integer,                            intent(in)  ::  nx, ny
  real, dimension(nx,ny),             intent(in)  ::  var
  real, dimension(0:nx/2,-ny/2:ny/2), intent(out) ::  psd
  real,                               intent(in),  optional ::  dx, dy
  real, dimension(0:nx/2),            intent(out), optional ::  k
  real, dimension(-ny/2:ny/2),        intent(out), optional ::  l

  integer                           ::  i,j, nk, nl
  real,    dimension(0:nx-1)        ::  k2
  real,    dimension(0:ny-1)        ::  l2
  real,    dimension(0:nx/2,0:ny-1) ::  psd0
  complex, dimension(0:nx-1,ny)     ::  fc1
  complex, dimension(0:nx/2,0:ny-1) ::  fc2

  nk = nx/2
  nl = ny/2

  if ( present(dx) .and. present(k) .and. &
       present(dy) .and. present(l) ) then
    k2(:) = waveno_fft(nx,dx)
    l2(:) = waveno_fft(ny,dy)
    k(:) = k2(0:nk)
    l(0:nl) = l2(0:nl)
    if ( ny/2*2 /= ny ) then
      l(-nl:-1) = l2(nl+1:)
    else
      l(-nl) = -l(nl)
      l(-nl+1:-1) = l2(nl+1:)
    end if
  end if

  do j=1, ny
    call fft1d_f(nx,var(:,j),fc1(:,j))
  enddo
  do i=0, nk
    call fft1d_f(ny,fc1(i,:),fc2(i,:))
  enddo

  psd0(:,:) = ( real (fc2(:,:))*real (fc2(:,:)) + &
                aimag(fc2(:,:))*aimag(fc2(:,:)) )*2.
  psd0(0,0) = 0.
  psd0(0,:) = 0.5*psd0(0,:)

  psd(:,0:nl) = psd0(:,0:nl)
  if ( ny/2*2 /= ny ) then
    psd(:,-nl:-1) = psd0(:,nl+1:)
  else
    psd(:,-nl+1:-1) = psd0(:,nl+1:)
    psd(:,nl) = 0.5*psd0(:,nl)
    psd(:,-nl) = psd(:,nl)
  end if

  if ( present(dx) .and. present(dy) )  &
     psd(:,:) = psd(:,:)/(real(two_pi/(nx*dx))*real(two_pi/(ny*dy)))

END subroutine f_psd2d

SUBROUTINE f_psd3d(nx,ny,nt,var,psd,dx,dy,dt,k,l,o)

  integer,                     intent(in)  ::  nx, ny, nt
  real, dimension(nx,ny,nt),   intent(in)  ::  var
  real, dimension(0:nx/2,-ny/2:ny/2,-nt/2:nt/2),                         &
                               intent(out) ::  psd
  real,                        intent(in),  optional ::  dx, dy, dt
  real, dimension(0:nx/2),     intent(out), optional ::  k
  real, dimension(-ny/2:ny/2), intent(out), optional ::  l
  real, dimension(-nt/2:nt/2), intent(out), optional ::  o

  integer                       ::  i,j,n, nk, nl, no
  real,    dimension(0:nx-1)                   ::  k2
  real,    dimension(0:ny-1)                   ::  l2
  real,    dimension(0:nt-1)                   ::  o2
  real,    dimension(0:nx/2,0:ny-1,0:nt-1)     ::  psd0
  real,    dimension(0:nx/2,-ny/2:ny/2,0:nt-1) ::  psd1
  complex, dimension(:,:,:), allocatable       ::  fc1, fc2

  nk = nx/2
  nl = ny/2
  no = nt/2

  if ( present(dx) .and. present(k) .and. &
       present(dy) .and. present(l) .and. &
       present(dt) .and. present(o) ) then
    k2(:) = waveno_fft(nx,dx)
    l2(:) = waveno_fft(ny,dy)
    o2(:) = waveno_fft(nt,dt)
    k(:) = k2(0:nk)
    l(0:nl) = l2(0:nl)
    if ( ny/2*2 /= ny ) then
      l(-nl:-1) = l2(nl+1:)
    else
      l(-nl) = -l(nl)
      l(-nl+1:-1) = l2(nl+1:)
    end if
    o(0:no) = o2(0:no)
    if ( nt/2*2 /= nt ) then
      o(-no:-1) = o2(no+1:)
    else
      o(-no) = -o(no)
      o(-no+1:-1) = o2(no+1:)
    end if
  end if

  allocate( fc1(0:nx-1,ny,nt) )
  allocate( fc2(0:nk,0:ny-1,nt) )
  do n=1, nt
  do j=1, ny
    call fft1d_f(nx,var(:,j,n),fc1(:,j,n))
  enddo
  enddo
  do n=1, nt
  do i=0, nk
    call fft1d_f(ny,fc1(i,:,n),fc2(i,:,n))
  enddo
  enddo
  deallocate( fc1 )
  allocate( fc1(0:nk,0:ny-1,0:nt-1) )
  do j=0, ny-1
  do i=0, nk
    call fft1d_f(nt,fc2(i,j,:),fc1(i,j,:))
  enddo
  enddo
  deallocate( fc2 )

  psd0(:,:,:) = ( real (fc1(:,:,:))*real (fc1(:,:,:)) + &
                  aimag(fc1(:,:,:))*aimag(fc1(:,:,:)) )*2.
  psd0(0,:,:) = 0.5*psd0(0,:,:)
  psd0(0,0,:) = 0.  ! The horizontal average is removed at every time.

  deallocate( fc1 )

  psd1(:,0:nl,:) = psd0(:,0:nl,:)
  if ( ny/2*2 /= ny ) then
    psd1(:,-nl:-1,:) = psd0(:,nl+1:,:)
  else
    psd1(:,-nl+1:-1,:) = psd0(:,nl+1:,:)
    psd1(:,nl,:) = 0.5*psd0(:,nl,:)
    psd1(:,-nl,:) = psd1(:,nl,:)
  end if
  psd(:,:,0:no) = psd1(:,:,0:no)
  if ( nt/2*2 /= nt ) then
    psd(:,:,-no:-1) = psd1(:,:,no+1:)
  else
    psd(:,:,-no+1:-1) = psd1(:,:,no+1:)
    psd(:,:,no) = 0.5*psd1(:,:,no)
    psd(:,:,-no) = psd(:,:,no)
  end if

  if ( present(dx) .and. present(dy) .and. present(dt) )  &
     psd(:,:,:) = psd(:,:,:)/ &
        (real(two_pi/(nx*dx))*real(two_pi/(ny*dy))*real(two_pi/(nt*dt)))

END subroutine f_psd3d

END module fpsd
