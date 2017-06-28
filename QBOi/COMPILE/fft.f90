MODULE fft

  implicit none

  interface fft1d_f
    module procedure fft1d_f_r, fft1d_f_c
  end interface
  interface fft1d_b
    module procedure fft1d_b_r, fft1d_b_c
  end interface

  double precision, parameter, private ::  two_pi = 2.d0*                &
                                                  3.14159265358979323846

  contains

FUNCTION waveno_fft(nx,dx)

  integer, intent(in) ::  nx
  real,    intent(in) ::  dx
  real, dimension(nx) ::  waveno_fft
  integer ::  i

  waveno_fft(1) = 0.
  do i=2, nx/2+1
    waveno_fft(i) = float(i-1)
    waveno_fft(nx+2-i) = -waveno_fft(i)
  enddo
  if (nx/2*2 == nx)  waveno_fft(nx/2+1) = abs(waveno_fft(nx/2+1))
  waveno_fft(:) = waveno_fft(:)*real(two_pi/(nx*dx))

END function waveno_fft

SUBROUTINE fft1d_f_r(nx,var,fc,dx,k)

  integer,                intent(in)  ::  nx
  real,    dimension(nx), intent(in)  ::  var
  complex, dimension(nx), intent(out) ::  fc
  real,                   intent(in),  optional ::  dx
  real,    dimension(nx), intent(out), optional ::  k

  if ( present(dx) .and. present(k) ) then
    call fft1d_f_c(nx,cmplx(var),fc,dx,k)
  else
    call fft1d_f_c(nx,cmplx(var),fc)
  end if

END subroutine fft1d_f_r

SUBROUTINE fft1d_f_c(nx,var,fc,dx,k)

  integer,                intent(in)  ::  nx
  complex, dimension(nx), intent(in)  ::  var
  complex, dimension(nx), intent(out) ::  fc
  real,                   intent(in),  optional ::  dx
  real,    dimension(nx), intent(out), optional ::  k

  real, dimension(nx+nx) ::  work
  real, dimension(:), allocatable ::  wsave

  integer ::  i, ier, lensav

  fc(:) = var(:)

  lensav = 2*nx + int(log(real(nx))/log(2.)) + 4
  allocate( wsave(lensav) )

  call cfft1i(nx,wsave,lensav,ier)
  call cfft1f(nx,1,fc,nx,wsave,lensav,work,nx+nx,ier)

  deallocate( wsave )

  if ( present(dx) .and. present(k) )  k(:) = waveno_fft(nx,dx)

END subroutine fft1d_f_c

SUBROUTINE fft1d_b_r(nx,fc,var)

  integer,                intent(in)  ::  nx
  complex, dimension(nx), intent(in)  ::  fc
  real,    dimension(nx), intent(out) ::  var

  complex, dimension(nx) ::  tmp

  call fft1d_b_c(nx,fc,tmp)
  var(:) = real(tmp(:))

END subroutine fft1d_b_r

SUBROUTINE fft1d_b_c(nx,fc,var)

  integer,                intent(in)  ::  nx
  complex, dimension(nx), intent(in)  ::  fc
  complex, dimension(nx), intent(out) ::  var

  real, dimension(nx+nx) ::  work
  real, dimension(:), allocatable ::  wsave

  integer ::  i, ier, lensav

  var(:) = fc(:)

  lensav = 2*nx + int(log(real(nx))/log(2.)) + 4
  allocate( wsave(lensav) )

  call cfft1i(nx,wsave,lensav,ier)
  call cfft1b(nx,1,var,nx,wsave,lensav,work,nx+nx,ier)

  deallocate( wsave )

END subroutine fft1d_b_c

END module fft

