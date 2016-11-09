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

FUNCTION waveno_fft(n,dx)

  integer, intent(in) ::  n
  real,    intent(in) ::  dx
  real, dimension(n)  ::  waveno_fft
  integer ::  i

  waveno_fft(1) = 0.
  do i=2, n/2+1
    waveno_fft(i) = float(i-1)
    waveno_fft(n+2-i) = -waveno_fft(i)
  enddo
  if (n/2*2 == n)  waveno_fft(n/2+1) = abs(waveno_fft(n/2+1))
  waveno_fft(:) = waveno_fft(:)*real(two_pi/(n*dx))

END function waveno_fft

SUBROUTINE fft1d_f_r(var,fc,dx,k)

  real,    dimension(:), intent(in)  ::  var
  complex, dimension(:), intent(out) ::  fc
  real,                  intent(in),  optional ::  dx
  real,    dimension(:), intent(out), optional ::  k

  if ( present(dx) .and. present(k) ) then
    call fft1d_f_c(cmplx(var),fc,dx,k)
  else
    call fft1d_f_c(cmplx(var),fc)
  end if

END subroutine fft1d_f_r

SUBROUTINE fft1d_f_c(var,fc,dx,k)

  complex, dimension(:), intent(in)  ::  var
  complex, dimension(:), intent(out) ::  fc
  real,                  intent(in),  optional ::  dx
  real,    dimension(:), intent(out), optional ::  k

  real, dimension(:), allocatable ::  work, wsave

  integer ::  n, i, ier, lensav

  n = size(var)

  fc(:) = var(:)

  lensav = 2*n + int(log(real(n))/log(2.)) + 4
  allocate( work(n+n), wsave(lensav) )

  call cfft1i(n,wsave,lensav,ier)
  call cfft1f(n,1,fc,n,wsave,lensav,work,n+n,ier)

  deallocate( work, wsave )

  if ( present(dx) .and. present(k) )  k(:) = waveno_fft(n,dx)

END subroutine fft1d_f_c

SUBROUTINE fft1d_b_r(fc,var)

  complex, dimension(:), intent(in)  ::  fc
  real,    dimension(:), intent(out) ::  var

  complex, dimension(size(fc)) ::  tmp

  call fft1d_b_c(fc,tmp)
  var(:) = real(tmp(:))

END subroutine fft1d_b_r

SUBROUTINE fft1d_b_c(fc,var)

  complex, dimension(:), intent(in)  ::  fc
  complex, dimension(:), intent(out) ::  var

  real, dimension(:), allocatable ::  work, wsave

  integer ::  n, i, ier, lensav

  n = size(fc)

  var(:) = fc(:)

  lensav = 2*n + int(log(real(n))/log(2.)) + 4
  allocate( work(n+n), wsave(lensav) )

  call cfft1i(n,wsave,lensav,ier)
  call cfft1b(n,1,var,n,wsave,lensav,work,n+n,ier)

  deallocate( work, wsave )

END subroutine fft1d_b_c

END module fft

