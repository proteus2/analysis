MODULE kl_filter

  implicit none

  public ::  kl_filter_case

  interface kl_filter_case
    module procedure kl_filter_case_c3d, kl_filter_case_r3d
  end interface

  private

  integer ::  nk, nl
  real    ::  dk, dl

  real,    dimension(:,:),   allocatable ::  det
  real,    dimension(:),     allocatable ::  kwn, lwn
  complex, dimension(:,:,:), allocatable ::  tmp3d

  integer ::  kk0, kk9, ll0, ll9
  integer ::  i,j
  integer ::  nd3d(3), nd4d(4)


  CONTAINS


SUBROUTINE kl_filter_case_c3d(nk_i,nl_i,dk_i,dl_i,icase,coef_kl)

  integer, intent(in) ::  nk_i, nl_i, icase
  real,    intent(in) ::  dk_i, dl_i

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nk = nk_i  ;  nl = nl_i
  dk = dk_i  ;  dl = dl_i

  call get_kl

  select case ( icase )
    case ( 1 )
      call kl_filter_w1(coef_kl)
    case ( 2 )
      call kl_filter_w2(coef_kl)
    case ( 3 )
      call kl_filter_w3(coef_kl)
    case ( 4 )
      call kl_filter_w4(coef_kl)
    case ( 5 )
      call kl_filter_w5(coef_kl)
    case ( 8 )
      call kl_filter_w8(coef_kl)
    case ( 9 )
      call kl_filter_w9(coef_kl)
  end select

  deallocate( tmp3d )
  deallocate( kwn, lwn )

END subroutine kl_filter_case_c3d

SUBROUTINE kl_filter_w8(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = 1
!!  kk9 = int(2.9e-2/dk) + 1
  kk9 = int(0.3e-2/dk) + 1
!!  ll0 = 1
  ll0 = int(0.2e-2/dk) + 1
!!  ll9 = int(2.4e-2/dl) + 1
  ll9 = int(0.6e-2/dl) + 1
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  kk0 = int(0.4e-2/dk) + 2
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = (lwn(j) - 0.2e-2)/(kwn(i) - 0.4e-2) - (1.8-0.2)/(2.9-0.4)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) < 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

  kk0 = 1
  kk9 = int(0.4e-2/dk) + 1
  ll0 = nl+1 - int(0.1e-2/dl)
  ll9 = nl
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)

END subroutine kl_filter_w8

SUBROUTINE kl_filter_w9(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = 1
  kk9 = int(1.8e-2/dk) + 1
  ll0 = 2
!!  ll9 = int(2.8e-2/dl) + 1
  ll9 = int(0.5e-2/dl) + 1
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  kk0 = int(0.3e-2/dk) + 2
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = (lwn(j) - 0.1e-2)/(kwn(i) - 0.3e-2) - (2.8-0.1)/(1.8-0.3)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) < 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

END subroutine kl_filter_w9

SUBROUTINE kl_filter_case_r3d(nk_i,nl_i,dk_i,dl_i,icase,var_kl)

  integer, intent(in) ::  nk_i, nl_i, icase
  real,    intent(in) ::  dk_i, dl_i

  real, dimension(:,:,:), intent(inout) ::  var_kl

  complex, dimension(:,:,:), allocatable ::  varc_kl

  nd3d(:) = shape(var_kl)
  allocate( varc_kl(nd3d(1),nd3d(2),nd3d(3)) )
  varc_kl(:,:,:) = cmplx(var_kl(:,:,:))

  call kl_filter_case_c3d(nk_i,nl_i,dk_i,dl_i,icase,varc_kl)

  var_kl(:,:,:) = real(varc_kl(:,:,:))

  deallocate( varc_kl )

END subroutine kl_filter_case_r3d

SUBROUTINE kl_filter_w1(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = int(0.6e-2/dk) + 2
  kk9 = int(2.9e-2/dk) + 1
  ll0 = nl+1 - int(2.4e-2/dl)
  ll9 = nl - int(0.2e-2/dl)
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  kk0 = int(1.0e-2/dk) + 2
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = (lwn(j) + 0.2e-2)/(kwn(i) - 1.0e-2) - (-1.5+0.2)/(2.9-1.0)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) > 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

END subroutine kl_filter_w1

SUBROUTINE kl_filter_w2(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = int(0.6e-2/dk) + 2
  kk9 = int(2.9e-2/dk) + 1
  ll0 = 2
  ll9 = int(0.8e-2/dl) + 1
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = lwn(j)/(kwn(i) - 0.3e-2) - 1.5/(2.9-0.3)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) > 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo

  kk0 = int(0.9e-2/dk) + 2
  kk9 = int(2.9e-2/dk) + 1
  ll0 = nl+1 - int(1.0e-2/dl)
  ll9 = nl
  coef_kl(kk0:kk9,1,:) = tmp3d(kk0:kk9,1,:)
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = lwn(j)/(kwn(i) - 0.9e-2) + 1.0/(2.9-0.9)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) < 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

END subroutine kl_filter_w2

SUBROUTINE kl_filter_w3(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = 1
  kk9 = int(2.9e-2/dk) + 1
  ll0 = 1
  ll9 = int(2.4e-2/dl) + 1
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  kk0 = int(0.4e-2/dk) + 2
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = (lwn(j) - 0.2e-2)/(kwn(i) - 0.4e-2) - (1.8-0.2)/(2.9-0.4)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) < 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

  kk0 = 1
  kk9 = int(0.4e-2/dk) + 1
  ll0 = nl+1 - int(0.1e-2/dl)
  ll9 = nl
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)

END subroutine kl_filter_w3

SUBROUTINE kl_filter_w4(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = 1
  kk9 = int(1.8e-2/dk) + 1
  ll0 = 2
  ll9 = int(2.8e-2/dl) + 1
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)
  kk0 = int(0.3e-2/dk) + 2
  allocate( det(nk,nl) )
  det(:,:) = 0.
  do j=ll0, ll9
  do i=kk0, kk9
    det(i,j) = (lwn(j) - 0.1e-2)/(kwn(i) - 0.3e-2) - (2.8-0.1)/(1.8-0.3)
  enddo
  enddo
  do j=ll0, ll9
  do i=kk0, kk9
    if ( det(i,j) < 0. )  coef_kl(i,j,:) = 0.
  enddo
  enddo
  deallocate( det )

END subroutine kl_filter_w4

SUBROUTINE kl_filter_w5(coef_kl)

  complex, dimension(:,:,:), intent(inout) ::  coef_kl

  nd3d(:) = shape(coef_kl)
  allocate( tmp3d(nk,nl,nd3d(3)) )

  tmp3d(:,:,:) = coef_kl(:,:,:)
  coef_kl(:,:,:) = 0.

  kk0 = int(0.25e-2/dk) + 2
  kk9 = int(2.0e-2/dk) + 1
  ll0 = nl+1 - int(1.5e-2/dk)
  ll9 = nl - int(0.05e-2/dk)
  coef_kl(kk0:kk9,ll0:ll9,:) = tmp3d(kk0:kk9,ll0:ll9,:)

END subroutine kl_filter_w5

SUBROUTINE get_kl

  allocate( kwn(nk), lwn(nl) )
  do i=1, nk
    kwn(i) = float(i-1)*dk
  enddo
  do j=1, nl/2+1
    lwn(j) = float(j-1)*dl
  enddo
  do j=nl/2+2, nl
    lwn(j) = (-1.)*lwn(nl-j+1)
  enddo

END subroutine get_kl

END module kl_filter

