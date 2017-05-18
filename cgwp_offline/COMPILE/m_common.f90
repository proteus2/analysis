MODULE subr_common

  USE netio,  ONLY: vset

  implicit none

  type(vset), dimension(:), allocatable ::  set

  integer           ::  ndim1d_put, ndim2d_put(2), ndim3d_put(3),        &
                        ndim4d_put(4)
  character(len=32) ::  axisname1d_put, axisname2d_put(2),               &
                        axisname3d_put(3), axisname4d_put(4)

!  integer                         ::  irec_put = 0
!  integer                         ::  nrec_put = 0
!  character(len=32)               ::  recname_put
!  real, dimension(:), allocatable ::  recval_put

  interface putset
    module procedure putset1d, putset2d, putset3d, putset4d
  end interface

CONTAINS


SUBROUTINE putset1d(ivl,varname,vout,axis1)

  character(len=*)     , intent(in) ::  varname
  real   , dimension(:), intent(in) ::  vout
  real   , dimension(:), intent(in) ::  axis1

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,(/axisname1d_put/),(/ndim1d_put/))

  set(ivl)%axis1(:) = axis1(:)

  set(ivl)%var_out(:,1,1,1) = vout(:)

END subroutine putset1d

SUBROUTINE putset2d(ivl,varname,vout,axis1,axis2)

  character(len=*)       , intent(in) ::  varname
  real   , dimension(:,:), intent(in) ::  vout
  real   , dimension(:)  , intent(in) ::  axis1, axis2

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axisname2d_put,ndim2d_put)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)

  set(ivl)%var_out(:,:,1,1) = vout(:,:)

END subroutine putset2d

SUBROUTINE putset3d(ivl,varname,vout,axis1,axis2,axis3)

  character(len=*)         , intent(in) ::  varname
  real   , dimension(:,:,:), intent(in) ::  vout
  real   , dimension(:)    , intent(in) ::  axis1, axis2, axis3

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axisname3d_put,ndim3d_put)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)
  set(ivl)%axis3(:) = axis3(:)

  set(ivl)%var_out(:,:,:,1) = vout(:,:,:)

END subroutine putset3d

SUBROUTINE putset4d(ivl,varname,vout,axis1,axis2,axis3,axis4)

  character(len=*)           , intent(in) ::  varname
  real   , dimension(:,:,:,:), intent(in) ::  vout
  real   , dimension(:)      , intent(in) ::  axis1, axis2, axis3, axis4

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axisname4d_put,ndim4d_put)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)
  set(ivl)%axis3(:) = axis3(:)
  set(ivl)%axis4(:) = axis4(:)

  set(ivl)%var_out(:,:,:,:) = vout(:,:,:,:)

END subroutine putset4d


SUBROUTINE setdim0(ivl,varname,axisnamexd,ndimxd)

  character(len=*)              , intent(in) ::  varname
  character(len=*), dimension(:), intent(in) ::  axisnamexd
  integer         , dimension(:), intent(in) ::  ndimxd

  integer, intent(inout) ::  ivl

  integer           ::  ndim(4), xd
  character(len=32) ::  axisname(4)

  ivl = ivl + 1

  ndim(:) = 1
  ndim(1:size(ndimxd)) = ndimxd(:)

  axisname(:) = ' '
  axisname(1:size(ndimxd)) = axisnamexd(:)

  set(ivl)%vname = trim(varname)

  set(ivl)%axis(:) = axisname(:)
  set(ivl)%nd  (:) = ndim(:)
  allocate( set(ivl)%axis1(ndim(1)) )
  allocate( set(ivl)%axis2(ndim(2)) )
  allocate( set(ivl)%axis3(ndim(3)) )
  allocate( set(ivl)%axis4(ndim(4)) )
  set(ivl)%axis1(:) = -999.
  set(ivl)%axis2(:) = -999.
  set(ivl)%axis3(:) = -999.
  set(ivl)%axis4(:) = -999.

  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )

END subroutine setdim0


SUBROUTINE column_to_grid3d(icol,jcol,vari,chkdim, varo)

  integer, dimension(:)  , intent(in) ::  icol, jcol
  integer                , intent(in) ::  chkdim
  real   , dimension(:,:), intent(in) ::  vari

  real, dimension(:,:,:), intent(out) ::  varo

  integer ::  ncol
  integer ::  k,l

  ncol = size(icol)

  if ( chkdim /= 0 ) then
    if ( size(vari,1) /= ncol .or. size(varo,3) /= size(vari,2) ) then
      write(6,*) 'ERROR found: column_to_grid3d, 1'  ;  STOP
    end if
    if ( size(varo,1) < maxval(icol) .or.                                &
         size(varo,2) < maxval(jcol) ) then
      write(6,*) 'ERROR found: column_to_grid3d, 2'  ;  STOP
    end if
  end if

  varo(:,:,:) = 0.
  do l=1, ncol
    varo(icol(l),jcol(l),:) = vari(l,:)
  enddo

END subroutine column_to_grid3d

SUBROUTINE column_to_grid2d(icol,jcol,vari,chkdim, varo)

  integer, dimension(:), intent(in) ::  icol, jcol
  integer              , intent(in) ::  chkdim
  real   , dimension(:), intent(in) ::  vari

  real, dimension(:,:), intent(out) ::  varo

  integer ::  ncol
  integer ::  l

  ncol = size(icol)

  if ( chkdim /= 0 ) then
    if ( size(vari) /= ncol ) then
      write(6,*) 'ERROR found: column_to_grid2d, 1'  ;  STOP
    end if
    if ( size(varo,1) < maxval(icol) .or.                                &
         size(varo,2) < maxval(jcol) ) then
      write(6,*) 'ERROR found: column_to_grid2d, 2'  ;  STOP
    end if
  end if

  varo(:,:) = 0.
  do l=1, ncol
    varo(icol(l),jcol(l)) = vari(l)
  enddo

END subroutine column_to_grid2d

SUBROUTINE column_to_grid3d_sp(icol,jcol,vari,chkdim, varo)

  integer, dimension(:)      , intent(in) ::  icol, jcol
  integer                    , intent(in) ::  chkdim
  real   , dimension(:,:,:,:), intent(in) ::  vari

  real, dimension(:,:,:,:,:), intent(out) ::  varo

  integer ::  ncol
  integer ::  i,j,k,l

  ncol = size(icol)

  if ( chkdim /= 0 ) then
    if ( size(vari,2) /= ncol .or. size(varo,4) /= size(vari,3) ) then
      write(6,*) 'ERROR found: column_to_grid3d_sp, 1'  ;  STOP
    end if
    if ( size(varo,2) < maxval(icol) .or.                                &
         size(varo,3) < maxval(jcol) ) then
      write(6,*) 'ERROR found: column_to_grid3d_sp, 2'  ;  STOP
    end if
  end if

  varo(:,:,:,:,:) = 0.
  do l=1, ncol
    i = icol(l)
    j = jcol(l)
    varo(:,i,j,:,:) = vari(:,l,:,:)
  enddo

END subroutine column_to_grid3d_sp

SUBROUTINE column_to_grid2d_sp(icol,jcol,vari,chkdim, varo)

  integer, dimension(:), intent(in) ::  icol, jcol
  integer              , intent(in) ::  chkdim
  real   , dimension(:,:,:), intent(in) ::  vari

  real, dimension(:,:,:,:), intent(out) ::  varo

  integer ::  ncol
  integer ::  i,j,l

  ncol = size(icol)

  if ( chkdim /= 0 ) then
    if ( size(vari,2) /= ncol ) then
      write(6,*) 'ERROR found: column_to_grid2d_sp, 1'  ;  STOP
    end if
    if ( size(varo,2) < maxval(icol) .or.                                &
         size(varo,3) < maxval(jcol) ) then
      write(6,*) 'ERROR found: column_to_grid2d_sp, 2'  ;  STOP
    end if
  end if

  varo(:,:,:,:) = 0.
  do l=1, ncol
    i = icol(l)
    j = jcol(l)
    varo(:,i,j,:) = vari(:,l,:)
  enddo

END subroutine column_to_grid2d_sp

END module subr_common

