MODULE subr_common

  use netio,  ONLY: vset

  implicit none

  type(vset), dimension(:), allocatable ::  set

  interface defset
    module procedure defset1d, defset2d, defset3d, defset4d
  end interface

CONTAINS


SUBROUTINE defset1d(ivl,varname,vout,axis,ndim,axis1)

  character(len=*)     , intent(in) ::  varname, axis(4)
  real   , dimension(:), intent(in) ::  vout
  integer, dimension(4), intent(in) ::  ndim
  real   , dimension(:), intent(in) ::  axis1

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)

  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,1,1,1) = vout(:)

END subroutine defset1d

SUBROUTINE defset2d(ivl,varname,vout,axis,ndim,axis1,axis2)

  character(len=*)       , intent(in) ::  varname, axis(4)
  real   , dimension(:,:), intent(in) ::  vout
  integer, dimension(4)  , intent(in) ::  ndim
  real   , dimension(:)  , intent(in) ::  axis1, axis2

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)

  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,:,1,1) = vout(:,:)

END subroutine defset2d

subroutine defset3d(ivl,varname,vout,axis,ndim,axis1,axis2,axis3)

  character(len=*)         , intent(in) ::  varname, axis(4)
  real   , dimension(:,:,:), intent(in) ::  vout
  integer, dimension(4)    , intent(in) ::  ndim
  real   , dimension(:)    , intent(in) ::  axis1, axis2, axis3

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)
  set(ivl)%axis3(:) = axis3(:)

  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,:,:,1) = vout(:,:,:)

end subroutine defset3d

subroutine defset4d(ivl,varname,vout,axis,ndim,axis1,axis2,axis3,axis4)

  character(len=*)           , intent(in) ::  varname, axis(4)
  real   , dimension(:,:,:,:), intent(in) ::  vout
  integer, dimension(4)      , intent(in) ::  ndim
  real   , dimension(:)      , intent(in) ::  axis1, axis2, axis3, axis4

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)
  set(ivl)%axis3(:) = axis3(:)
  set(ivl)%axis4(:) = axis4(:)

  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,:,:,:) = vout(:,:,:,:)

end subroutine defset4d


SUBROUTINE setdim0(ivl,varname,axis,ndim)

  character(len=*)     , intent(in) ::  varname, axis(4)
  integer, dimension(4), intent(in) ::  ndim

  integer, intent(inout) ::  ivl

  ivl = ivl + 1

  set(ivl)%vname = trim(varname)

  set(ivl)%axis(:) = axis(:)
  set(ivl)%nd  (:) = ndim(:)
  allocate( set(ivl)%axis1(ndim(1)) )
  allocate( set(ivl)%axis2(ndim(2)) )
  allocate( set(ivl)%axis3(ndim(3)) )
  allocate( set(ivl)%axis4(ndim(4)) )
  set(ivl)%axis1(:) = -999.
  set(ivl)%axis2(:) = -999.
  set(ivl)%axis3(:) = -999.
  set(ivl)%axis4(:) = -999.

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

