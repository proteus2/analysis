MODULE subr_common

  implicit none

  real, private, parameter :: pi = 3.14159265358979323846
  real, private, parameter :: deg2rad = pi/180.0

CONTAINS


SUBROUTINE basic_u_phi(u,v, u_phi)

  USE param_gwp,  ONLY: nphi, phi_deg

  real, dimension(:,:) , intent(in) ::  u, v

  real, dimension(:,:,:), intent(out) ::  u_phi

  integer ::  iphi

  do iphi=1, nphi
    u_phi(:,:,iphi) = u(:,:)*cos(phi_deg(iphi)*deg2rad) +                   &
                      v(:,:)*sin(phi_deg(iphi)*deg2rad)
  enddo

END subroutine basic_u_phi

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

