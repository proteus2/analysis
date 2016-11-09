    program mfx_dir 

!   To make 1 deg-intevral data from 10 degree-interval data

    use netcdfio

    implicit none

    integer, parameter ::  cf = 1

    integer :: i, j, k, ii, jj 
    integer :: ncid
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: nz=85
    integer, parameter :: nc=26
    integer, parameter :: nphi=36

    real, dimension(nz) :: z, mfxa, mfxa_p, mfxa_n, mfya, mfya_p, mfya_n
    real, dimension(nc) :: cp
    real, dimension(nphi) :: phi
    real, dimension(nphi,nc,nz) :: mfx, mfy
    real, dimension(91,nphi,nz) :: mfx1, mfy1
    real, dimension(91,nphi*10+1,nz) :: mfx2, mfy2
    real ::  delta
 
    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=2)  ::  speriod

    call getarg(1,speriod)

!== READ 3D MF SPECTRUM =================================================
 
    write(ncfilename,'(a,i0.0,a)') 'res/temp/mf_c_phi/mfx_c_phi_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfx'
    call get3d(ncid, varname, nphi, nc, nz, mfx)
    call get1d(ncid, 'direction', nphi, phi)
    call get1d(ncid, 'cp', nc, cp)
    call get1d(ncid, 'height', nz, z)
    call closenc(ncid)

    write(ncfilename,'(a,i0.0,a)') 'res/temp/mf_c_phi/mfy_c_phi_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfy'
    call get3d(ncid, varname, nphi, nc, nz, mfy)
    call closenc(ncid)



!== Make data

    mfxa_p(:) = 0.
    mfxa_n(:) = 0.
    mfya_p(:) = 0.
    mfya_n(:) = 0.

    do k=1, nz
    do j=1, nc
    do i=1, nphi
      if ( mfx(i,j,k) > 0. ) then
        mfxa_p(k) = mfxa_p(k) + mfx(i,j,k)
      else
        mfxa_n(k) = mfxa_n(k) + mfx(i,j,k)
      end if
      if ( mfy(i,j,k) > 0. ) then
        mfya_p(k) = mfya_p(k) + mfy(i,j,k)
      else
        mfya_n(k) = mfya_n(k) + mfy(i,j,k)
      end if
    enddo
    enddo
    enddo

    delta = (cp(2)-cp(1))*(360./nphi)

    mfxa_p(:) = mfxa_p(:) * delta
    mfxa_n(:) = mfxa_n(:) * delta
    mfya_p(:) = mfya_p(:) * delta
    mfya_n(:) = mfya_n(:) * delta

    mfxa(:) = mfxa_p(:) + mfxa_n(:)
    mfya(:) = mfya_p(:) + mfya_n(:)


!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp/mfa/mfa-c_phi_'//speriod//'_f',cf,'.nc'

    call out1d(ncfilename, 6, &
               (/'mfx  ','mfy','mfx_p','mfx_n','mfy_p','mfy_n'/), &
               (/mfxa,mfya,mfxa_p,mfxa_n,mfya_p,mfya_n/), &
               'height', nz, z, &
               'mf profile')


    end program  
