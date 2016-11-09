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

    real, dimension(nz) :: z
    real, dimension(nc) :: cp
    real, dimension(nphi) :: phi
    real, dimension(nphi*10+1) :: newphi
    real, dimension(91) :: newcp
    real, dimension(nphi,nc,nz) :: mfx, mfy
    real, dimension(91,nphi,nz) :: mfx1, mfy1
    real, dimension(91,nphi*10+1,nz) :: mfx2, mfy2
 
    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=2)  ::  speriod

    call getarg(1,speriod)

!== READ 3D MF SPECTRUM =================================================
 
    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_c_phi/mfx_c_phi_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfx'
    call get3d(ncid, varname, nphi, nc, nz, mfx)
    call get1d(ncid, 'direction', nphi, phi)
    call get1d(ncid, 'cp', nc, cp)
    call get1d(ncid, 'height', nz, z)
    call closenc(ncid)

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_c_phi/mfy_c_phi_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfy'
    call get3d(ncid, varname, nphi, nc, nz, mfy)
    call closenc(ncid)

    do i=1, nphi*10+1
      newphi(i)=float(i-1)
    end do
    do i=1, 91
      newcp(i)=float(i-1)
    end do


!== Make data

    do k=1, nz
      do j=1, nphi
        do ii=1, 75
          mfx1(ii,j,k)=0.
          mfy1(ii,j,k)=0.
        end do
        ii=91
        do i=1, 16
          mfx1(ii,j,k)=mfx(j,i,k)
          mfy1(ii,j,k)=mfy(j,i,k)
          ii=ii-1
        end do
      end do
    end do

    do k=1, nz
      do i=1, 91
        j=1
        do jj=1, 351, 10
          mfx2(i,jj,k)=mfx1(i,j,k)
          mfy2(i,jj,k)=mfy1(i,j,k)
          j=j+1
        end do
        mfx2(i,2,k)=mfx2(i,1,k)
        mfx2(i,3,k)=mfx2(i,1,k)
        mfx2(i,4,k)=mfx2(i,1,k)
        mfx2(i,5,k)=mfx2(i,1,k)
        mfx2(i,357,k)=mfx2(i,1,k)
        mfx2(i,358,k)=mfx2(i,1,k)
        mfx2(i,359,k)=mfx2(i,1,k)
        mfx2(i,360,k)=mfx2(i,1,k)
        mfx2(i,361,k)=mfx2(i,1,k)
        mfx2(i,6,k)=(mfx2(i,1,k)+mfx2(i,11,k))/2.
        mfy2(i,2,k)=mfy2(i,1,k)
        mfy2(i,3,k)=mfy2(i,1,k)
        mfy2(i,4,k)=mfy2(i,1,k)
        mfy2(i,5,k)=mfy2(i,1,k)
        mfy2(i,357,k)=mfy2(i,1,k)
        mfy2(i,358,k)=mfy2(i,1,k)
        mfy2(i,359,k)=mfy2(i,1,k)
        mfy2(i,360,k)=mfy2(i,1,k)
        mfy2(i,361,k)=mfy2(i,1,k)
        mfy2(i,6,k)=(mfy2(i,1,k)+mfy2(i,11,k))/2.
 
        do jj=11, 351, 10
          mfx2(i,jj-4,k)=mfx2(i,jj,k)
          mfx2(i,jj-3,k)=mfx2(i,jj,k)
          mfx2(i,jj-2,k)=mfx2(i,jj,k)
          mfx2(i,jj-1,k)=mfx2(i,jj,k)
          mfx2(i,jj+1,k)=mfx2(i,jj,k)
          mfx2(i,jj+2,k)=mfx2(i,jj,k)
          mfx2(i,jj+3,k)=mfx2(i,jj,k)
          mfx2(i,jj+4,k)=mfx2(i,jj,k)
          mfx2(i,jj+5,k)=0.5*(mfx2(i,jj,k)+mfx2(i,jj+10,k))
          mfy2(i,jj-4,k)=mfy2(i,jj,k)
          mfy2(i,jj-3,k)=mfy2(i,jj,k)
          mfy2(i,jj-2,k)=mfy2(i,jj,k)
          mfy2(i,jj-1,k)=mfy2(i,jj,k)
          mfy2(i,jj+1,k)=mfy2(i,jj,k)
          mfy2(i,jj+2,k)=mfy2(i,jj,k)
          mfy2(i,jj+3,k)=mfy2(i,jj,k)
          mfy2(i,jj+4,k)=mfy2(i,jj,k)
          mfy2(i,jj+5,k)=0.5*(mfy2(i,jj,k)+mfy2(i,jj+10,k))
        end do
      end do
      end do


!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/plotdata/mfx_c_phi_plot-data_'//speriod//'_f',cf,'.nc'
    varname='mfx'

    call out3d(ncfilename, varname, mfx2, &
               'cp', 91, newcp, &
               'direction', nphi*10+1, newphi, &
               'height', nz, z, &
               'mfx')

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/plotdata/mfy_c_phi_plot-data_'//speriod//'_f',cf,'.nc'
    varname='mfy'

    call out3d(ncfilename, varname, mfy2, &
               'cp', 91, newcp, &
               'direction', nphi*10+1, newphi, &
               'height', nz, z, &
               'mfy')

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/plotdata/mf_c_phi_plot-data_'//speriod//'_f',cf,'.nc'
    varname='mf'

    call out3d(ncfilename, varname, sqrt(mfx2*mfx2+mfy2*mfy2), &
               'cp', 91, newcp, &
               'direction', nphi*10+1, newphi, &
               'height', nz, z, &
               'mf')

    end program  
