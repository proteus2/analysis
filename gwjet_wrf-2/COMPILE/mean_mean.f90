PROGRAM FFT_x

  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2), nphi
  real               ::  c_itv(2)
  character(len=16)  ::  var_name(3)
  character(len=256) ::  f_namelist, file_i(3), file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, II, NPHI, C_ITV
  namelist /FILEIO/ FILE_I, FILE_O

  integer, parameter ::  nv = 3
  character(len=32), dimension(nv), parameter ::  var_o_name = (/'wu','wv','wuh'/)

  integer ::  i1, j1, n1, nx, ny, nt, nc
  integer ::  i,j,n, iphi, ic, iphi2, ic0, iv
  real    ::  phi_itv

  real, dimension(:,:,:,:), allocatable ::  var_in, mfs_phi_c
  real, dimension(:,:,:),   allocatable ::  phie_i, c_i
  real, dimension(:,:,:),   allocatable ::  mfs_phi_c0
  real, dimension(:),       allocatable ::  lat, t
  real, dimension(:),       allocatable ::  phi, c, cosphi, sinphi
  real, dimension(:),       allocatable ::  phi0, c0, phi1, c1, tmp1d

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  allocate( phie_i(nx,ny,nt), c_i(nx,ny,nt) )

  do n=1, nt
  do j=1, ny
  do i=1, nx
    phie_i(i,j,n) = atan2(var_in(i,j,n,3),var_in(i,j,n,2))
  enddo
  enddo
  enddo
  phie_i(:,:,:) = real(phie_i(:,:,:)*(45.d0/atan(1.d0)))

  c_i(:,:,:) = sqrt( var_in(:,:,:,2)*var_in(:,:,:,2) + &
                     var_in(:,:,:,3)*var_in(:,:,:,3) )

  where( abs(var_in(:,:,:,1)) < 10.e-3 )
    c_i(:,:,:) = 1.e20
  end where

  deallocate( var_in )

  mfs_phi_c0 = 0.

  do n=1, nt
  do j=1, ny
  do i=1, nx
    do iphi=2, nphi
      if ( phie_i(i,j,n) .ge. phi0(iphi) .and. phie_i(i,j,n) .lt. phi1(iphi) ) then
        do ic=1, nc
          if ( c_i(i,j,n) .ge. c0(ic) .and. c_i(i,j,n) .lt. c1(ic) ) then
            mfs_phi_c0(iphi,ic,j) = mfs_phi_c0(iphi,ic,j) + 1.
          end if
        enddo
      end if
    enddo
    iphi = 1
    if ( phie_i(i,j,n) .ge. phi1(nphi) .or. phie_i(i,j,n) .lt. phi1(1) ) then
      do ic=1, nc
        if ( c_i(i,j,n) .ge. c0(ic) .and. c_i(i,j,n) .lt. c1(ic) ) then
          mfs_phi_c0(iphi,ic,j) = mfs_phi_c0(iphi,ic,j) + 1.
        end if
      enddo
    end if
  enddo
  enddo
  enddo
  mfs_phi_c0(:,:,:) = mfs_phi_c0(:,:,:)/(float(nx)*float(nt))

  mfs_phi_c0(:,:,:) = mfs_phi_c0(:,:,:)*1.5e-3/(phi_itv*2.5)

  mfs_phi_c = 0.

  do j=1, ny
  do iphi=1, nphi
  do ic=1, nc
    tmp1d(:) = mfs_phi_c0(iphi,:,j)*exp(-((c(ic)-c(:))/30.)**2)
    do ic0=1, nc
      mfs_phi_c(iphi,ic,j,1) = mfs_phi_c(iphi,ic,j,1) + tmp1d(ic0)* &
          (cosphi(iphi)*sign(1.,c(ic)-c(ic0)))
      mfs_phi_c(iphi,ic,j,2) = mfs_phi_c(iphi,ic,j,2) + tmp1d(ic0)* &
          (sinphi(iphi)*sign(1.,c(ic)-c(ic0)))
      mfs_phi_c(iphi,ic,j,3) = mfs_phi_c(iphi,ic,j,3) + tmp1d(ic0)
    enddo
    mfs_phi_c(iphi,ic,j,1) = mfs_phi_c(iphi,ic,j,1) - tmp1d(ic)*cosphi(iphi)
    mfs_phi_c(iphi,ic,j,2) = mfs_phi_c(iphi,ic,j,2) - tmp1d(ic)*sinphi(iphi)
    mfs_phi_c(iphi,ic,j,3) = mfs_phi_c(iphi,ic,j,3) - tmp1d(ic)

    iphi2 = iphi - nphi/2
    if (iphi2 < 1)  iphi2 = iphi + nphi/2
    tmp1d(:) = mfs_phi_c0(iphi2,:,j)*exp(-((c(ic)+c(:))/30.)**2)
    do ic0=1, nc
      mfs_phi_c(iphi,ic,j,1) = mfs_phi_c(iphi,ic,j,1) + tmp1d(ic0)* &
          cosphi(iphi)
      mfs_phi_c(iphi,ic,j,2) = mfs_phi_c(iphi,ic,j,2) + tmp1d(ic0)* &
          sinphi(iphi)
      mfs_phi_c(iphi,ic,j,3) = mfs_phi_c(iphi,ic,j,3) + tmp1d(ic0)
    enddo
  enddo
  enddo
  enddo

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  i1 = ii(1)  ;  j1 = jj(1)  ;  n1 = nnn(1)
  nx = ii(2) - i1 + 1  ;  ny = jj(2) - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  allocate( lat(ny), t(nt) )
  do j=1, ny
    lat(j) = float(j1 + j - 1)
  enddo
  do n=1, nt
    t(n) = float(n1 + (n - 1)*nnn(3))
  enddo

  allocate( var_in(nx,ny,nt,3) )

  call get_var(file_i(1),var_name(1),var_in(:,:,:,1),start=(/i1,j1,n1/), &
               count=(/nx,ny,nt/), stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))
  call get_var(file_i(2),var_name(2),var_in(:,:,:,2),start=(/i1,j1,n1/), &
               count=(/nx,ny,nt/), stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))
  call get_var(file_i(3),var_name(3),var_in(:,:,:,3),start=(/i1,j1,n1/), &
               count=(/nx,ny,nt/), stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))

  nc = int(c_itv(2)/c_itv(1)) + 1
  allocate( phi(nphi), c(nc) )
  allocate( phi0(nphi), c0(nc) )
  allocate( phi1(nphi), c1(nc) )
  allocate( tmp1d(nc) )

  phi_itv = (360./float(nphi))
  do iphi=1, nphi
    phi(iphi) = float(iphi-1)*phi_itv - 180.  ! [deg]
  enddo
  do ic=1, nc
    c(ic) = float(ic-1)*c_itv(1)  ! [m/s]
  enddo

  phi0(:) = phi(:) - 0.5*phi_itv
  phi1(:) = phi(:) + 0.5*phi_itv
  c0  (:) = c  (:) - 0.5*c_itv(1)
  c1  (:) = c  (:) + 0.5*c_itv(1)

  allocate( cosphi(nphi), sinphi(nphi) )
  cosphi(:) = real(cos(phi(:)*(atan(1.d0)/45.d0)))
  sinphi(:) = real(sin(phi(:)*(atan(1.d0)/45.d0)))
  where( abs(cosphi) < 1.e-6 )
    cosphi = 0.
  end where
  where( abs(sinphi) < 1.e-6 )
    sinphi = 0.
  end where

  allocate( mfs_phi_c(nphi,nc,ny,nv), mfs_phi_c0(nphi,nc,ny) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'dir'   ,phi,axis='dir')
  call put_var('append'   ,file_o,'c'     ,c  ,axis='c'  )
  call put_var('append'   ,file_o,'lat'   ,lat,axis='lat')
  call put_var('append'   ,file_o,'t'     ,t  ,axis='t'  )

  do iv=1, nv
    call put_var('append',file_o,trim(var_o_name(iv))//'_phi_c', &
                 mfs_phi_c(:,:,:,iv), axes=(/'dir','c','lat'/))
  enddo

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( mfs_phi_c, mfs_phi_c0 )
  deallocate( lat, t )

  END subroutine finalize

END program FFT_x

