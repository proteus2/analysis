PROGRAM MFS_dir

  use easync

  implicit none

  integer            ::  k_rng(2), l_rng(2), o_rng(2), nphi
  real               ::  lengths(3), c_itv(2), inv_kh_itv, lat_c
  character(len=16)  ::  var_name(6)
  character(len=256) ::  f_namelist, file_i(3), file_o

  namelist /PARAM/ VAR_NAME, K_RNG, L_RNG, O_RNG, LENGTHS, NPHI, C_ITV, &
                   INV_KH_ITV, LAT_C
  namelist /FILEIO/ FILE_I, FILE_O

  integer, parameter ::  nv = 3
  character(len=32), dimension(nv), parameter ::  var_o_name = (/'wu','wv','wuh'/)

  integer ::  i1, j1, n1, nt, nk, nl, no, nc, nkh, n0
  integer ::  i,j,n, iphi, ic, ikh, iphi2, iv
  real    ::  phi_itv

  real, dimension(:,:,:,:), allocatable ::  var_in, mfs
  real, dimension(:,:,:),   allocatable ::  mfs_phi_c, mfs_phi_kh, c_i
  real, dimension(:,:),     allocatable ::  phie_i, phiw_i, kh_i
  real, dimension(:),       allocatable ::  kwn, lwn, ome, phi, c, kh, kwn2
  real, dimension(:),       allocatable ::  phi0, c0, kh0, phi1, c1, kh1

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

! CALCULATE PSD

  mfs_phi_c  = 0.
  mfs_phi_kh = 0.

  n0 = 1
  if (ome(1) == 0.)  n0 = 2

  do j=1, nl
  do i=1, nk
    do iphi=nphi/4+1, nphi*3/4+1  ! from -90 to 90
      if ( phie_i(i,j) .ge. phi0(iphi) .and. phie_i(i,j) .lt. phi1(iphi) ) then
        do n=(no+3)/2, no  ! positive omega
          do ic=1, nc
            if ( c_i(i,j,n) .ge. c0(ic) .and. c_i(i,j,n) .lt. c1(ic) ) then
              mfs_phi_c(iphi,ic,:) = mfs_phi_c(iphi,ic,:) + mfs(i,j,n,:)
            end if
          enddo
        enddo
        do ikh=1, nkh
          if ( kh_i(i,j) .ge. kh0(ikh) .and. kh_i(i,j) .lt. kh1(ikh) ) then
            mfs_phi_kh(iphi,ikh,:) = mfs_phi_kh(iphi,ikh,:) + sum(mfs(i,j,(no+3)/2:,:), dim=1)
          end if
        enddo
        iphi2 = iphi - nphi/2
        if (iphi2 < 1)  iphi2 = iphi + nphi/2
        do n=n0, (no+1)/2  ! negative omega
          do ic=1, nc
            if ( c_i(i,j,n) .ge. c0(ic) .and. c_i(i,j,n) .lt. c1(ic) ) then
              mfs_phi_c(iphi2,ic,:) = mfs_phi_c(iphi2,ic,:) + mfs(i,j,n,:)
            end if
          enddo
        enddo
        do ikh=1, nkh
          if ( kh_i(i,j) .ge. kh0(ikh) .and. kh_i(i,j) .lt. kh1(ikh) ) then
            mfs_phi_kh(iphi2,ikh,:) = mfs_phi_kh(iphi2,ikh,:) + sum(mfs(i,j,n0:(no+1)/2,:), dim=1)
          end if
        enddo
      end if
    enddo
  enddo
  enddo

  if (n0 == 2) then
    do j=1, nl
    do i=1, nk
      do iphi=nphi/4+1, nphi*3/4+1  ! from -90 to 90
        if ( phie_i(i,j) .ge. phi0(iphi) .and. phie_i(i,j) .lt. phi1(iphi) ) then
          iphi2 = iphi - nphi/2
          if (iphi2 < 1)  iphi2 = iphi + nphi/2
          mfs_phi_c(iphi ,1,:) = mfs_phi_c(iphi ,1,:) + 0.5*mfs(i,j,1,:)
          mfs_phi_c(iphi2,1,:) = mfs_phi_c(iphi2,1,:) + 0.5*mfs(i,j,1,:)
          do ikh=1, nkh
            if ( kh_i(i,j) .ge. kh0(ikh) .and. kh_i(i,j) .lt. kh1(ikh) ) then
              mfs_phi_kh(iphi ,ikh,:) = mfs_phi_kh(iphi ,ikh,:) + mfs(i,j,1,:)
              mfs_phi_kh(iphi2,ikh,:) = mfs_phi_kh(iphi2,ikh,:) + mfs(i,j,1,:)
            end if
          enddo
        end if
      enddo
    enddo
    enddo
  end if

  mfs_phi_c (:,:,:) = mfs_phi_c (:,:,:)/(phi_itv*c_itv(1)  )
  mfs_phi_kh(:,:,:) = mfs_phi_kh(:,:,:)/(phi_itv/inv_kh_itv)

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  real ::  deg_hr2m_s, tmp
  real, parameter :: a_earth = 6370.e3

  i1 = k_rng(1)+1  ;  nk = k_rng(2)+1 - i1 + 1
  j1 = l_rng(1)+1  ;  nl = l_rng(2)+1 - j1 + 1
  n1 = o_rng(1)+1  ;  no = o_rng(2)+1 - n1 + 1

  allocate( kwn(nk), kwn2(nk) )
  call get_var(file_i(1),'k_wn',kwn,start=(/i1/), count=(/nk/))
  if ( kwn(1) /= float(k_rng(1)) ) then
    print*, 'k_wn miss match'  ;  STOP
  end if

  allocate( lwn(nl) )
  call get_var(file_i(1),'l_wn',lwn,start=(/j1/), count=(/nl/))
  if ( lwn(1) /= float(l_rng(1)) ) then
    print*, 'l_wn miss match'  ;  STOP
  end if

  allocate( ome(no) )
  call get_var(file_i(1),'ome_fr',ome,start=(/n1/), count=(/no/))
  if ( ome(1) /= float(o_rng(1)) ) then
    print*, 'ome_fr miss match'  ;  STOP
  end if

  kwn(:) = kwn(:)/lengths(1)
  tmp = max(1.,lwn(1)) + lwn(nl)
  lwn((nl+3)/2:) = lwn((nl+3)/2:) - tmp
  lwn(:) = lwn(:)/lengths(2)
  tmp = max(1.,ome(1)) + ome(no)
  ome((no+3)/2:) = ome((no+3)/2:) - tmp
  ome(:) = ome(:)/lengths(3)

  allocate( phie_i(nk,nl), kh_i(nk,nl), c_i(nk,nl,no) )

  kwn2(:) = kwn(:)/cos(lat_c*0.017453293)

  do j=1, nl
  do i=1, nk
    phie_i(i,j) = atan2(lwn(j),kwn2(i))
    kh_i  (i,j) = sqrt(kwn2(i)*kwn2(i) + lwn(j)*lwn(j))
  enddo
  enddo
  phie_i(:,:) = real(phie_i(:,:)*(45.d0/atan(1.d0)))

  deg_hr2m_s = a_earth*real(atan(1.d0)/45.d0)/3600.
  do n=1, no
  do j=1, nl
  do i=1, nk
    c_i(i,j,n) = abs(ome(n))/kh_i(i,j)   ! [deg/hr]
  enddo
  enddo
  enddo
  c_i(:,:,:) = c_i(:,:,:)*deg_hr2m_s   ! [m/s]

  allocate( mfs(nk,nl,no,nv) )
  allocate( var_in(nk,nl,no,4) )

  call get_var(file_i(1),var_name(1),var_in(:,:,:,1),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))
  call get_var(file_i(1),var_name(2),var_in(:,:,:,2),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))

  call get_var(file_i(2),var_name(3),var_in(:,:,:,3),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))
  call get_var(file_i(2),var_name(4),var_in(:,:,:,4),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))

  mfs(:,:,:,1) = var_in(:,:,:,1)*var_in(:,:,:,3) + var_in(:,:,:,2)*var_in(:,:,:,4)

  call get_var(file_i(3),var_name(5),var_in(:,:,:,3),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))
  call get_var(file_i(3),var_name(6),var_in(:,:,:,4),start=(/i1,j1,n1/), &
               count=(/nk,nl,no/), map=(/1,nk,nk*nl/))

  mfs(:,:,:,2) = var_in(:,:,:,1)*var_in(:,:,:,3) + var_in(:,:,:,2)*var_in(:,:,:,4)

  mfs(:,:,:,3) = sqrt( mfs(:,:,:,1)**2 + mfs(:,:,:,2)**2 )

  if (kwn(1) == 0.) then
    mfs(2:nk-1,:,:,:) = mfs(2:nk-1,:,:,:)*2.  ! assume nx is even.
  else
    mfs(1:nk-1,:,:,:) = mfs(1:nk-1,:,:,:)*2.  ! assume nx is even.
  end if

  deallocate( var_in )

  nc = int(c_itv(2)/c_itv(1)) + 1
  nkh = max(kwn(nk),lwn(nl/2+1))*inv_kh_itv
  allocate( phi(nphi), c(nc), kh(nkh) )
  allocate( phi0(nphi), c0(nc), kh0(nkh) )
  allocate( phi1(nphi), c1(nc), kh1(nkh) )

  phi_itv = (360./float(nphi))
  do iphi=1, nphi
    phi(iphi) = float(iphi-1)*phi_itv - 180.  ! [deg]
  enddo
  do ic=1, nc
    c(ic) = float(ic-1)*c_itv(1)  ! [m/s]
  enddo
  do ikh=1, nkh
    kh(ikh) = float(ikh)/inv_kh_itv  ! [cyc/deg]
  enddo

  phi0(:) = phi(:) - 0.5*phi_itv
  phi1(:) = phi(:) + 0.5*phi_itv
  c0  (:) = c  (:) - 0.5*c_itv(1)
  c1  (:) = c  (:) + 0.5*c_itv(1)
  kh0 (:) = kh (:) - 0.5/inv_kh_itv
  kh1 (:) = kh (:) + 0.5/inv_kh_itv

  allocate( mfs_phi_c(nphi,nc,nv), mfs_phi_kh(nphi,nkh,nv) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'dir'   ,phi,axis='dir')
  call put_var('append'   ,file_o,'c'     ,c  ,axis='c'  )
  call put_var('append'   ,file_o,'kh'    ,kh ,axis='kh' )

  do iv=1, nv
    call put_var('append',file_o,trim(var_o_name(iv))//'_phi_c',mfs_phi_c(:,:,iv), &
                 axes=(/'dir','c'/))
    call put_var('append',file_o,trim(var_o_name(iv))//'_phi_kh',mfs_phi_kh(:,:,iv), &
                 axes=(/'dir','kh'/))
  enddo

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( mfs_phi_c, mfs_phi_kh, mfs )
  deallocate( kwn, lwn, ome, phie_i, kh_i, c_i, kwn2 )
  deallocate( phi, c, kh )
  deallocate( phi0, c0, kh0 )
  deallocate( phi1, c1, kh1 )

  END subroutine finalize

END program MFS_dir

