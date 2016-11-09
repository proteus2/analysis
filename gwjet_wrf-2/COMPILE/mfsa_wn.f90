PROGRAM MFS_A_dir

  use kl_filter
  use easync

  implicit none

  integer            ::  i_wave, nnn(3), k_rng(2), l_rng(2), nphi0
  real               ::  lengths(2), inv_kh_itv, lat_c
  character(len=16)  ::  var_name(6)
  character(len=256) ::  f_namelist, file_i(3), file_o

  namelist /PARAM/ I_WAVE, NNN, VAR_NAME, K_RNG, L_RNG, LENGTHS, NPHI0, &
                   INV_KH_ITV, LAT_C
  namelist /FILEIO/ FILE_I, FILE_O

  integer, parameter ::  nv = 3
  character(len=32), dimension(nv), parameter ::  var_o_name = (/'wu_a','wv_a','wuh_a'/)
  integer ::  i1, j1, n1, nt, nk, nl, nkh, nphi
  integer ::  i,j,n, iphi, ikh, iphi2, iv
  real    ::  phi_itv, dk, dl

  real, dimension(:,:,:,:), allocatable ::  var_in, mfs, mfs_phi_kh
  real, dimension(:,:),     allocatable ::  phie_i, phiw_i, kh_i
  real, dimension(:),       allocatable ::  t, kwn, lwn, phi, kh, kwn2
  real, dimension(:),       allocatable ::  phi0, kh0, phi1, kh1

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  if ( i_wave /= 0 .and. i_wave /= -999 ) then
    dk = 1./(lengths(1)*111.*cos(lat_c*0.017453293))  ! [cyc/km]
    dl = 1./(lengths(2)*111.)
    do iv=1, nv
      call kl_filter_case(nk,nl,dk,dl,i_wave,mfs(:,:,:,iv))
    enddo
  end if

! CALCULATE PSD

  mfs_phi_kh = 0.

  do j=1, nl
  do i=1, nk
    do iphi=2, nphi
      if ( phie_i(i,j) .ge. phi0(iphi) .and. phie_i(i,j) .lt. phi1(iphi) ) then
        do ikh=1, nkh
          if ( kh_i(i,j) .ge. kh0(ikh) .and. kh_i(i,j) .lt. kh1(ikh) ) then
            mfs_phi_kh(iphi,ikh,:,:) = mfs_phi_kh(iphi,ikh,:,:) + abs(mfs(i,j,:,:))
          end if
        enddo
      end if
    enddo
    if ( phie_i(i,j) .lt. phi1(1) .or. phie_i(i,j) .ge. phi0(nphi) ) then
      do ikh=1, nkh
        if ( kh_i(i,j) .ge. kh0(ikh) .and. kh_i(i,j) .lt. kh1(ikh) ) then
          mfs_phi_kh(1,ikh,:,:) = mfs_phi_kh(1,ikh,:,:) + abs(mfs(i,j,:,:))
        end if
      enddo
    end if
  enddo
  enddo

  mfs_phi_kh(:,:,:,:) = mfs_phi_kh(:,:,:,:)/(phi_itv/inv_kh_itv)

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  real ::  tmp
  real, parameter :: a_earth = 6370.e3

  nphi = nphi0/2

  i1 = k_rng(1)+1  ;  nk = k_rng(2)+1 - i1 + 1
  j1 = l_rng(1)+1  ;  nl = l_rng(2)+1 - j1 + 1
  n1 = nnn(1)  ;  nt = (nnn(2) - n1)/nnn(3) + 1

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

  allocate( t(nt) )
  call get_var(file_i(1),'t',t,start=(/n1/), count=(/nt/), stride=(/nnn(3)/))

  kwn(:) = kwn(:)/lengths(1)
  tmp = max(1.,lwn(1)) + lwn(nl)
  lwn((nl+3)/2:) = lwn((nl+3)/2:) - tmp
  lwn(:) = lwn(:)/lengths(2)

  allocate( phie_i(nk,nl), kh_i(nk,nl) )

  kwn2(:) = kwn(:)/cos(lat_c*0.017453293)

  do j=1, nl
  do i=1, nk
    phie_i(i,j) = atan2(lwn(j),kwn2(i))
    kh_i  (i,j) = sqrt(kwn2(i)*kwn2(i) + lwn(j)*lwn(j))
  enddo
  enddo
  phie_i(:,:) = real(phie_i(:,:)*(45.d0/atan(1.d0)))

  allocate( mfs(nk,nl,nt,nv) )
  allocate( var_in(nk,nl,nt,4) )

  call get_var(file_i(1),var_name(1),var_in(:,:,:,1),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))
  call get_var(file_i(1),var_name(2),var_in(:,:,:,2),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))

  call get_var(file_i(2),var_name(3),var_in(:,:,:,3),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))
  call get_var(file_i(2),var_name(4),var_in(:,:,:,4),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))

  mfs(:,:,:,1) = var_in(:,:,:,1)*var_in(:,:,:,3) + var_in(:,:,:,2)*var_in(:,:,:,4)

  call get_var(file_i(3),var_name(5),var_in(:,:,:,3),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))
  call get_var(file_i(3),var_name(6),var_in(:,:,:,4),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))

  mfs(:,:,:,2) = var_in(:,:,:,1)*var_in(:,:,:,3) + var_in(:,:,:,2)*var_in(:,:,:,4)

  mfs(:,:,:,3) = sqrt( mfs(:,:,:,1)**2 + mfs(:,:,:,2)**2 )

  if (kwn(1) == 0.) then
    mfs(2:nk-1,:,:,:) = mfs(2:nk-1,:,:,:)*2.  ! assume nx is even.
  else
    mfs(1:nk-1,:,:,:) = mfs(1:nk-1,:,:,:)*2.  ! assume nx is even.
  end if

  deallocate( var_in )

  nkh = max(kwn(nk),lwn(nl/2+1))*inv_kh_itv
  allocate( phi(nphi), kh(nkh) )
  allocate( phi0(nphi), kh0(nkh) )
  allocate( phi1(nphi), kh1(nkh) )

  phi_itv = (180./float(nphi))
  do iphi=1, nphi
    phi(iphi) = float(iphi-1)*phi_itv - 90.  ! [deg]
  enddo
  do ikh=1, nkh
    kh(ikh) = float(ikh)/inv_kh_itv  ! [cyc/deg]
  enddo

  phi0(:) = phi(:) - 0.5*phi_itv
  phi1(:) = phi(:) + 0.5*phi_itv
  kh0 (:) = kh (:) - 0.5/inv_kh_itv
  kh1 (:) = kh (:) + 0.5/inv_kh_itv

  allocate( mfs_phi_kh(nphi,nkh,nt,nv) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'dir',phi,axis='dir')
  call put_var('append'   ,file_o,'kh' ,kh ,axis='kh' )
  call put_var('append'   ,file_o,'t'  ,t  ,axis='t'  )

  do iv=1, nv
    call put_var('append',file_o,trim(var_o_name(iv))//'_phi_kh', &
                 mfs_phi_kh(:,:,:,iv), axes=(/'dir','kh','t'/))
  enddo

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( mfs_phi_kh, mfs )
  deallocate( kwn, lwn, phie_i, kh_i, kwn2 )
  deallocate( phi, kh )
  deallocate( phi0, kh0 )
  deallocate( phi1, kh1 )

  END subroutine finalize

END program MFS_A_dir

