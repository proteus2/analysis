program forcing_w
!=============================================================
!  Calculate the nonlinear vertical forcing (Fw)
!-------------------------------------------------------------
!  meanbase : 0;  when defined by w' = w
!             1;  when defined by w' = w - wbar(z,t)
!                      wbar is the zonal mean w
!=============================================================
  
  use netcdfio
 
  implicit none
  
!****************** hyun *************************
  integer, parameter :: nx =1002 , nxn = 802, xst = 100
  integer, parameter :: nz = 472 , nzn = 152, zst = 200  !100
  integer, parameter :: dt =  60   ! [s]
  integer, parameter :: init =     0
  integer, parameter :: fint = 14400
  integer, parameter :: nt = (fint - init) / dt
  integer, parameter :: fwopt = 1, fuopt = 0, meanbase = 0
  integer, parameter :: binopt = 1, ncopt = 0, mixopt = 1
  real,    parameter :: dx = 1000., dz = 300.
!*************************************************
  integer :: i,k,n
  integer :: tmstrln 
  integer :: time 
  real    :: umean, wmean, tmp
  character(len=255) :: fname
  character(len=128) :: timsnd 

  real, dimension(nx,nz) :: prt
  real*8,dimension(nx,nz):: ubar
  real, dimension(nxn) :: xme
  real, dimension(nzn) :: zme, rhobar
  real, dimension(nt)  :: t
  real, dimension(nz)  :: rhodata
  real, dimension(nxn,nzn) :: fw1, fw2, fw, fu1, fu2, fu, wtmix, wmix
  real, dimension(nxn,nzn) :: uprt, wprt, upwp, up2, wp2
  real, dimension(:,:,:), allocatable :: v1, v2, v3, v4, v5
!***************************************************
!***************************************************

  tmstrln = 6
  do i = 1,nxn
    xme(i) = (i+xst-100-1.5) * dx / 1000.
  end do 

  do k = 1, nzn
    zme(k) = (k+zst-1.5) * dz / 1000.
  end do

  do n = 1, nt
    t(n) = (init + n * dt)/60.
  end do

  open(10,file='/usr/users/kyh/ideal/work/idn/rhobar')
  do k=1, nz
    read(10,*) rhodata(k)
  enddo
  rhobar(:) = rhodata(zst+1:zst+nzn)
  open(20,file='/export23/kyh/ideal/idun/bin/ubar.bin', &
          form='unformatted',convert='big_endian')
  read(20) ubar
  close(10)  ;  close(20)


  if (ncopt .eq. 2) then
    allocate(v1(nxn,nzn,nt))
    allocate(v2(nxn,nzn,nt))
    allocate(v3(nxn,nzn,nt))
    if (mixopt .eq. 1) then
      allocate(v4(nxn,nzn,nt))
      allocate(v5(nxn,nzn,nt))
    end if
  end if


  do n=1, nt

    ! read
    time = init + dt * n  
    write(timsnd,'(i6.6)') time 
    open(1,file='/export23/kyh/ideal/idun/bin/uprt.bin'&
     //timsnd(1:tmstrln),form='unformatted',convert='big_endian')
    read(1) prt 
    close(1) 
    uprt(:,:) = prt(xst+1:xst+nxn,zst+1:zst+nzn)
    open(2,file='/export23/kyh/ideal/idun/bin/w.bin'&
     //timsnd(1:tmstrln),form='unformatted',convert='big_endian')
    read(2) prt
    close(2)
    wprt(:,:) = prt(xst+1:xst+nxn,zst+1:zst+nzn)
    if (mixopt .eq. 1) then
      open(3,file='/export23/kyh/ideal/idun/bin/wtmix.bin'&
       //timsnd(1:tmstrln),form='unformatted',convert='big_endian')
      read(3) prt
      close(3)
      wtmix(:,:) = prt(xst+1:xst+nxn,zst+1:zst+nzn)
      open(4,file='/export23/kyh/ideal/idun/bin/wmix.bin'&
       //timsnd(1:tmstrln),form='unformatted',convert='big_endian')
      read(4) prt
      close(4)
      wmix(:,:) = prt(xst+1:xst+nxn,zst+1:zst+nzn)
    end if


    if (meanbase .eq. 1) then
      do k=1, nzn
        umean = 0.
        wmean = 0.
        do i=1, nxn
          umean = umean + uprt(i,k)
          wmean = wmean + wprt(i,k)
        enddo
        umean = umean / nxn
        wmean = wmean / nxn
        do i=1, nxn
          uprt(i,k) = uprt(i,k) - umean
          wprt(i,k) = wprt(i,k) - wmean
        enddo
      enddo
    end if

    ! calculate terms
    do k=1, nzn
    do i=1, nxn
      upwp(i,k) = uprt(i,k) * wprt(i,k)
    enddo
    enddo

    if ( fwopt ) then
    do k=1, nzn
    do i=1, nxn
      wp2(i,k) = wprt(i,k)**2
    enddo
    enddo
    end if
    if ( fuopt ) then
    do k=1, nzn
    do i=1, nxn
      up2(i,k) = uprt(i,k)**2
    enddo
    enddo
    end if

    if ( fwopt ) then
      ! fw1
      do k=1, nzn
        do i=2, nxn-1
          fw1(i,k) = -(upwp(i+1,k)-upwp(i-1,k)) / (2.*dx)
        enddo
        fw1(1,k)   = fw1(2,k)
        fw1(nxn,k) = fw1(nxn-1,k)
      enddo

      ! fw2
      do k=1, nzn
        tmp = 0.
        if (meanbase .eq. 1) then
          do i=1, nxn
            tmp = tmp + wp2(i,k)
          enddo
          tmp = tmp / nxn
        end if
        do i=1, nxn
          wp2(i,k) = (wp2(i,k) - tmp) * rhobar(k)
        enddo
      enddo
      do k=2, nzn-1
      do i=1, nxn
        fw2(i,k) = -(wp2(i,k+1)-wp2(i,k-1)) / (2.*dz) / rhobar(k)
      enddo
      enddo
      fw2(:,1)   = fw2(:,2)
      fw2(:,nzn) = fw2(:,nzn-1)

      ! fw
      fw = fw1 + fw2
    end if

    if ( fuopt ) then
      ! fu1
      do k=1, nzn
        do i=2, nxn-1
          fu1(i,k) = -(up2(i+1,k)-up2(i-1,k)) / (2.*dx)
        enddo
        fu1(1,k)   = fu1(2,k)
        fu1(nxn,k) = fu1(nxn-1,k)
      enddo

      ! fu2
      do k=1, nzn
        tmp = 0.
        if (meanbase .eq. 1) then
          do i=1, nxn
            tmp = tmp + upwp(i,k)
          enddo
          tmp = tmp / nxn
        end if
        do i=1, nxn
          upwp(i,k) = (upwp(i,k) - tmp) * rhobar(k)
        enddo
      enddo
      do k=2, nzn-1
      do i=1, nxn
        fu2(i,k) = -(upwp(i,k+1)-upwp(i,k-1)) / (2.*dz) / rhobar(k)
      enddo
      enddo
      fu2(:,1)   = fu2(:,2)
      fu2(:,nzn) = fu2(:,nzn-1)

      ! fu
      fu = fu1 + fu2
    end if


    ! write
    if (binopt .eq. 1) then
      write(fname,'(a,i6.6)') '/export23/kyh/ideal/idun/fw/fw1.bin',time
      open(101,file=trim(fname), form='unformatted', convert='big_endian')
      write(101) fw1
      close(101)
      write(fname,'(a,i6.6)') '/export23/kyh/ideal/idun/fw/fw2.bin',time
      open(102,file=trim(fname), form='unformatted', convert='big_endian')
      write(102) fw2
      close(102)
      if (mixopt .eq. 1) then
        write(fname,'(a,i6.6)') '/export23/kyh/ideal/idun/fw/fwtmix.bin',time
        open(103,file=trim(fname), form='unformatted', convert='big_endian')
        write(103) wtmix
        close(103)
        write(fname,'(a,i6.6)') '/export23/kyh/ideal/idun/fw/fwmix.bin',time
        open(104,file=trim(fname), form='unformatted', convert='big_endian')
        write(104) wmix
        close(104)
      end if
    end if

    if (ncopt .eq. 1) then

      write(fname,'(a,i5.5,a)') '/export23/kyh/ideal/idun/fwnc/fw',time,'.nc'
      if (mixopt .eq. 0) then
        call out2d(trim(fname),2,(/'fw1 ','fw2'/),(/fw1,fw2/), &
                   'x',nxn,xme,'z',nzn,zme,'w-forcing')
      else
        fw = fw + wtmix
        call out2d(trim(fname),4,(/'fw1 ','fw2','tmix','mix'/), &
                   (/fw1,fw2,wtmix,wmix/),'x',nxn,xme,'z',nzn,zme,'w-forcing')
      end if

    else if (ncopt .eq. 2) then

      v1(:,:,n) = fw1(:,:)
      v2(:,:,n) = fw2(:,:)
      v3(:,:,n) = fw(:,:)
      if (mixopt .eq. 1) then
        v4(:,:,n) = wtmix(:,:)
        v5(:,:,n) = wmix(:,:)
        v3(:,:,n) = fw(:,:) + wtmix(:,:)
      end if

    end if


  enddo   ! n time


    if (ncopt .eq. 2) then
      write(fname,'(a)') '../fw.nc'
      if (mixopt .eq. 0) then
        call out3d(trim(fname),3,(/'fw1 ','fw2','fw'/),(/v1,v2,v3/), &
                   'x',nxn,xme,'z',nzn,zme,'t',nt,t,'w-forcing')
      else
        call out3d(trim(fname),5,(/'fw1 ','fw2','tmix','mix','fw'/), &
                   (/v1,v2,v4,v5,v3/),'x',nxn,xme,'z',nzn,zme,'t',nt,t,'w-forcing')
      end if
    end if


  STOP

end program forcing_w
