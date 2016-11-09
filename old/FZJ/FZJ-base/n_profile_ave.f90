    program n_profile

!------------------------------------------------------------------------
!
!    Program to calculated horizontal and temporal average N 
!
!------------------------------------------------------------------------

    use netcdfio

    implicit none

    integer :: i, j, k, t, t1, t2
    integer :: ncid
    integer :: hour, day

    integer, parameter :: nxt=186, nyt=186, nz=85
    integer, parameter :: nt1=1!6  ! number of time in each file
    integer, parameter :: nt2=1!29 ! number of files (except list file and spin-up)
    integer, parameter :: start_hour=0
    integer, parameter :: start_day=7
    integer, parameter :: sx=32, ex=116
    integer, parameter :: sy=32, ey=116
    integer, parameter :: nx=ex-sx+1
    integer, parameter :: ny=ey-sy+1  
    integer, parameter :: nt=nt1*nt2! + 1  ! +1 is for last file
 
    character(len=50) :: ncfilename
    character(len=50) :: varname
    character(len=50) :: d1name
    character(len=50) :: title

    real, parameter :: g=9.81

    real, dimension(nz) :: z
    real, allocatable, dimension(:) :: n_ave, n2_ave, t_ave
    real, allocatable, dimension(:,:,:,:) :: th_org,pb_org,pp_org
    real, allocatable, dimension(:,:,:,:) :: th4d,t4d


    allocate(th4d(nx,ny,nz,nt),t4d(nx,ny,nz,nt))

    day=start_day
    hour=start_hour

    t=0
    do t2=1, nt2  ! Loop: file 

    allocate(th_org(nxt,nyt,nz,nt1),pb_org(nxt,nyt,nz,nt1),pp_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_thp_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
!    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
    varname='perturbation_theta'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, 1,nz, 1,nt1, th_org)
    varname='height'
    call get1d(ncid, varname, nz, z)
    call closenc(ncid)

    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    varname='temperature'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, 1,nz, 1,nt1, pb_org)
    call closenc(ncid)

    do t1=1, nt1
      t=t+1
      do k=1, nz, 1 
      do j=sy, ey, 1 
      do i=sx, ex, 1 
        th4d(i-sx+1,j-sy+1,k,t)=th_org(i,j,k,t1)+300.
        t4d (i-sx+1,j-sy+1,k,t)=pb_org(i,j,k,t1)
      end do 
      end do 
      end do 
    end do

    deallocate(th_org, pb_org, pp_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!!-- For last file
!    allocate(th_org(nxt,nyt,nz,1))
!    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_thp_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
!!    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
!    varname='perturbation_theta'
!    call get4d(ncid, varname, nxt, nyt, nz, 1, th_org)
!    call closenc(ncid)
!
!    t=t+1
!    do k=1, nz, 1 
!    do j=sy, ey, 1 
!    do i=sx, ex, 1 
!      th4d(i-sx+1,j-sy+1,k,t)=th_org(i,j,k,1)+300.
!    end do 
!    end do 
!    end do 
!
!    deallocate(th_org)


!-- Calculate domain-average 
   
    allocate(n_ave(nz))
    allocate(n2_ave(nz))

    do k=2, nz-1
      n2_ave(k)=0.
      do t=1, nt
      do i=1, nx
      do j=1, ny
        n2_ave(k)=n2_ave(k)+(g/th4d(i,j,k,t)*(th4d(i,j,k+1,t)-th4d(i,j,k-1,t))/(z(k+1)-z(k-1)))
      end do
      end do
      end do
      n2_ave(k)=n2_ave(k)/float(nx)/float(ny)/float(nt)
    end do

!-- calculate N on upper most level from theta difference between nz-1 and nz 
!-- calculate N on lower most level from theta difference between 2 and 1

    n2_ave(nz)=0.
    n2_ave(1)=0.
      do t=1, nt
      do i=1, nx
      do j=1, ny
        n2_ave(nz)=n2_ave(nz)+(g/th4d(i,j,nz,t)&
                        *(th4d(i,j,nz,t)-th4d(i,j,nz-1,t))/(z(nz)-z(nz-1)))
        n2_ave(1)=n2_ave(1)+(g/th4d(i,j,1,t)&
                        *(th4d(i,j,2,t)-th4d(i,j,1,t))/(z(2)-z(1)))
      end do
      end do
      end do
    n2_ave(nz)=n2_ave(nz)/float(nx)/float(ny)/float(nt)
    n2_ave(1)=n2_ave(1)/float(nx)/float(ny)/float(nt)

    n_ave(:)=sqrt(n2_ave(:))

    allocate(t_ave(nz))
    do k=1, nz
      t_ave(k) = sum(t4d(:,:,k,:))/float(nx*ny*nt)
    enddo


    d1name='height'
    title='n_ave'
    ncfilename='n_profile_ave.nc'
    varname='n_ave'
    call out1d(ncfilename, 2, (/'t_ave',varname/), (/t_ave,n_ave/), &
               d1name, nz, z, &
               title)

    end program  
