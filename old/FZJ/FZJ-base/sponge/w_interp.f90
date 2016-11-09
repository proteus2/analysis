  program wrf_interp 

  implicit none
  include 'netcdf.inc'

  integer :: file
  integer :: day, hour
  integer :: i, j, k, t
  integer :: ncid

  real :: var1, var2

  integer, parameter :: nfile=30
  integer, parameter :: start_day=7
  integer, parameter :: start_hour=0

!  integer, parameter :: nx=186, ny=186, nz=76, nt=6  ! number of grids
  integer, parameter :: nx=186, ny=186, nz=131, nt=6  ! number of grids
  integer, parameter :: nxs=nx+1, nys=ny+1, nzs=nz+1  ! number of staggered grids

  integer, parameter :: new_nz_phy=101
!  integer, parameter :: new_nz_phy=121
!  integer, parameter :: new_nz_phy=95

!  real, parameter :: dz=450.   ! height-interval in new-coordinate
!  real, parameter :: zt=30000. ! top height
!  real, parameter :: dz=350.   ! height-interval in new-coordinate
  real, parameter :: zt=65000. ! top height
  real, parameter :: dz=500.   ! height-interval in new-coordinate
!  real, parameter :: zt=50000. ! top height

  real, parameter :: zs=15000.  ! bottom height

  integer, parameter :: new_nz=int((zt-zs)/dz)+1     !  number of new levels

  real, dimension(new_nz) :: new_z
  real, dimension(nx) :: axis1
  real, dimension(ny) :: axis2
  real, dimension(new_nz_phy) :: axis3
  real, dimension(nt) :: axis4
  real, dimension(nzs) :: var
  real, dimension(nx,ny,new_nz_phy,nt) :: var_out

  real, allocatable, dimension(:,:,:,:) :: w_org, w_full, w_new
  real, allocatable, dimension(:,:,:,:) :: ph_org, phb_org
  real, allocatable, dimension(:,:,:,:) :: z_half, z_full
  real, allocatable, dimension(:) :: z1, data1 ! temporary 1-D variable for spline

  character (len=40) :: ncfilename   
  character (len=20) :: d1name, d2name, d3name, d4name
  character (len=40) :: title
  character (len=30) :: varname   

!-- Make new levels ------ 
  do k=1, new_nz
    new_z(k)=zs+dz*(k-1)
  end do
!------------------------- 

 print *, '*** top of physical domain **', new_z(new_nz_phy)
  do i=1, nx
    axis1(i)=float(i)
  end do
  do j=1, ny
    axis2(j)=float(j)
  end do
  do k=1, new_nz_phy
    axis3(k)=new_z(k)
  end do
  do t=1, nt
    axis4(t)=float(t)
  end do


  day  = start_day
  hour = start_hour

  do file=1, nfile

!-- Read WRF output files: ex) wrfout_d01_2002-08-30_06:00:00

    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'wrfout_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/'//trim(ncfilename),ncid)

!   perturbation geopotential (m2/s2)
    allocate(ph_org(nx,ny,nzs,nt))
      varname='PH'
      call get4d(ncid, varname, nx, ny, nzs, nt, ph_org)
!   basic-state geopotential (m2/s2)
    allocate(phb_org(nx,ny,nzs,nt))
      varname='PHB'
      call get4d(ncid, varname, nx, ny, nzs, nt, phb_org)
!   geopotential height at half-levels (m)
    allocate(z_half(nx,ny,nz,nt))
    allocate(z_full(nx,ny,nzs,nt))
      do i=1, nx
      do j=1, ny
      do t=1, nt
        do k=1, nz
          z_half(i,j,k,t)=0.5*( ph_org(i,j,k,t)+phb_org(i,j,k,t) &
                               +ph_org(i,j,k+1,t)+phb_org(i,j,k+1,t) )/9.81 
        end do
        do k=1, nzs
          z_full(i,j,k,t)=(ph_org(i,j,k,t)+phb_org(i,j,k,t))/9.81
        end do
      end do
      end do
      end do
    deallocate(ph_org)
    deallocate(phb_org)

print*, minval(z_full(22:126,22:126,nzs,:))

!-- Vertical wind (m/s)

    allocate(w_org(nx,ny,nzs,nt))
      varname='W'
      call get4d(ncid, varname, nx, ny, nzs, nt, w_org)
    allocate(w_full(nx,ny,nzs,nt))
      do i=1, nx
      do j=1, ny
      do k=1, nzs
      do t=1, nt
          w_full(i,j,k,t)=w_org(i,j,k,t)
      end do
      end do
      end do
      end do
    deallocate(w_org)


    allocate(w_new(nx,ny,new_nz,nt))

    allocate(z1(nzs))
    allocate(data1(nzs))
    do i=1, nx
    do j=1, ny
    do t=1, nt
      do k=1, nzs
        data1(k)=w_full(i,j,k,t)
        z1(k)=z_full(i,j,k,t)
      end do
      call spline(z1,data1,nzs,1.e+30,1.e+30,var)
      do k=1, new_nz
        var1=new_z(k)
        call splint(z1,data1,var,nzs,var1,var2)
        w_new(i,j,k,t)=var2
      end do
    end do
    end do
    end do


    deallocate(z1)
    deallocate(data1)
    deallocate(w_full)

    varname='vertical_wind'
    d1name='west_east'
    d2name='south_north'
    d3name='height'
    d4name='time'
    title='Vertical wind of WRF in z-coordinate'

    do k=1, new_nz_phy
      var_out(:,:,k,:)=w_new(:,:,k,:)
    end do

    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call out4d('/data14/kyh/ewiniar/0700-0806/new_data/org/w/'//ncfilename,&
               varname, var_out, &
               d1name, nx, axis1, &
               d2name, ny, axis2, &
               d3name, new_nz_phy, axis3, &
               d4name, nt, axis4, &
               title)

    deallocate(w_new)

    deallocate(z_half)
    deallocate(z_full)


    hour=hour+1
    if(hour.eq.24) then
      hour=0
      day=day+1
    end if
  end do

  end program



!-----------------------------------------------------------------------
!   SUBROUTINE SPLINE
!-----------------------------------------------------------------------

    subroutine spline(x,y,n,yp1,ypn,y2)

    integer :: n
    integer,parameter :: nmax=1000
    integer :: i,k

    real :: yp1,ypn
    real,dimension(n) :: x,y,y2
    real :: p,qn,sig,un
    real,dimension(nmax) :: u

    if(yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
    else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if

    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if(ypn.gt..99e30) then
      qn=0.
      un=0.
    else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
    end do

    end subroutine


!-----------------------------------------------------------------------
!   SUBROUTINE SPLINT
!-----------------------------------------------------------------------

    subroutine splint(xa,ya,y2a,n,x,y)

    integer :: n
    real :: x,y
    real,dimension(n) :: xa,ya,y2a

    integer :: k,khi,klo
    real :: a,b,h
    klo=1
    khi=n

 1  if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x) then
        khi=k
      else
        klo=k
      end if
    goto 1
    end if

     h=xa(khi)-xa(klo)
     if(h.eq.0) pause 'bad xa input in splint'
     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h

     y=a*ya(klo)+b*ya(khi)+((a**3.-a)*y2a(klo)+(b**3.-b)*y2a(khi))*(h**2.)/6.

     end subroutine

!========================================================================

    subroutine opennc(fname,ncid)


    implicit none
    include 'netcdf.inc'
 
    character(len=*), intent(in) :: fname
    integer, intent(out) :: ncid
 
    integer :: st,l
 
    write(*,*) trim(fname)
  
    st = nf_open(trim(fname),nf_nowrite,ncid)
 
    write(*,*) 'NetCDF file id = ',ncid
 
    return
    end subroutine opennc
 
!---------------------------------------------------------------

    subroutine closenc(ncid)     
 
 
    implicit none
    include 'netcdf.inc'
  
    integer, intent(in) :: ncid
    integer :: st
 
    st = nf_close(ncid)
    call errhdr(st)
 
    return
    end subroutine closenc

!---------------------------------------------------------------

    subroutine out1d(wname,name,var,d1name,dim1,axis1,title)
 
    implicit none
 
    include 'netcdf.inc'
 
    character(len=*), intent(in) :: wname,name
    character(len=*), intent(in) :: d1name
    character(len=*), intent(in) :: title
    integer, intent(in) :: dim1
    real, dimension(dim1), intent(in) :: axis1
    real, dimension(dim1), intent(in) :: var
 
    integer :: i,j,k,it,st
    integer :: len1, len2, len3, len4, len5, len6, len7
    integer :: ncid
    integer :: d1id, d2id, d3id, d4id
    integer :: a1id, a2id, a3id, a4id
    integer :: varid
    integer, dimension(1) :: vardimid
 
    len1 = len_trim(wname)
    len2 = len_trim(name)
    len3 = len_trim(d1name)
    len4 = len_trim(title)
 
    st = nf_create(wname(1:len1),nf_share,ncid)
 
    st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
 
    st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
 
    vardimid(1) = d1id
    st = nf_def_var(ncid,name(1:len2),nf_real,1,vardimid,varid)
 
    st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
    st = nf_put_att_text(ncid,nf_global,'title',len4,title)
 
    st = nf_enddef(ncid)
 
    st = nf_put_var_real(ncid,a1id,axis1)
    st = nf_put_var_real(ncid,varid,var)
 
    st = nf_close(ncid)
 
    end subroutine out1d

!---------------------------------------------------------------

    subroutine get1d(ncid,varname,nx,var)
 
    implicit none
 
    include 'netcdf.inc'
 
    integer, intent(in) :: ncid,nx
    character(len=*), intent(in) :: varname
    real, dimension(nx), intent(out) :: var
 
    integer :: varid,st
 
    st = nf_inq_varid(ncid,trim(varname),varid)
    call errhdr(st)
    st = nf_get_var_real(ncid,varid,var)
    call errhdr(st)
 
    write(*,*) trim(varname)
    write(*,*) 'Max : ',maxval(var)
    write(*,*) 'Min : ',minval(var)
 
    return
    end subroutine get1d

!------------------------------------------------------------------------

    subroutine get2d(ncid,varname,nx,ny,var)
 
    implicit none
 
    include 'netcdf.inc'
 
    integer, intent(in) :: ncid,nx,ny
    character(len=*), intent(in) :: varname
    real, dimension(nx,ny), intent(out) :: var
 
    integer :: varid,st
 
    st = nf_inq_varid(ncid,trim(varname),varid)
    call errhdr(st)
    st = nf_get_var_real(ncid,varid,var)
    call errhdr(st)
 
    write(*,*) trim(varname)
    write(*,*) 'Max : ',maxval(var)
    write(*,*) 'Min : ',minval(var)
 
    return
    end subroutine get2d 

!------------------------------------------------------------------------

    subroutine get3d(ncid,varname,nx,ny,nz,var)

    implicit none
 
    include 'netcdf.inc'
 
    integer, intent(in) :: ncid,nx,ny,nz
    character(len=*), intent(in) :: varname
    real, dimension(nx,ny,nz), intent(out) :: var
     
    integer :: varid,st
 
    st = nf_inq_varid(ncid,trim(varname),varid)
    call errhdr(st)
    st = nf_get_var_real(ncid,varid,var)
    call errhdr(st)
 
    write(*,*) trim(varname)
    write(*,*) 'Max : ',maxval(var)
    write(*,*) 'Min : ',minval(var)
 
    return
    end subroutine get3d

!------------------------------------------------------------------------
 
    subroutine get4d(ncid,varname,nx,ny,nz,nt,var)
 
    implicit none
 
    include 'netcdf.inc'

    integer, intent(in) :: ncid,nx,ny,nz,nt
    character(len=*), intent(in) :: varname
    real, dimension(nx,ny,nz,nt), intent(out) :: var
 
    integer :: varid,st
 
    st = nf_inq_varid(ncid,trim(varname),varid)
    call errhdr(st)
    st = nf_get_var_real(ncid,varid,var)
    call errhdr(st)
 
    write(*,*) trim(varname)
    write(*,*) 'Max : ',maxval(var)
    write(*,*) 'Min : ',minval(var)
 
    return
    end subroutine get4d

!------------------------------------------------------------------------

    subroutine out2d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2,  &
                     title)

    implicit none
 
    include 'netcdf.inc'
 
    character(len=*), intent(in) :: wname,name
    character(len=*), intent(in) :: d1name, d2name
    character(len=*), intent(in) :: title
    integer, intent(in) :: dim1, dim2
    real, dimension(dim1), intent(in) :: axis1
    real, dimension(dim2), intent(in) :: axis2
    real, dimension(dim1,dim2), intent(in) :: var
 
    integer :: i,j,k,it,st
    integer :: len1, len2, len3, len4, len5, len6, len7
    integer :: ncid
    integer :: d1id, d2id, d3id, d4id
    integer :: a1id, a2id, a3id, a4id
    integer :: varid
    integer, dimension(2) :: vardimid
 
    len1 = len_trim(wname)
    len2 = len_trim(name)
    len3 = len_trim(d1name)
    len4 = len_trim(d2name)
    len5 = len_trim(title)
 
    st = nf_create(wname(1:len1),nf_share,ncid)
 
    st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
    st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
 
    st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
    st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)

    vardimid(1) = d1id
    vardimid(2) = d2id
    st = nf_def_var(ncid,name(1:len2),nf_real,2,vardimid,varid)
 
    st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
    st = nf_put_att_text(ncid,nf_global,'title',len5,title)
    st = nf_enddef(ncid)
 
    st = nf_put_var_real(ncid,a1id,axis1)
    st = nf_put_var_real(ncid,a2id,axis2)
    st = nf_put_var_real(ncid,varid,var)
 
    st = nf_close(ncid)
    end subroutine out2d

!---------------------------------------------------------------

    subroutine out3d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2, &
                                   d3name,dim3,axis3,title)

!  wname: The file name of the new netCDF dataset
!  name: variable name  
!  var: variable
!  d1name, d2name, d3name: dimension name
!  dim1, dim2, dim3: lengh of dimension
!  axis1, axis2, axis3: The block of data values to be written
 
 
    implicit none
    include 'netcdf.inc'
 
    character(len=*), intent(in) :: wname,name
    character(len=*), intent(in) :: d1name, d2name, d3name 
    character(len=*), intent(in) :: title
    integer, intent(in) :: dim1, dim2, dim3
    real, dimension(dim1), intent(in) :: axis1
    real, dimension(dim2), intent(in) :: axis2
    real, dimension(dim3), intent(in) :: axis3
    real, dimension(dim1,dim2,dim3), intent(in) :: var
 
    integer :: i,j,k,it,st
    integer :: len1, len2, len3, len4, len5, len6, len7
    integer :: ncid
    integer :: d1id, d2id, d3id, d4id
    integer :: a1id, a2id, a3id, a4id
    integer :: varid
    integer, dimension(3) :: vardimid
 
    len1 = len_trim(wname)
    len2 = len_trim(name)
    len3 = len_trim(d1name)
    len4 = len_trim(d2name)
    len5 = len_trim(d3name)
    len6 = len_trim(title)
 
    st = nf_create(wname(1:len1),nf_share,ncid)
 
    st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
    st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
    st = nf_def_dim(ncid,d3name(1:len5),dim3,d3id)
 
    st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
    st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)
    st = nf_def_var(ncid,d3name(1:len5),nf_real,1,d3id,a3id)
 
    vardimid(1) = d1id
    vardimid(2) = d2id
    vardimid(3) = d3id
    st = nf_def_var(ncid,name(1:len2),nf_real,3,vardimid,varid) 
    
    st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
    st = nf_put_att_text(ncid,nf_global,'title',len6,title)
 
    st = nf_enddef(ncid)
 
    st = nf_put_var_real(ncid,a1id,axis1)
    st = nf_put_var_real(ncid,a2id,axis2)
    st = nf_put_var_real(ncid,a3id,axis3)
    st = nf_put_var_real(ncid,varid,var)
 
    st = nf_close(ncid)
    end subroutine out3d


    subroutine out4d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2, &
                                    d3name,dim3,axis3,d4name,dim4,axis4,  &
                                    title)
 
    implicit none
 
    include 'netcdf.inc'
 
    character(len=*), intent(in) :: wname,name
    character(len=*), intent(in) :: d1name, d2name, d3name, d4name
    character(len=*), intent(in) :: title
    integer, intent(in) :: dim1, dim2, dim3, dim4
    real, dimension(dim1), intent(in) :: axis1
    real, dimension(dim2), intent(in) :: axis2
    real, dimension(dim3), intent(in) :: axis3
    real, dimension(dim4), intent(in) :: axis4
    real, dimension(dim1,dim2,dim3,dim4), intent(in) :: var
 
    integer :: i,j,k,it,st
    integer :: len1, len2, len3, len4, len5, len6, len7
    integer :: ncid
    integer :: d1id, d2id, d3id, d4id
    integer :: a1id, a2id, a3id, a4id
    integer :: varid
    integer, dimension(4) :: vardimid
 
    len1 = len_trim(wname)
    len2 = len_trim(name)
    len3 = len_trim(d1name)
    len4 = len_trim(d2name)
    len5 = len_trim(d3name)
    len6 = len_trim(d4name)
    len7 = len_trim(title)
 
    st = nf_create(wname(1:len1),nf_share,ncid)
 
    st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
    st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
    st = nf_def_dim(ncid,d3name(1:len5),dim3,d3id)
    st = nf_def_dim(ncid,d4name(1:len6),dim4,d4id)
 
    st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
    st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)
    st = nf_def_var(ncid,d3name(1:len5),nf_real,1,d3id,a3id)
    st = nf_def_var(ncid,d4name(1:len6),nf_real,1,d4id,a4id)
 
    vardimid(1) = d1id
    vardimid(2) = d2id
    vardimid(3) = d3id
    vardimid(4) = d4id
    st = nf_def_var(ncid,name(1:len2),nf_real,4,vardimid,varid)
 
    st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
    st = nf_put_att_text(ncid,nf_global,'title',len7,title)
 
    st = nf_enddef(ncid)
 
    st = nf_put_var_real(ncid,a1id,axis1)
    st = nf_put_var_real(ncid,a2id,axis2)
    st = nf_put_var_real(ncid,a3id,axis3)
    st = nf_put_var_real(ncid,a4id,axis4)
    st = nf_put_var_real(ncid,varid,var)
 
    st = nf_close(ncid)
    end subroutine out4d
 
 
    subroutine errhdr(ier,loc)
 
    implicit none
 
    include 'netcdf.inc'
 
     integer, intent(in) :: ier
    integer, intent(in), optional :: loc
 
    if (ier /= NF_NOERR) then
      print *,loc,':',NF_STRERROR(ier),' in reading daily mean data file'
      stop
    end if
 
    return
    end subroutine errhdr
