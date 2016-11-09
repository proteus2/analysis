PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  1
  integer, parameter ::  nyear  =  7, ystart = 2009
  integer, parameter ::  nivar = 1, nvar = 3
  real,    parameter ::  missv = 1.e32
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'xf150'/
  data                   vartype /'var_xypt'/
  data                   ofdir   /'/data5/kyh/umres/'/ 
  data                   outname /'gphspec'/
  ! input variables
  data                   ivarname(1) /'z'/
  ! output variables
  data                   ovarname(1) /'amp1_z'/
  data                   ovarname(2) /'amp2_z'/
  data                   ovarname(3) /'amp3_z'/


  real, dimension(:,:,:,:),     allocatable ::  var, amp
  real, dimension(:),           allocatable ::  lon, lat, p, lonu, latv, pt
  double complex, dimension(:), allocatable ::  coef

  integer ::  year, month
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, imn, iyr
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz
  integer ::  i,j,k,n,iv, ncid, tag
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=128) ::  fname

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  N_MON:   DO imn=1, nmonth
  N_YEAR:  DO iyr=1, nyear


  year  = ystart + (iyr-1)
  month = mstart + (imn-1)
  if (month > 12)  month = month - 12
  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(vartype)//'.'//cyear//'.'//cmonth//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (imn == 1 .and. iyr == 1) then
    call diminfop(ncid,.FALSE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
  end if
  call timeinfo(ncid, nt)


  ! get var
  do iv=1, nvar

    nxr = NX   ;   nxa =  1
    nyr = NY   ;   nya = NY
    nzr = NZ   ;   nza = NZ  ! ;   iz = 28
    ntr = NT   ;   nta =  1

    allocate( set(iv)%var_out(nxa,nya,nza,nta) )

    if (iv == 1) then

      allocate( var(nxr,nyr,nzr,ntr) )
      call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var)

      allocate( coef(nxr) )
      allocate( amp(nvar,nyr,nzr,ntr) )
      do n=1, ntr
      do k=1, nzr
      do j=1, nyr
        tag = 0
        do i=1, nxr
          if (var(i,j,k,n) == missv)  tag = 1
        enddo
        if (tag == 0) then
          call fft1df(nxr,var(:,j,k,n),coef)
          do i=1, nvar
            amp(i,j,k,n) = 2.*cdabs(coef(i+1))/nxr
          enddo
        else
          amp(:,j,k,n) = missv
        end if
      enddo
      enddo
      enddo
      deallocate ( var )

    end if

    call tempo_avg(1,nyr,nzr,ntr,missv,amp(iv,:,:,:), set(iv)%var_out)

    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta

  enddo

  deallocate ( amp, coef )

  if (imn == 1 .and. iyr == 1) then
    do iv=1, nvar
      set(iv)%axis = (/'     ','lat ','p ','t'/)
      set(iv)%vname = ovarname(iv)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      set(iv)%axis1 = -999.
      set(iv)%axis2 = lat
      set(iv)%axis3 = p
    enddo
  end if
  do iv=1, nvar
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis4 = (year-2000)*100.+month
  enddo

  call closenc(ncid)

  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  print*, fname
  call outnc(trim(fname),nvar,set,'')

  do iv=1, nvar
    deallocate( set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program

