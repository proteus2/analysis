PROGRAM GPH_monthly_mean_amplitude

  use netio
  use fft

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  7
  integer, parameter ::  nyear  = 10, ystart = 1999
  integer, parameter ::  nwav = 3
  integer, parameter ::  nt_per_day = 4
  real,    parameter ::  missv = 1.e32
  character(len=128) ::  ifdir, ofdir, outname
  character(len=64)  ::  ivarname, ovarname(nwav)
  character(len=32)  ::  c_axis(4)

  ! files
  data                   ifdir   /'/export30/NCEPR2/GPH_6h/'/
  data                   ofdir   /'./out/'/ 
  data                   outname /'gph_amp'/
  ! input variables
  data                   ivarname/'hgt'/
  ! output variables
  data                   ovarname(1) /'amp1_z'/
  data                   ovarname(2) /'amp2_z'/
  data                   ovarname(3) /'amp3_z'/
  ! axis
  data                   c_axis(1) /'lon'/
  data                   c_axis(2) /'lat'/
  data                   c_axis(3) /'level'/
  data                   c_axis(4) /'time'/

  real, dimension(:,:,:,:),     allocatable ::  var, amp
  real, dimension(:),           allocatable ::  lon, lat, lev
  double complex, dimension(:), allocatable ::  coef

  integer ::  year, month, nday(12)
  integer ::  nx, ny, nz, nt, nxa, nya, nza, nta, imn, iyr
  integer ::  i,j,k,n,iw,it, ncid, tag
  logical ::  l_err
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=128) ::  fname

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nwav) ::  set


  nday = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  N_MON:   DO imn=1, nmonth
  N_YEAR:  DO iyr=1, nyear


  year  = ystart + (iyr-1)
  month = mstart + (imn-1)
  if (month > 12)  month = month - 12

  if ( mod(year,4) == 0 ) then
    nday(2) = 29
  else
    nday(2) = 28
  end if

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//'hgt.'//cyear//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (imn == 1 .and. iyr == 1) then
    call dilen(ncid,trim(c_axis(1)), nx,l_err)
    call dilen(ncid,trim(c_axis(2)), ny,l_err)
    call dilen(ncid,trim(c_axis(3)), nz,l_err)
    allocate( lon(nx), lat(ny), lev(nz) )
    call get1d(ncid,trim(c_axis(2)),ny, lat)
    call get1d(ncid,trim(c_axis(3)),nz, lev)
  end if


  ! get var
  nt = nday(month)*nt_per_day
  it = (sum(nday(1:month))-nday(month))*nt_per_day + 1

  allocate( amp(ny,nz,nt,nwav) )
  allocate( coef(nx) )

  allocate( var(nx,ny,nz,nt) )
  call geta4d(ncid,trim(ivarname),1,nx,1,ny,1,nz,it,nt, var)

  ! calculate wave amplitude
  do n=1, nt
  do k=1, nz
  do j=1, ny
    tag = 0
    do i=1, nx
      if (var(i,j,k,n) == missv)  tag = 1
    enddo
    if ( tag ) then
      amp(:,j,k,n) = missv
    else
      call fft1df(nx,var(:,j,k,n),coef)
      do iw=1, nwav
        amp(j,k,n,iw) = 2.*cdabs(coef(iw+1))/nx
      enddo
    end if
  enddo
  enddo
  enddo

  deallocate ( var )
  deallocate ( coef )


  ! monthly average
  do iw=1, nwav
    nxa =  1
    nya = NY
    nza = NZ
    nta =  1

    allocate( set(iw)%var_out(nxa,nya,nza,nta) )
    call tempo_avg(1,ny,nz,nt,missv,amp(:,:,:,iw), set(iw)%var_out)

    set(iw)%nd(1) = nxa
    set(iw)%nd(2) = nya
    set(iw)%nd(3) = nza
    set(iw)%nd(4) = nta
  enddo

  deallocate ( amp )

  ! axis information for dump
  if (imn == 1 .and. iyr == 1) then
    do iw=1, nwav
      set(iw)%axis = (/'     ',trim(c_axis(2)),trim(c_axis(3)),trim(c_axis(4))/)
      set(iw)%vname = ovarname(iw)
      allocate( set(iw)%axis1(set(iw)%nd(1)) )
      allocate( set(iw)%axis2(set(iw)%nd(2)) )
      allocate( set(iw)%axis3(set(iw)%nd(3)) )
      set(iw)%axis1 = -999.
      set(iw)%axis2 = lat
      set(iw)%axis3 = lev
    enddo
  end if
  do iw=1, nwav
    allocate( set(iw)%axis4(set(iw)%nd(4)) )
    set(iw)%axis4 = year*100.+month
  enddo

  call closenc(ncid)

  ! dump
  fname = trim(ofdir)//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  print*, fname
  call outnc(trim(fname),nwav,set,'')

  do iw=1, nwav
    deallocate( set(iw)%axis4 )
    deallocate( set(iw)%var_out )
  enddo


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program


SUBROUTINE tempo_avg(nx,ny,nz,nt,missv,vari,varo)
  
  implicit none
  
  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(nx,ny,nz),    intent(out) ::  varo
  
  integer ::  i,j,k,n 
  integer ::  num(nx,ny,nz)
    

  varo(:,:,:) = 0.
  num (:,:,:) = 0.
  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (vari(i,j,k,n) /= missv) then
      num(i,j,k) = num(i,j,k) + 1
      varo(i,j,k) = varo(i,j,k) + vari(i,j,k,n)
    end if
  enddo
  enddo
  enddo
  enddo

  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (num(i,j,k) /= 0) then
      varo(i,j,k) = varo(i,j,k) / num(i,j,k)
    else
      varo(i,j,k) = missv
    end if
  enddo
  enddo
  enddo


  RETURN

END subroutine tempo_avg

