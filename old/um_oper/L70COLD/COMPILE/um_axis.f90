MODULE um_axis

  use netio

  implicit none

  character(len=32) ::  axis(3)
  character(len=4)  ::  empty


  contains


SUBROUTINE diminfo(ncid,raw,zh_int,nx,ny,nz,nxu,nyv,nzt,c_axis)

  implicit none

  integer,           intent(in)  ::  ncid
  real,              intent(in)  ::  zh_int
  logical,           intent(in)  ::  raw
  integer,           intent(out) ::  nx, ny, nz, nxu, nyv, nzt
  character(len=32), intent(out) ::  c_axis(3,2)

  integer          ::  tag, temi
  real             ::  temp(1)
  logical          ::  l_err
  character(len=2) ::  c_tail


  if ( raw ) then
    axis(1) = 'longitude'
    axis(2) = 'latitude'
    axis(3) = 'hybrid_ht'
    c_tail  = '_1'
  else
    axis(1) = 'lon'
    axis(2) = 'lat'
    axis(3) = 'zh'
    c_tail  = '_s'
  end if

  empty = 'XXXX'

  c_axis(:,:) = empty

  nxu = 1
  call dilen(ncid,trim(axis(1)), nx,l_err)
  if ( l_err ) then
    nx = 1
  else
    call geta1d(ncid,trim(axis(1)),1,1,temp)
    if ( temp(1) == 0.0 ) then
      c_axis(1,1) = trim(axis(1))
      tag = 2
    else
      c_axis(1,2) = trim(axis(1))
      nxu = nx
      nx  = 1
      tag = 1
    end if
    call dilen(ncid,trim(axis(1))//c_tail, temi,l_err)
    if ( .not. l_err ) then
      c_axis(1,tag) = trim(axis(1))//c_tail
      if (nx == 1) then
        nx = nxu
      else
        nxu = nx
      end if
    end if
  end if

  nyv = 1
  call dilen(ncid,trim(axis(2)), ny,l_err)
  if ( l_err ) then
    ny = 1
  else
    call geta1d(ncid,trim(axis(2)),1,1,temp)
    if ( temp(1) == -90.0 ) then
      c_axis(2,1) = trim(axis(2))
      tag = 2
    else
      c_axis(2,2) = trim(axis(2))
      nyv = ny
      ny  = 1
      tag = 1
    end if
    call dilen(ncid,trim(axis(2))//c_tail, temi,l_err)
    if ( .not. l_err ) then
      c_axis(2,tag) = trim(axis(2))//c_tail
      if (ny == 1) then
        ny = nyv + 1
      else
        nyv = ny - 1
      end if
    end if
  end if

  nzt = 1
  call dilen(ncid,trim(axis(3)), nz,l_err)
  if ( l_err ) then
    nz = 1
  else
    call geta1d(ncid,trim(axis(3)),1,1,temp)
    if ( temp(1) <= zh_int .and. temp(1) > 0. ) then
      c_axis(3,1) = trim(axis(3))
      tag = 2
    else
      c_axis(3,2) = trim(axis(3))
      nzt = nz
      nz  = 1
      tag = 1
    end if
    call dilen(ncid,trim(axis(3))//c_tail, temi,l_err)
    if ( .not. l_err ) then
      c_axis(3,tag) = trim(axis(3))//c_tail
      if (nz == 1) then
        nz = nzt
      else
        nzt = nz
      end if
    end if
  end if


  RETURN

END subroutine diminfo


SUBROUTINE diminfop(ncid,raw,nx,ny,np,nxu,nyv,c_axis)

  implicit none

  integer,           intent(in)  ::  ncid
  logical,           intent(in)  ::  raw
  integer,           intent(out) ::  nx, ny, np, nxu, nyv
  character(len=32), intent(out) ::  c_axis(3,2)

  integer          ::  tag, temi
  real             ::  temp(1)
  logical          ::  l_err
  character(len=2) ::  c_tail


  if ( raw ) then
    axis(1) = 'longitude'
    axis(2) = 'latitude'
    axis(3) = 'p'
    c_tail  = '_1'
  else
    axis(1) = 'lon'
    axis(2) = 'lat'
    axis(3) = 'p'
    c_tail  = '_s'
  end if

  empty = 'XXXX'

  c_axis(:,:) = empty

  nxu = 1
  call dilen(ncid,trim(axis(1)), nx,l_err)
  if ( l_err ) then
    nx = 1
  else
    call geta1d(ncid,trim(axis(1)),1,1,temp)
    if ( temp(1) == 0.0 ) then
      c_axis(1,1) = trim(axis(1))
      tag = 2
    else
      c_axis(1,2) = trim(axis(1))
      nxu = nx
      nx  = 1
      tag = 1
    end if
    call dilen(ncid,trim(axis(1))//c_tail, temi,l_err)
    if ( .not. l_err ) then
      c_axis(1,tag) = trim(axis(1))//c_tail
      if (nx == 1) then
        nx = nxu
      else
        nxu = nx
      end if
    end if
  end if

  nyv = 1
  call dilen(ncid,trim(axis(2)), ny,l_err)
  if ( l_err ) then
    ny = 1
  else
    call geta1d(ncid,trim(axis(2)),1,1,temp)
    if ( temp(1) == -90.0 ) then
      c_axis(2,1) = trim(axis(2))
      tag = 2
    else
      c_axis(2,2) = trim(axis(2))
      nyv = ny
      ny  = 1
      tag = 1
    end if
    call dilen(ncid,trim(axis(2))//c_tail, temi,l_err)
    if ( .not. l_err ) then
      c_axis(2,tag) = trim(axis(2))//c_tail
      if (ny == 1) then
        ny = nyv + 1
      else
        nyv = ny - 1
      end if
    end if
  end if

  call dilen(ncid,trim(axis(3)), np,l_err)
  if ( l_err ) then
    np = 1
  else
    c_axis(3,1) = trim(axis(3))
  end if


  RETURN

END subroutine diminfop


SUBROUTINE axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis,lon,lat,zh,lonu,latv,zht)

  implicit none

  integer,           intent(in)  ::  ncid, nx, ny, nz, nxu, nyv, nzt
  character(len=32), intent(in)  ::  c_axis(3,2)
  real,              intent(out) ::  lon(nx), lat(ny), zh(nz)
  real,              intent(out) ::  lonu(nxu), latv(nyv), zht(nzt)


  if (trim(c_axis(1,1)) /= empty)  call get1d(ncid,trim(c_axis(1,1)),nx, lon)
  if (trim(c_axis(2,1)) /= empty)  call get1d(ncid,trim(c_axis(2,1)),ny, lat)
  if (trim(c_axis(3,1)) /= empty)  call get1d(ncid,trim(c_axis(3,1)),nz, zh )
  if (trim(c_axis(1,2)) /= empty)  call get1d(ncid,trim(c_axis(1,2)),nxu, lonu)
  if (trim(c_axis(2,2)) /= empty)  call get1d(ncid,trim(c_axis(2,2)),nyv, latv)
  if (trim(c_axis(3,2)) /= empty)  call get1d(ncid,trim(c_axis(3,2)),nzt, zht )


  RETURN

END subroutine axisinfo


SUBROUTINE timeinfo(ncid,nt)

  implicit none

  integer, intent(in)  ::  ncid
  integer, intent(out) ::  nt

  logical ::  l_err


  call dilen(ncid,'t',nt, l_err)
  if ( l_err )  nt = 1


  RETURN

END subroutine timeinfo


END module um_axis
