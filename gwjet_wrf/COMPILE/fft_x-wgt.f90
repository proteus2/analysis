PROGRAM FFT_x_window

  use fft
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2)
  real               ::  yy(2), xx(2), yb(4), xb(4)
  logical            ::  x_periodic
  character(len=16)  ::  var_name
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, II, YY, XX, YB, XB, X_PERIODIC
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, nx, ny, nt, nk, nxx
  integer ::  ib0(2), ib9(2), jb0(2), jb9(2)
  integer ::  i,j,n
  real    ::  lx, ly, x0, y0, xbs(4), ybs(4)

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:),   allocatable ::  var_in
  real,    dimension(:,:),     allocatable ::  wgt
  real,    dimension(:),       allocatable ::  kwn, lat, t, x, y
  complex, dimension(:),       allocatable ::  fc_var

  real, parameter ::  pi = 3.141592

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  allocate( wgt(nxx,ny) )

  do n=1, nt

    xbs(1) = xb(1) + float(n)*xb(4)
    xbs(4) = xb(2) + float(n)*xb(4)
    xbs(2) = xbs(1) + xb(3)
    xbs(3) = xbs(4) - xb(3)
    ybs(1) = yb(1) + float(n)*yb(4)
    ybs(4) = yb(2) + float(n)*yb(4)
    ybs(2) = ybs(1) + yb(3)
    ybs(3) = ybs(4) - yb(3)

    if ( x_periodic ) then
      do while ( xbs(1) < x(1) )
        xbs(:) = xbs(:) + lx
      enddo
      do while ( xbs(4) > x(nxx) )
        xbs(:) = xbs(:) - lx
      enddo
    end if

    wgt(:,:) = 0.

    ib0(:) = 1  ;  ib9(:) = nxx
    jb0(:) = 1  ;  jb9(:) = ny
    do i=1, nxx
      if ( x(i) < xbs(1) ) then
        ib0(1) = i + 1
      end if
      if ( x(i) < xbs(2) ) then
        ib0(2) = i + 1
      end if
      if ( x(i) <= xbs(4) ) then
        ib9(1) = i
      end if
      if ( x(i) <= xbs(3) ) then
        ib9(2) = i
      end if
    enddo
    do j=1, ny
      if ( y(j) < ybs(1) ) then
        jb0(1) = j + 1
      end if
      if ( y(j) < ybs(2) ) then
        jb0(2) = j + 1
      end if
      if ( y(j) <= ybs(4) ) then
        jb9(1) = j
      end if
      if ( y(j) <= ybs(3) ) then
        jb9(2) = j
      end if
    enddo

    wgt(ib0(1):ib9(1),jb0(1):jb9(1)) = 1.

    do j=jb0(1), jb9(1)
      do i=ib0(1), ib0(2)
        wgt(i,j) = abs(cos((x(i)-x(ib0(2)))/(xb(3)*2.)*pi))
      enddo
      do i=ib9(2), ib9(1)
        wgt(i,j) = abs(cos((x(i)-x(ib9(2)))/(xb(3)*2.)*pi))
      enddo
    enddo

    do j=jb0(1), jb0(2)
    do i=ib0(1), ib9(1)
      wgt(i,j) = wgt(i,j)*abs(cos((y(j)-y(jb0(2)))/(yb(3)*2.)*pi))
    enddo
    enddo
    do j=jb9(2), jb9(1)
    do i=ib0(1), ib9(1)
      wgt(i,j) = wgt(i,j)*abs(cos((y(j)-y(jb9(2)))/(yb(3)*2.)*pi))
    enddo
    enddo

    if ( x_periodic ) then
      wgt(1:nx,:) = wgt(1:nx,:) + wgt(nx+1:nxx,:)
!      wgt(nx+1:nxx,:) = wgt(1:nx,:)  ! not used afterward
    end if

    var_in(:,:,n) = var_in(:,:,n)*wgt(1:nx,:)

  enddo

  deallocate( wgt )

  ! calculate FFT
  allocate( fc_var(nx) )
  do n=1, nt
  do j=1, ny
    call fft1d_f(nx,var_in(:,j,n),fc_var)
    var4d(:,j,n,1) = real (fc_var(1:nk))
    var4d(:,j,n,2) = aimag(fc_var(1:nk))
  enddo
  enddo
  deallocate( fc_var )

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

  nk = nx/2 + 1
  allocate( kwn(nk) )
  do i=1, nk
    kwn(i) = float(i-1)
  enddo

  nxx = nx
  if ( x_periodic )  nxx = nx+nx

  allocate( x(nxx), y(ny) )
  lx = xx(2) - xx(1)
  ly = yy(2) - yy(1)
  x0 = xx(1) + 0.5*lx/float(nx)
  y0 = yy(1) - 0.5*ly/float(ny-2)
  do i=1, nxx
    x(i) = x0 + lx*(float(i-1)/float(nx))
  enddo
  do j=1, ny
    y(j) = y0 + ly*(float(j-1)/float(ny-2))
  enddo

  allocate( var_in(nx,ny,nt) )

  call get_var(file_i,var_name,var_in,start=(/i1,j1,n1/), count=(/nx,ny,nt/), &
               stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))

  allocate( var4d(nk,ny,nt,2) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'k_wn',kwn,axis='k_wn')
  call put_var('append'   ,file_o,'lat' ,lat,axis='lat' )
  call put_var('append'   ,file_o,'t'   ,t  ,axis='t'   )

  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_r',var4d(:,:,:,1), &
               axes=(/'k_wn','lat','t'/))
  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_i',var4d(:,:,:,2), &
               axes=(/'k_wn','lat','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var4d, var_in )
  deallocate( kwn )
  deallocate( lat, t, x, y )

  END subroutine finalize

END program FFT_x_window

