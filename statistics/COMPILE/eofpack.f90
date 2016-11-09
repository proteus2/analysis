MODULE eofpack

!-----------------------------------------------------------------------
!  Routines for EOF analysis
!-----------------------------------------------------------------------
!  OUTPUT : NetCDF format
!
!  INPUT
!       - ofname  : output file name                        [CHAR]
!       - nx      : no. of the 1st dimension                [INT]
!       - x(nx)   : axis of the 1st dimension               [REAL]
!       - ny      : no. of the 2nd dimension                [INT]
!       - y(ny)   : axis of the 2nd dimension               [REAL]
!       - nt      : no. of time dimension                   [INT]
!       - t(nt)   : time                                    [REAL]
!       - data3d(nx,ny,nt)                                  [REAL]
!                 : data used to analysis
!       - missn   : no. of missing values in x,y domain     [INT]
!       - missv   : missing value                           [REAL]
!       - nout    : no. of component for output             [INT]
!       - normopt : option for normalization                [INT]
!                  == 1, data is normalized by norm(x,y)
!       - norm(nx,ny)                                       [REAL]
!                 : normalization factor
!
!-----------------------------------------------------------------------

  use netcdfio

  public :: eofnc, eeofnc, seofnc, reofnc, ceofnc, cond_eofnc

  contains

subroutine eofnc(ofname,nx,x,ny,y,nt,t,data3d,missn,missv,nout,normopt,norm)
!-----------------------------------------------------------------------
!  EOF analysis
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, nout, normopt
  real,    intent(in) :: x(nx), y(ny), t(nt), data3d(nx,ny,nt)
  real,    intent(in) :: missv, norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, nn, ij, k
  integer :: mij, ier

  real :: datanorm(nx,ny,nt)
  real, dimension((nx*ny-missn)*(nx*ny-missn+1)/2) :: sxy, sx, sy, sx2, sy2
  real, dimension((nx*ny-missn)*(nx*ny-missn+1)) :: wk
  real, dimension(nx*ny-missn) :: zcm, zc, zz, d, evalue100
  real, dimension(nx*ny-missn,nx*ny-missn) :: z
  real, dimension(nx*ny-missn,nout) :: vt
  real, dimension(nx,ny,nout) :: eof

  real :: com(nout), pc(nt,nout), err100(nout), xx, yy, sumd


  mij = nx*ny-missn

  sx = 0.  ;  sy = 0.  ;  sx2 = 0.  ;  sy2 = 0.  ;  sxy = 0.
  zcm = 0.  ;  zz = 0.

  ! normalize
  if (normopt .eq. 1) then
    do n=1, nt
    do j=1, ny
    do i=1, nx
      if (data3d(i,j,n) .ne. missv) then
        datanorm(i,j,n) = data3d(i,j,n) * norm(i,j)
      else
        datanorm(i,j,n) = missv
      end if
    enddo
    enddo
    enddo
  else
    datanorm = data3d
  end if

  do n=1, nt
    ij = 0
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n)
        zcm(ij) = zcm(ij) + zc(ij)/nt
      end if
    enddo
    enddo
    if (mij .ne. ij)  print*, ij, 'error (no. of missing values)'
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx(nn) = sx(nn) + zc(i)
      sy(nn) = sy(nn) + zc(j)
    enddo
    enddo
  enddo

  do n=1, nt
    ij = 0
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n)
        zz(ij) = zz(ij) + (zc(ij)-zcm(ij))**2/nt
      end if
    enddo
    enddo
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx2(nn) = sx2(nn) + (zc(i)-sx(nn)/nt)**2
      sy2(nn) = sy2(nn) + (zc(j)-sy(nn)/nt)**2
      sxy(nn) = sxy(nn) + (zc(i)-sx(nn)/nt)*(zc(j)-sy(nn)/nt)
    enddo
    enddo
  enddo

  do n=1, (mij*(mij+1))/2
    xx = sqrt(sx2(n)/nt)
    yy = sqrt(sy2(n)/nt)
! covar. or correl. ---------------------
    sxy(n) = sxy(n)/nt
!    sxy(n) = sxy(n)/nt / (xx*yy)
! ---------------------------------------
  enddo

  ! calculate eigen-function
  call symtrx(sxy,mij,d,z,mij,wk,ier)

  ! eigen-values (%)
  sumd = 0.
  do k=1, mij
    sumd = sumd + d(k)
  enddo
  do k=1, mij
    evalue100(k) = d(mij-k+1)*100./sumd
  enddo

  ! EOF
  do i=1, mij
  do k=1, nout
    vt(i,k) = z(i,mij-k+1)
  enddo
  enddo

  ij = 0
  do j=1, ny
  do i=1, nx
    if (datanorm(i,j,1) .ne. missv) then
      ij = ij + 1
      eof(i,j,:) = vt(ij,:)
      if (normopt .eq. 1)  eof(i,j,:) = eof(i,j,:) / norm(i,j)
    else
      eof(i,j,:) = missv
    end if
  enddo
  enddo

  ! PC
  pc = 0.
  do n=1, nt
    ij = 0
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n) .ne. missv) then
        ij = ij + 1
        zc(ij)= datanorm(i,j,n) - zcm(ij)
!       zc(ij)= zc(ij)/sqrt(zz(ij))
      end if
    enddo
    enddo

    do k=1, nout
    do i=1, mij
      pc(n,k) = pc(n,k) + zc(i)*vt(i,k)
    enddo
    enddo
  enddo

  do k=1, nout
    com(k) = real(k)
    err100(k) = evalue100(k) * sqrt(2./nt)
  enddo

  ! output
  call out3d2d1d(trim(ofname)//'.nc',1,(/'EOF'/),eof, &
                 'x',nx,x,'y',ny,y,'com_eof',nout,com, &
                 1,(/'PC'/),pc,'t',nt,t,'com_pc',nout,com, &
                 2,(/'e_value_100','err_100'/),(/evalue100(:nout),err100/), &
                 'com_ev',nout,com, 'EOF analysis')


  return

end subroutine eofnc


subroutine eeofnc(ofname,nx,x,ny,y,nt,t,data3d,missn,missv, &
                  nd,dtlag,nout,normopt,norm)
!-----------------------------------------------------------------------
!  Extended EOF analysis
!-----------------------------------------------------------------------
!  INPUT
!       - nd    : no. of dimensions for lag                 [INT]
!       - dtlag : grid interval for time lag                [INT]
!
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, nd, dtlag, nout, normopt
  real,    intent(in) :: x(nx), y(ny), t(nt), data3d(nx,ny,nt)
  real,    intent(in) :: missv, norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, nn, ij, k, nt2
  integer :: mij, ier

  real :: datanorm(nx,ny,nt)
  real, dimension((nx*ny-missn)*nd*((nx*ny-missn)*nd+1)/2) &
                                             :: sxy, sx, sy, sx2, sy2
  real, dimension((nx*ny-missn)*nd*((nx*ny-missn)*nd+1)) :: wk
  real, dimension((nx*ny-missn)*nd) :: zcm, zc, zz, d, evalue100
  real, dimension((nx*ny-missn)*nd,(nx*ny-missn)*nd) :: z
  real, dimension((nx*ny-missn)*nd,nout) :: vt
  real, dimension(nx,ny,nout,nd) :: eof

  real :: com(nout), pc(nt-(nd-1)*dtlag,nout), err100(nout)
  real :: xx, yy, sumd
  character(len=8), dimension(nd) :: varname


  mij = (nx*ny-missn)*nd
  nt2 = nt-(nd-1)*dtlag

  sx = 0.  ;  sy = 0.  ;  sx2 = 0.  ;  sy2 = 0.  ;  sxy = 0.
  zcm = 0.  ;  zz = 0.

  ! normalize
  if (normopt .eq. 1) then
    do n=1, nt
    do j=1, ny
    do i=1, nx
      if (data3d(i,j,n) .ne. missv) then
        datanorm(i,j,n) = data3d(i,j,n) * norm(i,j)
      else
        datanorm(i,j,n) = missv
      end if
    enddo
    enddo
    enddo
  else
    datanorm = data3d
  end if

  do n=1, nt2
    ij = 0
    do k=1, nd
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+(k-1)*dtlag) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n+(k-1)*dtlag)
        zcm(ij) = zcm(ij) + zc(ij)/nt2
      end if
    enddo
    enddo
    enddo
    if (mij .ne. ij)  print*, ij, 'error (no. of missing values)'
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx(nn) = sx(nn) + zc(i)
      sy(nn) = sy(nn) + zc(j)
    enddo
    enddo
  enddo

  do n=1, nt2
    ij = 0
    do k=1, nd
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+(k-1)*dtlag) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n+(k-1)*dtlag)
        zz(ij) = zz(ij) + (zc(ij)-zcm(ij))**2/nt2
      end if
    enddo
    enddo
    enddo
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx2(nn) = sx2(nn) + (zc(i)-sx(nn)/nt2)**2
      sy2(nn) = sy2(nn) + (zc(j)-sy(nn)/nt2)**2
      sxy(nn) = sxy(nn) + (zc(i)-sx(nn)/nt2)*(zc(j)-sy(nn)/nt2)
    enddo
    enddo
  enddo

  do n=1, (mij*(mij+1))/2
    xx = sqrt(sx2(n)/nt2)
    yy = sqrt(sy2(n)/nt2)
! covar. or correl. ---------------------
    sxy(n) = sxy(n)/nt2
!    sxy(n) = sxy(n)/nt2 / (xx*yy)
! ---------------------------------------
  enddo

  ! calculate eigen-function
  call symtrx(sxy,mij,d,z,mij,wk,ier)

  ! eigen-values (%)
  sumd = 0.
  do k=1, mij
    sumd = sumd + d(k)
  enddo
  do k=1, mij
    evalue100(k) = d(mij-k+1)*100./sumd
  enddo

  ! EOF
  do i=1, mij
  do k=1, nout
    vt(i,k) = z(i,mij-k+1)
  enddo
  enddo

  ij = 0
  do k=1, nd
  do j=1, ny
  do i=1, nx
    if (datanorm(i,j,1) .ne. missv) then
      ij = ij + 1
      eof(i,j,:,k) = vt(ij,:)
      if (normopt .eq. 1)  eof(i,j,:,k) = eof(i,j,:,k) / norm(i,j)
    else
      eof(i,j,:,k) = missv
    end if
  enddo
  enddo
  enddo

  ! PC
  pc = 0.
  do n=1, nt2
    ij = 0
    do k=1, nd
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+(k-1)*dtlag) .ne. missv) then
        ij = ij + 1
        zc(ij)= datanorm(i,j,n+(k-1)*dtlag) - zcm(ij)
!       zc(ij)= zc(ij)/sqrt(zz(ij))
      end if
    enddo
    enddo
    enddo

    do k=1, nout
    do i=1, mij
      pc(n,k) = pc(n,k) + zc(i)*vt(i,k)
    enddo
    enddo
  enddo

  do k=1, nout
    com(k) = real(k)
    err100(k) = evalue100(k) * sqrt(2./nt2)
  enddo
  do k=1, nd
    if (k .ge. 10) then
      write(varname(k),'(a,i2)') 'EEOF_t', k
    else
      write(varname(k),'(a,i1)') 'EEOF_t', k
    end if
  enddo

  ! output
  call out3d2d1d(trim(ofname)//'.nc',nd,varname,eof, &
                 'x',nx,x,'y',ny,y,'com_eof',nout,com, &
                 1,(/'PC'/),pc,'t',nt2,t(:nt2),'com_pc',nout,com, &
                 2,(/'e_value_100','err_100'/),(/evalue100(:nout),err100/), &
                 'com_ev',nout,com, 'Extended EOF analysis')


  return

end subroutine eeofnc


subroutine seofnc(ofname,nx,x,ny,y,nt,tinf,data3d,missn,missv, &
                  nout,normopt,norm)
!-----------------------------------------------------------------------
!  Season-dependent EOF analysis
!-----------------------------------------------------------------------
!  INPUT
!       - tinf(3)                                           [INT]
!              // tinf(1) : time interval of data in month
!                 tinf(2) : 1st month of data
!                 tinf(3) : 1st year of data
!   
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, nout, normopt
  integer, intent(in) :: tinf(3)
  real,    intent(in) :: x(nx), y(ny), data3d(nx,ny,nt)
  real,    intent(in) :: missv, norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, nn, ij, k, nt2
  integer :: mij, ier, tmp

  real :: datanorm(nx,ny,nt), t(nt*tinf(1)/12)
  real, dimension((nx*ny-missn)*4*((nx*ny-missn)*4+1)/2) &
                                             :: sxy, sx, sy, sx2, sy2
  real, dimension((nx*ny-missn)*4*((nx*ny-missn)*4+1)) :: wk
  real, dimension((nx*ny-missn)*4) :: zcm, zc, zz, d, evalue100
  real, dimension((nx*ny-missn)*4,(nx*ny-missn)*4) :: z
  real, dimension((nx*ny-missn)*4,nout) :: vt
  real, dimension(nx,ny,nout,4) :: eof

  real :: com(nout), pc(nt*tinf(1)/12,nout), err100(nout), xx, yy, sumd
  character(len=8), dimension(4) :: varname


  if (tinf(1).ne.1 .and. tinf(1).ne.3) then
    print'(//,a,/)', 'ERROR: tinf'
    STOP
  end if

  mij = (nx*ny-missn)*4
  nt2 = nt*tinf(1)/12

  sx = 0.  ;  sy = 0.  ;  sx2 = 0.  ;  sy2 = 0.  ;  sxy = 0.
  zcm = 0.  ;  zz = 0.

  ! normalize
  if (normopt .eq. 1) then
    do n=1, nt
    do j=1, ny
    do i=1, nx
      if (data3d(i,j,n) .ne. missv) then
        datanorm(i,j,n) = data3d(i,j,n) * norm(i,j)
      else
        datanorm(i,j,n) = missv
      end if
    enddo
    enddo
    enddo
  else
    datanorm = data3d
  end if

  do n=1, nt2*12/tinf(1), 12/tinf(1)
    ij = 0
    do k=0, 9/tinf(1), 3/tinf(1)
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+k) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n+k)
        zcm(ij) = zcm(ij) + zc(ij)/nt2
      end if
    enddo
    enddo
    enddo
    if (mij .ne. ij)  print*, ij, 'error (no. of missing values)'
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx(nn) = sx(nn) + zc(i)
      sy(nn) = sy(nn) + zc(j)
    enddo
    enddo
  enddo

  do n=1, nt2*12/tinf(1), 12/tinf(1)
    ij = 0
    do k=0, 9/tinf(1), 3/tinf(1)
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+k) .ne. missv) then
        ij = ij + 1
        zc(ij) = datanorm(i,j,n+k)
        zz(ij) = zz(ij) + (zc(ij)-zcm(ij))**2/nt2
      end if
    enddo
    enddo
    enddo
    nn = 0
    do j=1, mij
    do i=1, j
      nn = nn + 1
      sx2(nn) = sx2(nn) + (zc(i)-sx(nn)/nt2)**2
      sy2(nn) = sy2(nn) + (zc(j)-sy(nn)/nt2)**2
      sxy(nn) = sxy(nn) + (zc(i)-sx(nn)/nt2)*(zc(j)-sy(nn)/nt2)
    enddo
    enddo
  enddo

  do n=1, (mij*(mij+1))/2
    xx = sqrt(sx2(n)/nt2)
    yy = sqrt(sy2(n)/nt2)
! covar. or correl. ---------------------
    sxy(n) = sxy(n)/nt2
!    sxy(n) = sxy(n)/nt2 / (xx*yy)
! ---------------------------------------
  enddo

  ! calculate eigen-function
  call symtrx(sxy,mij,d,z,mij,wk,ier)

  ! eigen-values (%)
  sumd = 0.
  do k=1, mij
    sumd = sumd + d(k)
  enddo
  do k=1, mij
    evalue100(k) = d(mij-k+1)*100./sumd
  enddo

  ! EOF
  do i=1, mij
  do k=1, nout
    vt(i,k) = z(i,mij-k+1)
  enddo
  enddo

  ij = 0
  do k=1, 4
  do j=1, ny
  do i=1, nx
    if (datanorm(i,j,1) .ne. missv) then
      ij = ij + 1
      eof(i,j,:,k) = vt(ij,:)
      if (normopt .eq. 1)  eof(i,j,:,k) = eof(i,j,:,k) / norm(i,j)
    else
      eof(i,j,:,k) = missv
    end if
  enddo
  enddo
  enddo

  ! PC
  pc = 0.
  do n=1, nt2*12/tinf(1), 12/tinf(1)
    ij = 0
    do k=0, 9/tinf(1), 3/tinf(1)
    do j=1, ny
    do i=1, nx
      if (datanorm(i,j,n+k) .ne. missv) then
        ij = ij + 1
        zc(ij)= datanorm(i,j,n+k) - zcm(ij)
!       zc(ij)= zc(ij)/sqrt(zz(ij))
      end if
    enddo
    enddo
    enddo

    nn = (n-1)*tinf(1)/12+1
    do k=1, nout
    do i=1, mij
      pc(nn,k) = pc(nn,k) + zc(i)*vt(i,k)
    enddo
    enddo
  enddo

  do k=1, nout
    com(k) = real(k)
    err100(k) = evalue100(k) * sqrt(2./nt2)
  enddo
  do k=1, 4
    tmp = tinf(2)+(k-1)*3
    if (tmp .gt. 12)  tmp = tmp - 12
    if (tmp .ge. 10) then
      write(varname(k),'(a,i2)') 'SEOF_m', tmp
    else
      write(varname(k),'(a,i1)') 'SEOF_m', tmp
    end if
  enddo
  do n=1, nt2
    t(n) = tinf(3) + (n-1)
  enddo

  ! output
  call out3d2d1d(trim(ofname)//'.nc',4,varname,eof, &
                 'x',nx,x,'y',ny,y,'com_eof',nout,com, &
                 1,(/'PC'/),pc,'t',nt2,t,'com_pc',nout,com, &
                 2,(/'e_value_100','err_100'/),(/evalue100(:nout),err100/), &
                 'com_ev',nout,com, 'S-EOF analysis')


  return

end subroutine seofnc


subroutine reofnc(ofname,nx,x,ny,y,nt,t,data3d,missn,missv,neof,eof, &
                  it_max,normopt,norm)
!-----------------------------------------------------------------------
!  Rotational EOF analysis
!-----------------------------------------------------------------------
!  INPUT
!       - neof   : no. of input EOFs which are used         [INT]
!       - eof(nx,ny,neof)                                   [REAL]
!                : input EOFs
!                  EOFs must be orthonormal.
!       - ev100(neof)                                       [REAL]
!                : input eigenvalues in percent
!       - it_max : maximum iteration no. for rotation       [INT]
!   
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, neof, it_max, normopt
  real,    intent(in) :: x(nx), y(ny), t(nt), data3d(nx,ny,nt)
  real,    intent(in) :: missv, eof(nx,ny,neof), norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, ij, k, km1
  integer :: mij, ier, tmpi, irpos(neof)
  real    :: normv(nx,ny), reof(nx,ny,neof), pc(nt,neof)
  real    :: reofs(nx,ny,neof), pcs(nt,neof)
  real    :: zcs(nx*ny-missn,nt), zcm(nx*ny-missn)
  real    :: com(neof), summ(neof), var100(neof)
  real    :: wk(nx*ny-missn), tmp1d(neof), tmp2d(neof,neof), tmp
  real, dimension(nx*ny-missn,neof) :: vt, rvt


  mij = nx*ny-missn

  if (normopt .eq. 1) then
    normv = norm
  else
    normv = 1.
  end if

  ij = 0
  do j=1, ny
  do i=1, nx
    if (data3d(i,j,1) .ne. missv) then
      ij = ij + 1
      zcs(ij,:) = data3d(i,j,:) * normv(i,j)
      vt(ij,:) = eof(i,j,:) * normv(i,j) !* sqrt(ev100(:))
    end if
  enddo
  enddo
  if (mij .ne. ij)  print*, ij, 'error (no. of missing values)'

  summ = 0.
  do k=1, neof
  do i=1, mij
    summ(k) = summ(k) + vt(i,k)**2
  enddo
  enddo
  print*, 'Input eigen-vectors must be orthonormal.'
  print*, sqrt(summ)

  zcm = 0.
  do ij=1, mij
    do n=1, nt
      zcm(ij) = zcm(ij) + zcs(ij,n)/nt
    enddo
    zcs(ij,:) = zcs(ij,:) - zcm(ij)
  enddo

  call ofrota(vt,mij,mij,neof,0,0,it_max,1.,0.0001,0.001, &
              rvt,mij,tmp2d,neof,tmp1d,wk,ier)

  summ = 0.
  do k=1, neof
    do ij=1, mij
      summ(k) = summ(k) + rvt(ij,k)**2
    enddo
    rvt(:,k) = rvt(:,k) / sqrt(summ(k))
  enddo

  ij = 0
  do j=1, ny
  do i=1, nx 
    if (data3d(i,j,1) .ne. missv) then 
      ij = ij + 1
      reof(i,j,:) = rvt(ij,:) / normv(i,j)
    else
      reof(i,j,:) = missv
    end if
  enddo
  enddo

  pc = 0.
  do n=1, nt
  do k=1, neof
  do ij=1, mij
    pc(n,k) = pc(n,k) + zcs(ij,n)*rvt(ij,k)
  enddo
  enddo
  enddo

  tmp = 0.
  var100 = 0.
  do k=1, neof
    do n=1, nt
      var100(k) = var100(k) + pc(n,k)**2
    enddo
    tmp = tmp + var100(k)
  enddo
  var100(:) = var100(:) / tmp * 100.

 
  do j=1, neof
    irpos(j) = j
  enddo
  ! sort -------------------- not be tested ----------
  do j=1, neof
  do k=neof, j+1, -1
    km1 = k - 1
    if (abs(var100(k)) .gt. abs(var100(km1))) then
      tmp         = var100(km1)
      var100(km1) = var100(k)
      var100(k)   = tmp
      tmpi       = irpos(km1)
      irpos(km1) = irpos(k)
      irpos(k)   = tmpi
    end if
  enddo
  enddo

  do k=1, neof
    reofs(:,:,k) = reof(:,:,irpos(k))
    pcs(:,k)     = pc(:,irpos(k))
  enddo
  ! --------------------------------------------------

  do k=1, neof
    com(k) = real(k)
  enddo

  call out3d2d1d(trim(ofname)//'.nc',1,(/'REOF'/),reofs, &
                 'x',nx,x,'y',ny,y,'com_reof',neof,com, &
                 1,(/'PC'/),pcs,'t',nt,t,'com_pc',neof,com, &
                 1,(/'e_value_100'/),var100,'com_ev',neof,com, &
                 'Rotated EOF analysis')


  return

end subroutine reofnc


subroutine ceofnc(ofname,nx,x,ny,y,nt,t,data3d,missn,missv,normopt,norm)
!-----------------------------------------------------------------------
!  Complex EOF analysis - not completely programed yet
!-----------------------------------------------------------------------

  use fft

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, normopt
  real,    intent(in) :: x(nx), y(ny), t(nt), data3d(nx,ny,nt)
  real,    intent(in) :: missv, norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, ij
  integer :: mij, ier
  real    :: sf(nx*ny-missn,nt), sfm(nx*ny-missn), hilbt(nt)
  real    :: normv(nx,ny), wei(nt)
  real    :: wk((nx*ny-missn)*2)
  complex :: iii, cf(nx*ny-missn,nt), covm(nx*ny-missn,nx*ny-missn)
  complex :: eval(nx*ny-missn), evec(nx*ny-missn,nx*ny-missn)
  double complex :: coef(nt)

  data       iii/(0.,1.)/

  mij = nx*ny-missn

  if (normopt .eq. 1) then
    normv = norm
  else
    normv = 1.
  end if

  ij = 0
  do j=1, ny
  do i=1, nx
    if (data3d(i,j,1) .ne. missv) then
      ij = ij + 1
      sf(ij,:) = data3d(i,j,:) * normv(i,j)
    end if
  enddo
  enddo
  if (mij .ne. ij)  print*, ij, 'error (no. of missing values)'

  sfm = 0.
  do ij=1, mij
    do n=1, nt
      sfm(ij) = sfm(ij) + sf(ij,n)/nt
    enddo
    sf(ij,:) = sf(ij,:) - sfm(ij)
  enddo

  wei = 1.
  do ij=1, mij
    call fft1df(nt,sf(ij,:),coef)
    coef(:) = (0.d0,1.d0)*coef(:) * wei(:)
    coef(nt/2+2:) = -coef(nt/2+2:)
    call fft1db(nt,coef,hilbt)

    cf(ij,:) = sf(ij,:) + iii*hilbt(:)
  enddo

  covm = (0.,0.)
  do i=1, mij
  do j=1, i
    do n=1, nt
      covm(i,j) = covm(i,j) + (real(cf(i,n))-iii*aimag(cf(i,n))) * cf(j,n)
    enddo
  enddo
  enddo
  covm = covm / nt
  do i=1, mij
  do j=i+1, mij
    covm(i,j) = real(covm(j,i)) - iii*aimag(covm(j,i))
  enddo
  enddo
print*, covm(1,1)
print*, covm(1,10), covm(10,1)

!  covm(1,:) = (/1.,1./)
!  covm(2,:) = (/1.,3./)

!  call eigrf(covm,2,2,1,eval,evec,2,wk,ier)

!  print*, eval
!  print*, evec

  call out1d(trim(ofname)//'.nc',2,(/'REOF','2'/),(/sf(1,:),hilbt/), &
               't',nt,t,'')

  return

end subroutine ceofnc


subroutine cond_eofnc(ofname,nx,x,ny,y,nt,t,data3d,missn,missv,ts, &
                      nout,normopt,norm)
!-----------------------------------------------------------------------
!  Conditional EOF analysis
!-----------------------------------------------------------------------
!  INPUT
!       - ts(nt)  : time series to removed                  [REAL]
!   
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nx, ny, nt, missn, nout, normopt
  real,    intent(in) :: x(nx), y(ny), t(nt), data3d(nx,ny,nt), ts(nt)
  real,    intent(in) :: missv, norm(nx,ny)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, tag
  real :: dat_rem(nx,ny,nt)
  real :: tsavg, varnt, covnt, tmp

  tsavg = 0.
  do n=1, nt
    tsavg = tsavg + ts(n)/nt
  enddo

  varnt = 0.
  do n=1, nt
    varnt = varnt + (ts(n)-tsavg)**2
  enddo

  ! remove the time series
  do j=1, ny
  do i=1, nx

    tag = 0
    do n=1, nt
      if (data3d(i,j,n) .eq. missv)  tag = 1
    enddo

    if (tag .eq. 0) then
      tmp = 0.
      do n=1, nt
        tmp = tmp + data3d(i,j,n)/nt
      enddo
    
      covnt = 0.
      do n=1, nt
        covnt = covnt + (data3d(i,j,n)-tmp)*(ts(n)-tsavg)
      enddo

      do n=1, nt
        dat_rem(i,j,n) = (data3d(i,j,n)-tmp) - (ts(n)-tsavg)*covnt/varnt
      enddo
    else
      dat_rem(i,j,:) = missv
    end if

  enddo
  enddo

  ! EOF 
  call eofnc(ofname,nx,x,ny,y,nt,t,dat_rem,missn,missv,nout,normopt,norm)

  call out1d(trim(ofname)//'_ts.nc',1,(/'TS'/),ts,'t',nt,t, &
             'time series of conditional EOF')


  return

end subroutine cond_eofnc


END module eofpack

