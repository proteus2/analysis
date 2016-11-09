MODULE svdpack
!-----------------------------------------------------------------------
!  Routines for SVD
!-----------------------------------------------------------------------
!  OUTPUT : NetCDF format
!
!  INPUT
!       - ofname   : output file name                       [CHAR]
!       - nt       : no. of time dimension                  [INT]
!       - t(nt)    : time                                   [REAL]
!       - nx?      : no. of the 1st dimension of dat3d?     [INT]
!       - x?(nx?)  : axis of the 1st dimension of dat3d?    [REAL]
!       - ny?      : no. of the 2nd dimension of dat3d?     [INT]
!       - y?(ny?)  : axis of the 2nd dimension of dat3d?    [REAL]
!       - dat3d1(nx1,ny1,nt)                                [REAL]
!                  : 1st data used to analysis 
!                    It is decomposed into left singular vectors.
!       - dat3d2(nx2,ny2,nt)                                [REAL]
!                  : 2nd data used to analysis
!                    It is decomposed into right singular vectors.
!       - missn?   : no. of missing values of dat3d?        [INT]
!                    in x,y domain
!       - missv?   : missing value of dat3d?                [REAL]
!       - normopt? : option for normalization of dat3d?     [INT]
!                   == 1, data is normalized by norm?(x?,y?)
!       - norm?(nx?,ny?)                                    [REAL]
!                  : normalization factor for dat3d?
!       - coropt   : option for matrix to be used           [INT]
!                   == 1, correlation matrix
!                   /= 1, covariance matrix
!       - nout     : no. of component for output            [INT]
!
!  NOTICE
!       - (nx1 * ny1 - missn1) must be larger than (nx2 * ny2 - missn2).
!       - coropt is not valid yet.
!
!-----------------------------------------------------------------------

  use netcdfio

  public :: svdnc, jsvdnc

  contains

subroutine svdnc(ofname,nt,t, &
                 nx1,x1,ny1,y1,dat3d1,missn1,missv1,normopt1,norm1, &
                 nx2,x2,ny2,y2,dat3d2,missn2,missv2,normopt2,norm2, &
                 coropt,nout)
!-----------------------------------------------------------------------
!  SVD
!-----------------------------------------------------------------------


  implicit none

  integer, intent(in) :: nt, nx1, ny1, nx2, ny2, missn1, missn2
  integer, intent(in) :: normopt1, normopt2, coropt, nout
  real,    intent(in) :: t(nt), x1(nx1), y1(ny1), x2(nx2), y2(ny2)
  real,    intent(in) :: dat3d1(nx1,ny1,nt), dat3d2(nx2,ny2,nt)
  real,    intent(in) :: missv1, missv2, norm1(nx1,ny1), norm2(nx2,ny2)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, nn, ij, k
  integer :: mij1, mij2, ier, tag, jp1, km1, tmpi
  integer :: irpos(nx2*ny2-missn2)

  real :: normv1(nx1,ny1), normv2(nx2,ny2)
  real :: f1(nx1*ny1-missn1,nt), f2(nx2*ny2-missn2,nt)
  real :: var1(nx1*ny1-missn1), var2(nx2*ny2-missn2)
  real :: covm(nx1*ny1-missn1,nx2*ny2-missn2)
  real :: tmp, tmp1(nx1*ny1-missn1), tmp2(nx2*ny2-missn2)
  real :: sv(nx2*ny2-missn2), resm(nx2*ny2-missn2,nx2*ny2-missn2)
  real :: scf100(nx2*ny2-missn2)
  real :: svecl(nx1*ny1-missn1,nout), svecr(nx2*ny2-missn2,nout)
  real :: svec1(nx1,ny1,nout), svec2(nx2,ny2,nout)
  real :: sumvar1, sumvar2, com(nout)

  real, dimension(nt,nout) :: pc1, pc2
  real, dimension(nout)    :: var1k100, var2k100, corr


  mij1 = nx1*ny1-missn1
  mij2 = nx2*ny2-missn2

  ! normalize
  normv1 = 1.
  normv2 = 1.
  if (normopt1 .eq. 1)  normv1 = norm1
  if (normopt2 .eq. 1)  normv2 = norm2

  ij = 0
  do j=1, ny1
  do i=1, nx1
    tag = 0
    do n=1, nt
      if (dat3d1(i,j,n) .eq. missv1)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      f1(ij,:) = dat3d1(i,j,:) * normv1(i,j)
    end if
  enddo
  enddo
  if (mij1 .ne. ij)  print*, ij, 'data1 error (no. of missing values)'

  ij = 0
  do j=1, ny2
  do i=1, nx2
    tag = 0
    do n=1, nt 
      if (dat3d2(i,j,n) .eq. missv2)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      f2(ij,:) = dat3d2(i,j,:) * normv2(i,j)
    end if
  enddo
  enddo
  if (mij2 .ne. ij)  print*, ij, 'data2 error (no. of missing values)'

  ! remove mean
  tmp1 = 0.  ;  tmp2 = 0.
  do n=1, nt
    do ij=1, mij1
      tmp1(ij) = tmp1(ij) + f1(ij,n)
    enddo
    do ij=1, mij2
      tmp2(ij) = tmp2(ij) + f2(ij,n)
    enddo
  enddo
  do n=1, nt
    f1(:,n) = f1(:,n) - tmp1(:)/nt
    f2(:,n) = f2(:,n) - tmp2(:)/nt
  enddo

  ! covariance matrix
  var1 = 0.  ;  var2 = 0.  ;  covm = 0.
  do n=1, nt
    do ij=1, mij1
      var1(ij) = var1(ij) + f1(ij,n)**2
    enddo
    do ij=1, mij2
      var2(ij) = var2(ij) + f2(ij,n)**2
    enddo
    do j=1, mij2
    do i=1, mij1
      covm(i,j) = covm(i,j) + f1(i,n)*f2(j,n)
    enddo
    enddo
  enddo
  var1 = var1 / nt
  var2 = var2 / nt
  covm = covm / nt

  ! correlation matrix
  if (coropt .eq. 1) then
    do j=1, mij2
    do i=1, mij1
      covm(i,j) = covm(i,j) / sqrt(var1(i)*var2(j))
    enddo
    enddo
  end if

  ! SVD
  print*, ' Start SVDCMP'
  call svdcmp(covm,mij1,mij2,mij1,mij2,sv,resm)
  print*, ' End SVDCMP'
  scf100(:) = sv(:)**2
  tmp = 0.
  do j=1, mij2
    tmp = tmp + scf100(j)
  enddo
  scf100 = scf100 / tmp * 100.

  print*, scf100(:100)
  do j=1, mij2
    irpos(j) = j
  enddo
  ! sort -------------------- not be tested ----------
  do j=1, nout
    jp1 = j + 1
    do k=mij2, jp1, -1
      km1 = k - 1
      if (abs(scf100(k)) .gt. abs(scf100(km1))) then
        tmp         = scf100(km1)
        scf100(km1) = scf100(k)
        scf100(k)   = tmp
        tmpi       = irpos(km1)
        irpos(km1) = irpos(k)
        irpos(k)   = tmpi
      end if
    enddo
  enddo
  ! --------------------------------------------------

  do k=1, nout
    print*, k, irpos(k)
    do ij=1, mij1
      svecl(ij,k) = covm(ij,irpos(k))
    enddo
    do ij=1, mij2
      svecr(ij,k) = resm(ij,irpos(k))
    enddo
  enddo

  ! PC
  pc1 = 0.  ;  pc2 = 0.
  var1k100 = 0.  ;  var2k100 = 0.
  do n=1, nt
  do k=1, nout
    do ij=1, mij1
      pc1(n,k) = pc1(n,k) + svecl(ij,k)*f1(ij,n)
    enddo
    do ij=1, mij2
      pc2(n,k) = pc2(n,k) + svecr(ij,k)*f2(ij,n)
    enddo
    var1k100(k) = var1k100(k) + pc1(n,k)**2
    var2k100(k) = var2k100(k) + pc2(n,k)**2
  enddo
  enddo
  var1k100 = var1k100 / nt
  var2k100 = var2k100 / nt

  ! PC correlation
  corr = 0.
  do k=1, nout
    do n=1, nt
      corr(k) = corr(k) + pc1(n,k)*pc2(n,k)
    enddo
    corr(k) = corr(k)/nt / sqrt(var1k100(k)*var2k100(k))
  enddo

  ! variance of each mode (%)
  sumvar1 = 0.  ;  sumvar2 = 0.
  do ij=1, mij1
    sumvar1 = sumvar1 + var1(ij)
  enddo
  do ij=1, mij2
    sumvar2 = sumvar2 + var2(ij)
  enddo
  var1k100 = var1k100 / sumvar1 * 100.
  var2k100 = var2k100 / sumvar2 * 100.

  ! singular vectors
  ij = 0
  do j=1, ny1
  do i=1, nx1
    tag = 0
    do n=1, nt
      if (dat3d1(i,j,n) .eq. missv1)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      svec1(i,j,:) = svecl(ij,:) / normv1(i,j)
    else
      svec1(i,j,:) = missv1
    end if
  enddo
  enddo
  ij = 0
  do j=1, ny2
  do i=1, nx2
    tag = 0
    do n=1, nt
      if (dat3d2(i,j,n) .eq. missv2)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      svec2(i,j,:) = svecr(ij,:) / normv2(i,j)
    else
      svec2(i,j,:) = missv2
    end if
  enddo
  enddo

  ! output
  do k=1, nout
    com(k) = real(k)
  enddo

  call out3d2d1d(trim(ofname)//'_field1.nc',1,(/'EOF'/),svec1, &
                 'x',nx1,x1,'y',ny1,y1,'com_eof',nout,com,     &
                 1,(/'PC'/),pc1,'t',nt,t,'com_pc',nout,com,    &
                 1,(/'var100'/),var1k100,'com_var',nout,com,   &
                 'SVD EOF field1')
  call out3d2d1d(trim(ofname)//'_field2.nc',1,(/'EOF'/),svec2, &
                 'x',nx2,x2,'y',ny2,y2,'com_eof',nout,com,     &
                 1,(/'PC'/),pc2,'t',nt,t,'com_pc',nout,com,    &
                 1,(/'var100'/),var2k100,'com_var',nout,com,   &
                 'SVD EOF field2')
  call out1d(trim(ofname)//'_scf.nc',2,(/'SCF100','corr'/),    &
             (/scf100(:nout),corr/),'com',nout,com,            &
             'SVD squared covariance fraction and corr.')


  return

end subroutine svdnc


subroutine jsvdnc(ofname,nt,t, &
              nx1a,x1a,ny1a,y1a,dat3d1a,missn1a,missv1a,normopt1a,norm1a,&
              nx1b,x1b,ny1b,y1b,dat3d1b,missn1b,missv1b,normopt1b,norm1b,&
              nx2,x2,ny2,y2,dat3d2,missn2,missv2,normopt2,norm2, &
              coropt,nout)
!-----------------------------------------------------------------------
!  Joint SVD
!-----------------------------------------------------------------------
!  INPUT
!       - dat3d1a and dat3d1b                               [REAL]
!                  : input datum to be jointed
!
!  NOTICE
!       - (nx1a * ny1a - missn1a + nx1b * ny1b - missn1b)
!         must be larger than (nx2 * ny2 - missn2).
!       - coropt is not valid yet.
!
!-----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nt, nx1a, ny1a, missn1a, nx1b, ny1b, missn1b
  integer, intent(in) :: nx2, ny2, missn2
  integer, intent(in) :: normopt1a, normopt1b, normopt2, coropt, nout
  real,    intent(in) :: t(nt), x1a(nx1a), y1a(ny1a), x1b(nx1b), y1b(ny1b)
  real,    intent(in) :: x2(nx2), y2(ny2)
  real,    intent(in) :: dat3d1a(nx1a,ny1a,nt), dat3d1b(nx1b,ny1b,nt)
  real,    intent(in) :: dat3d2(nx2,ny2,nt)
  real,    intent(in) :: missv1a, missv1b, missv2
  real,    intent(in) :: norm1a(nx1a,ny1a), norm1b(nx1b,ny1b), norm2(nx2,ny2)
  character(len=*), intent(in) :: ofname

  integer :: i,j,n, nn, ij, k
  integer :: mij1, mij1a, mij1b, mij2, ier, tag, jp1, km1, tmpi
  integer :: irpos(nx2*ny2-missn2)

  real :: normv1a(nx1a,ny1a), normv1b(nx1b,ny1b), normv2(nx2,ny2)
  real :: f1(nx1a*ny1a-missn1a+nx1b*ny1b-missn1b,nt)
  real :: f2(nx2*ny2-missn2,nt)
  real :: var1(nx1a*ny1a-missn1a+nx1b*ny1b-missn1b)
  real :: var2(nx2*ny2-missn2)
  real :: covm(nx1a*ny1a-missn1a+nx1b*ny1b-missn1b,nx2*ny2-missn2)
  real :: tmp1(nx1a*ny1a-missn1a+nx1b*ny1b-missn1b)
  real :: tmp2(nx2*ny2-missn2), tmp
  real :: sv(nx2*ny2-missn2), resm(nx2*ny2-missn2,nx2*ny2-missn2)
  real :: scf100(nx2*ny2-missn2)
  real :: svecl(nx1a*ny1a-missn1a+nx1b*ny1b-missn1b,nout)
  real :: svecr(nx2*ny2-missn2,nout)
  real :: svec1a(nx1a,ny1a,nout), svec1b(nx1b,ny1b,nout)
  real :: svec2(nx2,ny2,nout)
  real :: sumvar1, sumvar2, com(nout)

  real, dimension(nt,nout) :: pc1, pc2
  real, dimension(nout)    :: var1k100, var2k100, corr


  mij1a = nx1a*ny1a-missn1a
  mij1b = nx1b*ny1b-missn1b
  mij1 = mij1a + mij1b
  mij2 = nx2*ny2-missn2

  normv1a = 1.
  normv1b = 1.
  normv2  = 1.
  if (normopt1a .eq. 1)  normv1a = norm1a
  if (normopt1b .eq. 1)  normv1b = norm1b
  if (normopt2 .eq. 1)   normv2  = norm2


  ij = 0
  do j=1, ny1a
  do i=1, nx1a
    tag = 0
    do n=1, nt
      if (dat3d1a(i,j,n) .eq. missv1a)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      f1(ij,:) = dat3d1a(i,j,:) * normv1a(i,j)
    end if
  enddo
  enddo
  if (mij1a .ne. ij)  print*, ij, 'data1a error (no. of missing values)'
  do j=1, ny1b
  do i=1, nx1b
    tag = 0
    do n=1, nt
      if (dat3d1b(i,j,n) .eq. missv1b)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      f1(ij,:) = dat3d1b(i,j,:) * normv1b(i,j)
    end if
  enddo
  enddo
  if (mij1 .ne. ij)  print*, ij, 'data1 error (no. of missing values)'

  ij = 0
  do j=1, ny2
  do i=1, nx2
    tag = 0
    do n=1, nt
      if (dat3d2(i,j,n) .eq. missv2)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      f2(ij,:) = dat3d2(i,j,:) * normv2(i,j)
    end if
  enddo
  enddo
  if (mij2 .ne. ij)  print*, ij, 'data2 error (no. of missing values)'

  ! remove mean
  tmp1 = 0.  ;  tmp2 = 0.
  do n=1, nt
    do ij=1, mij1
      tmp1(ij) = tmp1(ij) + f1(ij,n)
    enddo
    do ij=1, mij2
      tmp2(ij) = tmp2(ij) + f2(ij,n)
    enddo
  enddo
  do n=1, nt
    f1(:,n) = f1(:,n) - tmp1(:)/nt
    f2(:,n) = f2(:,n) - tmp2(:)/nt
  enddo

  ! covariance matrix
  var1 = 0.  ;  var2 = 0.  ;  covm = 0.
  do n=1, nt
    do ij=1, mij1
      var1(ij) = var1(ij) + f1(ij,n)**2
    enddo
    do ij=1, mij2
      var2(ij) = var2(ij) + f2(ij,n)**2
    enddo
    do j=1, mij2
    do i=1, mij1
      covm(i,j) = covm(i,j) + f1(i,n)*f2(j,n)
    enddo
    enddo
  enddo
  var1 = var1 / nt
  var2 = var2 / nt
  covm = covm / nt

  ! correlation matrix
  if (coropt .eq. 1) then
    do j=1, mij2
    do i=1, mij1
      covm(i,j) = covm(i,j) / sqrt(var1(i)*var2(j))
    enddo
    enddo
  end if
  
  ! SVD
  print*, ' Start SVDCMP'
  call svdcmp(covm,mij1,mij2,mij1,mij2,sv,resm)
  print*, ' End SVDCMP'
  scf100(:) = sv(:)**2
  tmp = 0.
  do j=1, mij2
    tmp = tmp + scf100(j)
  enddo
  scf100 = scf100 / tmp * 100.
  
  print*, scf100(:100)
  do j=1, mij2 
    irpos(j) = j
  enddo
  ! sort -------------------- not be tested ----------
  do j=1, nout
    jp1 = j + 1
    do k=mij2, jp1, -1
      km1 = k - 1
      if (abs(scf100(k)) .gt. abs(scf100(km1))) then
        tmp         = scf100(km1)
        scf100(km1) = scf100(k)
        scf100(k)   = tmp
        tmpi       = irpos(km1)
        irpos(km1) = irpos(k)
        irpos(k)   = tmpi
      end if
    enddo
  enddo
  ! --------------------------------------------------
  
  do k=1, nout
    print*, k, irpos(k)
    do ij=1, mij1
      svecl(ij,k) = covm(ij,irpos(k))
    enddo 
    do ij=1, mij2
      svecr(ij,k) = resm(ij,irpos(k))
    enddo 
  enddo

  ! PC
  pc1 = 0.  ;  pc2 = 0.
  var1k100 = 0.  ;  var2k100 = 0.
  do n=1, nt
  do k=1, nout
    do ij=1, mij1
      pc1(n,k) = pc1(n,k) + svecl(ij,k)*f1(ij,n)
    enddo
    do ij=1, mij2
      pc2(n,k) = pc2(n,k) + svecr(ij,k)*f2(ij,n)
    enddo
    var1k100(k) = var1k100(k) + pc1(n,k)**2
    var2k100(k) = var2k100(k) + pc2(n,k)**2
  enddo
  enddo 
  var1k100 = var1k100 / nt
  var2k100 = var2k100 / nt

  ! PC correlation
  corr = 0.
  do k=1, nout
    do n=1, nt
      corr(k) = corr(k) + pc1(n,k)*pc2(n,k)
    enddo
    corr(k) = corr(k)/nt / sqrt(var1k100(k)*var2k100(k))
  enddo
  
  ! variance of each mode (%)
  sumvar1 = 0.  ;  sumvar2 = 0.
  do ij=1, mij1
    sumvar1 = sumvar1 + var1(ij)
  enddo
  do ij=1, mij2
    sumvar2 = sumvar2 + var2(ij)
  enddo
  var1k100 = var1k100 / sumvar1 * 100.
  var2k100 = var2k100 / sumvar2 * 100.

  ! singular vectors
  ij = 0
  do j=1, ny1a
  do i=1, nx1a
    tag = 0
    do n=1, nt
      if (dat3d1a(i,j,n) .eq. missv1a)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      svec1a(i,j,:) = svecl(ij,:) / normv1a(i,j)
    else
      svec1a(i,j,:) = missv1a
    end if
  enddo
  enddo
  do j=1, ny1b
  do i=1, nx1b
    tag = 0
    do n=1, nt
      if (dat3d1b(i,j,n) .eq. missv1b)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      svec1b(i,j,:) = svecl(ij,:) / normv1b(i,j)
    else
      svec1b(i,j,:) = missv1b
    end if
  enddo
  enddo

  ij = 0
  do j=1, ny2
  do i=1, nx2
    tag = 0
    do n=1, nt
      if (dat3d2(i,j,n) .eq. missv2)  tag = 1
    enddo
    if (tag .ne. 1) then
      ij = ij + 1
      svec2(i,j,:) = svecr(ij,:) / normv2(i,j)
    else
      svec2(i,j,:) = missv2
    end if
  enddo
  enddo

  ! output
  do k=1, nout
    com(k) = real(k)
  enddo

  call out3d2d1d(trim(ofname)//'_field1a.nc',1,(/'EOF'/),svec1a, &
                 'x',nx1a,x1a,'y',ny1a,y1a,'com_eof',nout,com,     &
                 1,(/'PC'/),pc1,'t',nt,t,'com_pc',nout,com,    &
                 1,(/'var100'/),var1k100,'com_var',nout,com,   &
                 'J-SVD EOF field1a')
  call out3d(trim(ofname)//'_field1b.nc',1,(/'EOF'/),svec1b, &
             'x',nx1b,x1b,'y',ny1b,y1b,'com_eof',nout,com,     &
             'J-SVD EOF field1b')
  call out3d2d1d(trim(ofname)//'_field2.nc',1,(/'EOF'/),svec2, &
                 'x',nx2,x2,'y',ny2,y2,'com_eof',nout,com,     &
                 1,(/'PC'/),pc2,'t',nt,t,'com_pc',nout,com,    &
                 1,(/'var100'/),var2k100,'com_var',nout,com,   &
                 'J-SVD EOF field2')
  call out1d(trim(ofname)//'_scf.nc',2,(/'SCF100','corr'/),    &
             (/scf100(:nout),corr/),'com',nout,com,            &
             'J-SVD squared covariance fraction and corr.')


  return

end subroutine jsvdnc


END module svdpack

