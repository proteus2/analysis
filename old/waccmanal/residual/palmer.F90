    subroutine elliptic_palmer(nny1,nnz1,n0sqin,frcin,chiout,duout)
!
    integer,                    intent(in) :: nny1  ,nnz1
    real, dimension(nnz1),      intent(in) :: n0sqin
    real, dimension(nny1,nnz1), intent(in) :: frcin
    real, dimension(nny1,nnz1), intent(out) :: chiout,duout
!   
!-----------------------------------------------------------------------
!   Coefficients of elliptic partial differential equation (cxxusr, cxyusr,
!   cyyusr, cxusr, cyusr, ceusr) are declared in module mud2crf90.
!-----------------------------------------------------------------------
!
    real, dimension(nny1,nnz1) :: forcing
    real, dimension(nny1,nnz1) :: vs
    real, dimension(nny1,nnz1) :: tmp   ,tmp1
!
    real(r8) :: fcor,tanova,cosf,fsign,dx,dy,eps
    real(r8), dimension(nny1)  :: ymr8
    real(r8), dimension(nnz1)  :: zmr8
    real(r8), dimension(nny1,nnz1) :: frcr8,chir8
!
    integer :: j     ,k
    integer :: j1    ,j2    ,k1    ,k2
    integer :: jj    ,kk    ,l
    integer :: npos
!
    integer  :: ixp,iex,jyq,jey
    integer  :: nxa,nxb,nyc,nyd
    integer  :: iguess,maxcy,method
    integer, dimension(4) :: mgopt
!
!   xusr, yusr, cxxusr, cxyusr, cyyusr, cxusr, cyusr, ceusr are declared in
!   mud2crf90, and they should be allocated here.
!
    allocate(xusr  (1:nny1))
    allocate(yusr  (1:nnz1))
    allocate(cxxusr(1:nny1,1:nnz1))
    allocate(cxyusr(1:nny1,1:nnz1))
    allocate(cyyusr(1:nny1,1:nnz1))
    allocate(cxusr (1:nny1,1:nnz1))
    allocate(cyusr (1:nny1,1:nnz1))
    allocate(ceusr (1:nny1,1:nnz1))

    do j=1,nny1
      xusr(j) = yp(j)
    end do
    do k=1,nnz1
      yusr(k) = zp(k)
    end do

    eps = 1.e-20
    do k=1,nnz1
      do j=1,nny1
        fcor   = 2._r8*omega*dsin(latp(j)*pi/180._r8)
        tanova = dtan(latp(j)*pi/180._r8)/arad
        cxxusr(j,k) = n0sqin(k)
        cxyusr(j,k) = eps
        cyyusr(j,k) = fcor*fcor
        cxusr (j,k) = n0sqin(k)*tanova
        cyusr (j,k) = -fcor*fcor/hscal
        ceusr (j,k) = 0._r8
      end do
    end do

!-----------------------------------------------------------------------
!   Determine forcing of the Elliptic equation
!   Initialize mass stream function
!-----------------------------------------------------------------------

    tmp(1:nny1,1:nnz1) = 0.0

    do k=1,nnz1-1
      do j=1,nny1
        tmp(j,k) = (2.0*omega*sin(latp(j)*pi/180.))*  &
                    (frcin(j,k+1)-frcin(j,k))/(zp(k+1)-zp(k))*  &
                     cos(latp(j)*pi/180.)
      end do
    end do

    do k=2,nnz1-1
      do j=1,nny1
        forcing(j,k) = (tmp(j,k)+tmp(j,k-1))*0.5
      end do
    end do
    
    do j=1,nny1
      forcing(j,1) = forcing(j,2)
      forcing(j,nnz1) = forcing(j,nnz-1)
    end do

    do k=1,nnz1
      do j=1,nny1
        forcing(j,k) = forcing(j,k)
      end do   
    end do 

    do k=1,nnz1
      do j=1,nny1
        chiout(j,k) = 0.0
      end do 
    end do 

    do k=1,nnz1
      zmr8(k) = zm(k)*1._r8
    end do
    do j=1,nny1
      ymr8(j) = ym(j)*1._r8
    end do
    do k=1,nnz1
      do j=1,nny1
        frcr8(j,k) = forcing(j,k)*1._r8
        chir8(j,k) = chiout(j,k)*1._r8
      end do
    end do

!-----------------------------------------------------------------------
!   Solve Elliptic equation using MULTIGRID PACKAGE
!-----------------------------------------------------------------------

    ixp = 2; jyq = 3; iex = 6; jey = 5
    nxa = 1; nxb = 1; nyc = 1; nyd = 2
    iguess = 0; maxcy = 1; method = 0
    mgopt = (/2,2,1,3/)

    call sol2cr(ixp,jyq,iex,jey,nny1,nnz1,  &
                nxa,nxb,nyc,nyd,iguess,maxcy,method,  &
                mgopt,ymr8,zmr8,frcr8,chir8,20)

    do k=1,nnz1
      do j=1,nny1
        if ( abs(chir8(j,k)) < 1.e-20 ) then
          chir8(j,k) = 0.0
        end if
        chiout(j,k) = chir8(j,k)
      end do
    end do

    deallocate(xusr)
    deallocate(yusr)
    deallocate(cxxusr)
    deallocate(cxyusr)
    deallocate(cyyusr)
    deallocate(cxusr)
    deallocate(cyusr)
    deallocate(ceusr)

!-----------------------------------------------------------------------
!   meridional flow
!-----------------------------------------------------------------------

    tmp (1:nny1,1:nnz1) = 0.0
    tmp1(1:nny1,1:nnz1) = 0.0
    vs  (1:nny1,1:nnz1) = 0.0

    do k=1,nnz1-1
    do j=1,nny1
      tmp(j,k) = (chiout(j,k+1)-chiout(j,k))/(zp(k+1)-zp(k))
    end do
    end do

    do k=2,nnz1-1
    do j=1,nny1
      tmp1(j,k) = (tmp(j,k)+tmp(j,k-1))*0.5
    end do
    end do

    do k=2,nnz1-1
    do j=1,nny1
      vs (j,k) = (-1./cos(latp(j)*pi/180.))*(tmp1(j,k)-chiout(j,k)/hscal)
    end do
    end do   

    do k=2,nnz1-1
    do j=1,nny1
      duout(j,k) = (2.0*omega*sin(latp(j)*pi/180.))*vs(j,k)+forcing(j,k)
    end do
    end do

    do j=1,nny1
      duout(j,1)    = 1.e+20
      duout(j,nnz1) = 1.e+20
    end do

    return
    end subroutine elliptic_palmer
