!2345
    module meanwavnum
!----------------------------------------------------------------------------------
!   
!   Calculate mean wave numbers (i.e., m and k) using the dispersion relation
!   of inertia gravity waves
!
!----------------------------------------------------------------------------------

    use window

    contains

    subroutine meanmk(nz,dz,f,wf,nbas,ustprt,vstprt,mbar  ,mbaran,kbar   ,kbaran)
 
    implicit none

    integer,             intent(in)  :: nz
    real,                intent(in)  :: dz    ,f     ,wf     ,nbas
    real, dimension(nz), intent(in)  :: ustprt,vstprt
    real,                intent(out) :: mbar  ,mbaran,kbar   ,kbaran

    integer                          :: k
    real                             :: pi    ,s1    ,s2
    real                             :: umean ,vmean ,wss  
    real, dimension(nz)              :: winfcn
    real, dimension(nz)              :: us    ,vs    ,usvs
    real, dimension(nz/2+1)          :: m 
    double precision, dimension(4*nz+15) :: work
    double complex                   :: ci
    double complex, dimension(nz)    :: seq, coef 

    pi = 4.0*atan(1.0)
 
    do k=1,nz/2+1
      m(k) = float(k-1)/(nz*dz)
    end do

    call welch(nz,winfcn,wss)

    umean = 0.0
    vmean = 0.0
    do k=1,nz
      umean = umean + ustprt(k)
      vmean = vmean + vstprt(k)
    end do
    umean = umean/float(nz)
    vmean = vmean/float(nz)

    do k=1,nz
      seq(k) = dble((ustprt(k)-umean)*winfcn(k)) + ci*0.d0
      coef(k) = seq(k)
    end do

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,coef,work)

    do k=1,nz
      us(k) = real(cdabs(coef(k)*conjg(coef(k))))
    end do

    do k=1,nz
      seq(k) = dble((vstprt(k)-vmean)*winfcn(k)) + ci*0.d0
      coef(k) = seq(k)
    end do

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,coef,work)

    do k=1,nz
      vs(k) = real(cdabs(coef(k)*conjg(coef(k))))
    end do
 
    do k=1,nz
      usvs(k) = (us(k)+vs(k))*dz/wss
    end do

    s1 = 0.0
    s2 = 0.0
    do k=2,nz/2
      s1 = s1 + usvs(k)*m(k)
      s2 = s2 + usvs(k)
    end do

    mbar   = s1/s2
    mbaran = 2.0*pi*mbar
    kbaran = f*mbaran*sqrt(wf**2-1)/nbas
    kbar   = kbaran/(2.0*pi)

    return
 
    end subroutine meanmk

    end module meanwavnum

    

