!2345
    module pwrspd

    use window

    contains 

    subroutine psd(ndata,ma,npc,x,dx,data,wavnum,pwsd)

!-------------------------------------------------------------------------------
!
!   Calculate power spectrum density.
!
!   1. Mean of given time series is removed
!   2. Trend is removed in the time series
!   3. Welch window function is applied to the modified time series to prevent
!      energy leakage in each wavenumber bin
!   4. Calculate power spectral density
!
!   This program needs cfftpack.f, Double complex Fast Fourier Transform routine  
!   lfit, covsrt, and gaussj subroutines are imported from Numerical Recipes.  
!
!   Author  
!   Song, In-Sun
!
!   Date
!   November 19, 2001
!
!   Usage
!
!   INPUT DATA
!
!   ndata : The number of elements of time series vector to be analysized.
!   ma    : The number of coefficients of polynomial to be used to remove
!           a trend in time series
!           ma = 2 : Linear trend
!           ma = 3 : Second order polynomial trend
!           ma cannot exceed 50
!   npc   : The number of points of covariance (npcxnpc) matrix
!           npc should be larger than or equal to ma
!   x     : Time ( or distance )
!   dx    : dt ( or dx )
!   data  : The time series vector
!
!   OUTPUT DATA
!
!   wavnum : Wavenumber vector      : wavnum(ndata/2+1)
!   pwsd   : Power spectral density : pwsd(ndata/2+1)
!
!-------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: ndata,ma,npc
    real, intent(in) :: dx
    real, dimension(ndata), intent(inout) :: x, data
    real, dimension(ndata/2+1), intent(out) :: wavnum,pwsd 

    integer :: i,j
    integer, dimension(ma) :: ia
    real :: chisq
    real :: sum, an, bn, wss, sqr, sqf
    real, dimension(ma) :: a
    real, dimension(ndata) :: data1, sig, func, trend
    real, dimension(npc,npc) :: covar
    real, dimension(50) :: afunc
    real, dimension(ndata) :: windowfcn

    double precision :: pi
    double precision, dimension(4*ndata+15) :: wsave
    double complex :: ci
    double complex, dimension(ndata) :: seq, coef
   
    pi = 3.141592d0
    ci = (0.d0,1.d0) 
       
    sum = 0.0
    do i=1,ndata
     sum  = sum  + data(i)
    end do 
    sum = sum/float(ndata)

    an = ( ndata - 1.0 )/2.
    bn = ( ndata + 1.0 )/2.

    data1 = data

    call welch(ndata,windowfcn,wss)

    do i=1,ndata
      data1(i) = data1(i)  - sum
      sig(i) = 1.0
      trend(i) = 0.0
    end do 

    do i=1,ma
     ia(i) = 1
      a(i) = 1.0
    end do 
   
    call lfit(x,data1,sig,ndata,a,ia,ma,covar,npc,chisq)

    do i=1,ndata
      call funcs(x(i),afunc,ma)
      do j=1,ma
        trend(i) = trend(i) + a(j)*afunc(j)
      end do 
      data1(i) = data1(i) - trend(i)
    end do 

    do i=1,ndata
      data1(i) = data1(i)*windowfcn(i)
      seq(i) = dble(data1(i)) + ci*0.0 
    end do 

    coef(:) = seq(:)

    call cffti(ndata,wsave)
    call cfftf(ndata,coef,wsave)

    pwsd(1) = (dx/wss)*real(cdabs(coef(1))**2)
    do i=2,ndata/2
      pwsd(i) = (dx/wss)*real(cdabs(coef(i))**2 +    & 
                              cdabs(coef(ndata-i+2))**2 )
    end do 
    pwsd(ndata/2+1) = 2*(dx/wss)*real(cdabs(coef(ndata/2+1))**2)

    sqr = 0.
    do i=1,ndata
      sqr = sqr + real(cdabs(seq(i))**2)*dx
    end do
    sqf = 0.
    do i=1,ndata
      sqf = sqf + (dx/ndata)*real(cdabs(coef(i))**2)
    end do

    write(6,'(a,e17.10)') 'SQR = ',sqr
    write(6,'(a,e17.10)') 'SQF = ',sqf

    do i=1,ndata/2+1
      wavnum(i) = (float(i-1)/ndata)*(1./dx)
    end do

    return
    end subroutine psd 

!-------------------------------------------------------------------------------
!
!   SUBROUTINE SMOOTH
!
!-------------------------------------------------------------------------------

    subroutine smooth(n   ,wavn,psd ,psdm)

    integer,            intent(in)  :: n
    real, dimension(n), intent(in)  :: wavn, psd
    real, dimension(n), intent(out) :: psdm

    integer                         :: k    ,l
    integer                         :: nreg 
    real                            :: mnlg ,mxlg
    real                            :: mnwn ,mxwn
    real                            :: x0   ,x1   ,x2
    real                            :: sum1 ,sum2 ,sum3
    real, dimension(n)              :: wavnlog
    real, allocatable, dimension(:) :: xreg ,psdreg

    do k=1,n
      if ( wavn(k) == 0.0 ) then
        write(6,'(a)') '(SMOOTH): Wavenumber should not have 0.'
        write(6,'(a)') '(SMOOTH): Program stops'
        stop
      end if
    end do

    do k=1,n-1
      if ( wavn(k+1) < wavn(k) ) then
        write(6,'(a)') '(SMOOTH): Wavenumber should monotonically increase.'
        write(6,'(a)') '(SMOOTH): Program stops.'
        stop
      end if
    end do
      
    wavnlog = log10(wavn)
    mnlg    = wavnlog(1)
    mxlg    = wavnlog(n)

!-------------------------------------------------------------------------------
!
!   In case of mnlg = -3.848, 
!   minw = -3.0 + (-3.848-(-3))/0.1 - 0.1 = -3.0 - 0.8 - 0.1 = -3.9
!
!   In case of mxlg = -2.323,
!   maxw = -2.0 + (-2.323-(-2))/0.1 + 0.1 = -2.0 - 0.3 + 0.1 = -2.2
!
!-------------------------------------------------------------------------------

    mnwn = float(int(mnlg)) + int((mnlg-int(mnlg))/0.1)*0.1 - 0.1
    mxwn = float(int(mxlg)) + int((mxlg-int(mxlg))/0.1)*0.1 + 0.1

    nreg = int((mxwn-mnwn)/0.1) + 1

    allocate(xreg  (1:nreg)) 
    allocate(psdreg(1:nreg)) 

    do k=1,nreg
      psdreg(k) = 0.0
    end do

    do k=1,nreg
      xreg(k) = mnwn + float(k-1)*0.1 
    end do

    do k=1,nreg
      do l=1,n
        if ( wavnlog(l) >= xreg(k) - 0.05 .and. wavnlog(l) < xreg(k) + 0.05 ) then
          psdreg(k) = psdreg(k) + psd(l)
        end if        
      end do
    end do

!-------------------------------------------------------------------------------
!   Back to original wavenumber grid
!-------------------------------------------------------------------------------

    do l=1,n
      psdm(l) = 0.0
    end do

    do l=1,n
      do k=1,nreg-1
        x2 = xreg(k+1)
        x1 = xreg(k)
        x0 = wavnlog(l)
        if ( x0 >= x1 .and. x0 < x2 ) then
          psdm(l) = psdreg(k)*(x2-x0)/(x2-x1) + psdreg(k+1)*(x0-x1)/(x2-x1)
        end if
      end do
    end do

!-------------------------------------------------------------------------------
!   Check consistency
!-------------------------------------------------------------------------------

    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0

    do k=1,n
      sum1 = sum1 + psd(k)
      sum2 = sum2 + psdm(k)
    end do 
    do k=1,nreg
      sum3 = sum3 + psdreg(k)
    end do

    do k=1,n
      psdm(k) = psdm(k)/sum2*sum1
    end do

    write(6,'(a,f12.5)') '(SMOOTH): SUM OF ORIGINAL PSD : ',sum1
    write(6,'(a,f12.5)') '(SMOOTH): SUM OF INTERPOL PSD : ',sum2
    write(6,'(a,f12.5)') '(SMOOTH): SUM OF SMOOTHED PSD : ',sum3

    deallocate(xreg)
    deallocate(psdreg)

    return

    end subroutine smooth
    
!-------------------------------------------------------------------------------
!
!   SUBROUTINE LFIT
!
!-------------------------------------------------------------------------------

    subroutine lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq)
      
    integer max,ia(ma),npc,ndat,MMAX
    real chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
    parameter (MMAX=50)
    integer i,j,k,l,m,mfit
    real sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)
    mfit=0
       
    do j=1,ma
      if (ia(j).ne.0) mfit=mfit+1
    enddo

    if(mfit.eq.0) pause 'lfit: no parameters to be fitted'

    do j=1,mfit
      do k=1,mfit
        covar(j,k)=0.
      end do
      beta(j)=0.
    end do

    do i=1,ndat
      call funcs(x(i),afunc,ma)
      ym=y(i)
      if(mfit.lt.ma) then
        do j=1,ma
          if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
        end do
      endif
      sig2i=1./sig(i)**2
      j=0
      do l=1,ma
        if (ia(l).ne.0) then
          j=j+1
          wt=afunc(l)*sig2i
          k=0
          do m=1,l
            if (ia(m).ne.0) then
              k=k+1
              covar(j,k)=covar(j,k)+wt*afunc(m)
            end if
          end do
          beta(j)=beta(j)+ym*wt
        end if
      end do
    end do

    do j=2,mfit
      do k=1,j-1
        covar(k,j)=covar(j,k)
      end do
    end do

    call gaussj(covar,mfit,npc,beta,1,1)

    j=0
    do l=1,ma
      if(ia(l).ne.0) then
        j=j+1
        a(l)=beta(j)
      end if
    end do
    chisq=0.
    do i=1,ndat
      call funcs(x(i),afunc,ma)
      sum=0.
      do j=1,ma
        sum=sum+a(j)*afunc(j)
      end do
      chisq=chisq+((y(i)-sum)/sig(i))**2
    end do
    call covsrt(covar,npc,ma,ia,mfit)

    return
    end subroutine lfit

!-------------------------------------------------------------------------------
!
!   SUBROUTINE COVSRT
!
!-------------------------------------------------------------------------------

    subroutine covsrt(covar,npc,ma,ia,mfit)
     
    integer ma,mfit,npc,ia(ma)
    real covar(npc,npc)
    integer i,j,k
    real swap

    do i=mfit+1,ma
      do j=1,i
        covar(i,j)=0.
        covar(j,i)=0.
      end do
    end do
   
    k=mfit

    do j=ma,1,-1
      if (ia(j).ne.0) then
        do i=1,ma
          swap=covar(i,k)
          covar(i,k)=covar(i,j)
          covar(i,j)=swap
        end do
        do i=1,ma
          swap=covar(k,i)
          covar(k,i)=covar(j,i)
          covar(j,i)=swap
        end do
        k=k-1
      endif
    end do

    return
    end subroutine covsrt

!-------------------------------------------------------------------------------
!
!   SUBROUTINE GAUSSJ
!
!-------------------------------------------------------------------------------

    subroutine gaussj(a,n,np,b,m,mp)
       
    integer mp,n,np,NMAX
    real a(np,np),b(np,mp)
    parameter (NMAX=50)
    integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
    real big,dum,pivinv
          
    do j=1,n
      ipiv(j)=0
    end do
    do i=1,n
      big=0.
      do j=1,n
        if (ipiv(j).ne.1) then
          do k=1,n
            if (ipiv(k).eq.0) then
              if (abs(a(j,k)).ge.big) then
                big=abs(a(j,k))
                irow=j
                icol=k
              end if
            else if (ipiv(k).gt.1) then
              pause 'sigular matrix in gaussj'
            end if
          end do
        end if
      end do
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
        do l=1,n
          dum=a(irow,l)
          a(irow,l)=a(icol,l)
          a(icol,l)=dum
        end do
        do l=1,m
          dum=b(irow,l)
          b(irow,l)=b(icol,l)
          b(icol,l)=dum
        end do
      end if
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
        a(icol,l)=a(icol,l)*pivinv
      end do
      do l=1,m
        b(icol,l)=b(icol,l)*pivinv
      end do
      do ll=1,n
        if (ll.ne.icol) then
          dum=a(ll,icol)
          a(ll,icol)=0.
          do l=1,n
            a(ll,l)=a(ll,l)-a(icol,l)*dum
          end do
          do l=1,m
            b(ll,l)=b(ll,l)-b(icol,l)*dum
          end do
        end if  
      end do
    end do
    do l=n,1,-1
      if(indxr(l).ne.indxc(l)) then
        do k=1,n
          dum=a(k,indxr(l))
          a(k,indxr(l))=a(k,indxc(l))
          a(k,indxc(l))=dum
        end do
      end if
    end do

    return
    end subroutine gaussj

    subroutine funcs(x,afunc,ma)
      
    integer ma,MMAX
    parameter(MMAX=50)
    real x,product,afunc(MMAX)
      
    afunc(1)=1.
    do i=2,ma
      product=1
      do j=2,i
        product=product*x        
      end do
      afunc(i)=product
    end do

    return
    end subroutine funcs

    end module pwrspd
