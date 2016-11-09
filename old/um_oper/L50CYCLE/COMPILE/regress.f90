    module regress

!----------------------------------------------------------------------------------
!
!   PURPOSE:
!
!   To perform one-dimensional regression using polynomial function
!
!
!   AUTHOR:
!
!   In-Sun Song.
!   Laboratory for Mesoscale Dynamics.
!   Department of Atmospheric Sciences, Yonsei University, Seoul, Korea.
!
!
!   TESTED SYSTEMS:
!
!   Compaq workstation: Compaq Unix 4.0 or 5.1.
!
!
!   VERSION HISTORY:
!
!   Imported from Numerical Recipes in Fortran 77.
!
!----------------------------------------------------------------------------------

    public :: lfit  ,covsrt, gaussj, funcs

    contains

!----------------------------------------------------------------------------------
!
!   Subroutine LFIT
!
!----------------------------------------------------------------------------------

    SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq)

    INTEGER max,ia(ma),npc,ndat,MMAX
    REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
    PARAMETER (MMAX=50)
    INTEGER i,j,k,l,m,mfit
    REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)
    mfit=0

    DO j=1,ma
      if (ia(j).ne.0) mfit=mfit+1
    ENDDO

    if(mfit.eq.0) then
!      print *, 'lfit: no parameters to be fitted'
    else
!      print *,' mfit =',mfit
    endif

    DO j=1,mfit
     DO k=1,mfit
       covar(j,k)=0.
     ENDDO
     beta(j)=0.
    ENDDO

    DO i=1,ndat
     call funcs(x(i),afunc,ma)
!     write(6,'(f8.1,11e10.2)') x(i),(afunc(ll),ll=1,11)
     ym=y(i)
     if(mfit.lt.ma) then
       DO j=1,ma
        if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
       ENDDO
     endif
     sig2i=1./sig(i)**2
     j=0
     DO l=1,ma
      if (ia(l).ne.0) then
       j=j+1
       wt=afunc(l)*sig2i
       k=0
       DO m=1,l
        if (ia(m).ne.0) then
         k=k+1
         covar(j,k)=covar(j,k)+wt*afunc(m)
        endif
       ENDDO
       beta(j)=beta(j)+ym*wt
      endif
     ENDDO
    ENDDO

    DO j=2,mfit
     DO k=1,j-1
      covar(k,j)=covar(j,k)
     ENDDO
    ENDDO

    CALL gaussj(covar,mfit,npc,beta,1,1)
    j=0
    DO l=1,ma
     if(ia(l).ne.0) then
      j=j+1
      a(l)=beta(j)
     endif
    ENDDO
    chisq=0.
    DO i=1,ndat
     CALL funcs(x(i),afunc,ma)
     sum=0.
     DO j=1,ma
      sum=sum+a(j)*afunc(j)
     ENDDO
     chisq=chisq+((y(i)-sum)/sig(i))**2
    ENDDO
    CALL covsrt(covar,npc,ma,ia,mfit)
    return
    END subroutine lfit

!----------------------------------------------------------------------------------
!
!   Subroutine COVSRT
!
!----------------------------------------------------------------------------------

    SUBROUTINE covsrt(covar,npc,ma,ia,mfit)

    INTEGER ma,mfit,npc,ia(ma)
    REAL covar(npc,npc)
    INTEGER i,j,k
    REAL swap

    DO i=mfit+1,ma
     DO j=1,i
       covar(i,j)=0.
       covar(j,i)=0.
     ENDDO
    ENDDO

    k=mfit
    DO j=ma,1,-1
     if (ia(j).ne.0) then
      DO i=1,ma
       swap=covar(i,k)
       covar(i,k)=covar(i,j)
       covar(i,j)=swap
      ENDDO
      DO i=1,ma
       swap=covar(k,i)
       covar(k,i)=covar(j,i)
       covar(j,i)=swap
      ENDDO
      k=k-1
     endif
    ENDDO
    return
    END subroutine covsrt

!----------------------------------------------------------------------------------
!
!   Subroutine GAUSSJ
!
!----------------------------------------------------------------------------------

    SUBROUTINE gaussj(a,n,np,b,m,mp)

    INTEGER mp,n,np,NMAX
    REAL a(np,np),b(np,mp)
    PARAMETER (NMAX=50)
    INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
    REAL big,dum,pivinv

    DO j=1,n
     ipiv(j)=0
    ENDDO
    DO i=1,n
     big=0.
     DO j=1,n
      if (ipiv(j).ne.1) then
       DO k=1,n
        if (ipiv(k).eq.0) then
         if (abs(a(j,k)).ge.big) then
          big=abs(a(j,k))
          irow=j
          icol=k
         endif
        else if (ipiv(k).gt.1) then
         pause 'sigular matrix in gaussj'
        endif
       ENDDO
      endif
     ENDDO
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
      DO l=1,n
       dum=a(irow,l)
       a(irow,l)=a(icol,l)
       a(icol,l)=dum
      ENDDO
      DO l=1,m
       dum=b(irow,l)
       b(irow,l)=b(icol,l)
       b(icol,l)=dum
      ENDDO
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     DO l=1,n
      a(icol,l)=a(icol,l)*pivinv
     ENDDO
     DO l=1,m
      b(icol,l)=b(icol,l)*pivinv
     ENDDO
     DO ll=1,n
      if (ll.ne.icol) then
       dum=a(ll,icol)
       a(ll,icol)=0.
       DO l=1,n
        a(ll,l)=a(ll,l)-a(icol,l)*dum
       ENDDO
       DO l=1,m
        b(ll,l)=b(ll,l)-b(icol,l)*dum
       ENDDO
      endif
     ENDDO
    ENDDO
    DO l=n,1,-1
     if(indxr(l).ne.indxc(l)) then
      DO k=1,n
       dum=a(k,indxr(l))
       a(k,indxr(l))=a(k,indxc(l))
       a(k,indxc(l))=dum
      ENDDO
     endif
    ENDDO
    return
    END subroutine gaussj

!----------------------------------------------------------------------------------
!
!   Subroutine FUNCS
!
!----------------------------------------------------------------------------------

    SUBROUTINE funcs(x,afunc,ma)

    INTEGER ma,MMAX
    PARAMETER(MMAX=50)
    REAL x,product,afunc(MMAX)

    afunc(1)=1.
    DO i=2,ma
     product=1
     DO j=2,i
      product=product*x
     ENDDO
     afunc(i)=product
    ENDDO
!   print *,' in regress.f '
    return
    END subroutine funcs

    end module regress
