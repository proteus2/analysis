
    module fft

!----------------------------------------------------------------------------------
!
!   PURPOSE:
!
!   To take forward and backward Fourier tranform of one-dimensional or
!   two-dimensional data. 
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
!    1 AUG, 2002: First wrote this module.
!   15 Aug, 2002: Bugs in fft2df is found and fixed by Mr. Jung-Suk Koh.
!   16 Feb, 2004: Reviewed and added descriptions.
!
!----------------------------------------------------------------------------------

    public :: fft1df,fft1db
    public :: fft2df,fft2db
    public :: wavnum
    public :: rearr_coef1d 
    public :: rearr_coef2d

    contains

!----------------------------------------------------------------------------------
!
!   SUBROUTINE fft1df
!
!   PURPOSE:
!
!   To calculate one-dimensional forward Fourier transform.
!   (real space -> spectral space)
!
!   INPUT:
!  
!   nx      : The number of one-dimensional grid points.
!   data    : Data to be taken forward Fourier transformed (real).
!             Data should be grided at a constant grid spacing (regular grid).
!
!   OUTPUT:
!
!   coef    : One-dimensional Fourier coefficients (double complex).
!
!----------------------------------------------------------------------------------

    subroutine fft1df(nx    ,data  ,coef  )

    implicit none

    integer                            ,intent(in)  :: nx
    real,            dimension(nx)     ,intent(in)  :: data
    double complex  ,dimension(nx)     ,intent(out) :: coef

    integer                              :: i
    double precision, dimension(4*nx+15) :: work

    work(1:4*nx+15) = 0.d0
    coef(1:nx) = (0.d0,0.d0)
 
    do i=1,nx
      coef(i) = dble(data(i))
    end do

    call CFFTI(nx,work)
    call CFFTF(nx,coef,work)

    return 
    end subroutine fft1df

!----------------------------------------------------------------------------------
!
!   SUBROUTINE fft1db
!
!   PURPOSE:
!
!   To calculate one-dimensional backward Fourier transform 
!   (spectral space -> real space)
!
!   INPUT:
!  
!   nx      : The number of one-dimensional grid points.
!   coef    : One-dimensional Fourier coefficients (double complex)
!
!   OUTPUT:
!
!   data    : Reconstructed data in real space (real).
!
!----------------------------------------------------------------------------------

    subroutine fft1db(nx    ,coef  ,data  )

    implicit none

    integer                        ,intent(in)  :: nx
    double complex  ,dimension(nx) ,intent(in)  :: coef
    real,            dimension(nx) ,intent(out) :: data

    integer                              :: i
    double precision, dimension(4*nx+15) :: work

    work(1:4*nx+15) = 0.d0
    data(1:nx) = 0.0

    call CFFTI(nx,work)
    call CFFTB(nx,coef,work)

    do i=1,nx
      data(i) = (1./float(nx))*real(coef(i))  ! real part only
    end do

    return
    end subroutine fft1db

!----------------------------------------------------------------------------------
!
!   SUBROUTINE fft2df
!
!   PURPOSE:
!
!   To calculate two-dimensional forward Fourier transform
!   (real space -> spectral space)
!
!   INPUT:
! 
!   nx      : The number of grid points in x direction.
!   ny      : The number of grid points in y direction.
!   data    : Two-dimensional data in real space (real).
!             Data should be grided at a constant grid spacing (regular grid).
!
!   OUTPUT:
!
!   coef    : Two-dimensional Fourier coefficients (double complex).
!
!----------------------------------------------------------------------------------

    subroutine fft2df(nx    ,ny    ,data  ,coef)

    implicit none

    integer                            ,intent(in)  :: nx     ,ny
    real,            dimension(nx,ny)  ,intent(in)  :: data
    double complex  ,dimension(nx,ny)  ,intent(out) :: coef

    integer                              :: i      ,j
    double precision, dimension(4*nx+15) :: workx
    double precision, dimension(4*ny+15) :: worky
    double complex  , dimension(nx)      :: coefx
    double complex  , dimension(ny)      :: coefy

    do j=1,ny
      do i=1,nx
        coef(i,j) = dble(data(i,j))
      end do
    end do

    do j=1,ny
      workx(1:4*nx+15) = 0.d0
      coefx(1:nx) = coef(1:nx,j)
      call cffti(nx,workx)
      call cfftf(nx,coefx,workx)
      coef(1:nx,j) = coefx(1:nx)
    end do

    do i=1,nx
      worky(1:4*ny+15) = 0.d0
      coefy(1:ny) = coef(i,1:ny)
      call cffti(ny,worky)
      call cfftf(ny,coefy,worky)
      coef(i,1:ny) = coefy(1:ny)
    end do 
 
    return 
    end subroutine fft2df

!----------------------------------------------------------------------------------
!
!   SUBROUTINE fft2db
!
!   PURPOSE:
!
!   To calculate two-dimensional backward Fourier transform.
!   (spectral space -> real space)
!
!   INPUT:
! 
!   nx      : The number of grid points in x direction.
!   ny      : The number of grid points in y direction.
!   coef    : Two-dimensional Fourier coefficients (double complex).
!
!   OUTPUT:
!
!   data    : Two-dimensional data in real space (real).
!
!----------------------------------------------------------------------------------

    subroutine fft2db(nx    ,ny    ,coef  ,data  )

    implicit none

    integer                            ,intent(in)  :: nx     ,ny
    double complex  ,dimension(nx,ny)  ,intent(in)  :: coef
    real            ,dimension(nx,ny)  ,intent(out) :: data

    integer                              :: i      ,j
    double precision, dimension(4*nx+15) :: workx
    double precision, dimension(4*ny+15) :: worky
    double complex  , dimension(nx)      :: coefx
    double complex  , dimension(ny)      :: coefy
    double complex  , dimension(nx,ny)   :: coef2

    coef2 = coef

    do j=1,ny
      workx(1:4*nx+15) = 0.d0
      coefx(1:nx) = coef2(1:nx,j)
      call cffti(nx,workx)
      call cfftb(nx,coefx,workx)
      coef2(1:nx,j) = coefx(1:nx)
    end do

    do i=1,nx
      worky(1:4*ny+15) = 0.d0
      coefy(1:ny) = coef2(i,1:ny)
      call cffti(ny,worky)
      call cfftb(ny,coefy,worky)
      coef2(i,1:ny) = coefy(1:ny)
    end do

    do i=1,nx
      do j=1,ny
        data(i,j) = 1./(float(nx)*float(ny))*real(coef2(i,j))
      end do
    end do 

    return
    end subroutine fft2db

!----------------------------------------------------------------------------------
!
!   SUBROUTINE wavnum
!
!   PURPOSE:
!
!   To calculate wavenumber or frequencies where the Fourier coefficients are
!   defined irrespective of whether the number of grids is even or odd.
!
!   INPUT:
!
!   nx      : The number of grid points in x direction.
!   dx      : Grid spacing
!   
!   OUTPUT:
!
!   freq    : Frequency or wavenumber
!   
!   DESCRIPTION:
!
!   When nx = 100                         When nx = 101                 
!  
!   freq(1)   =  0                        freq(1)   =  0  
!   freq(2)   =  1/(100*dx)               freq(2)   =  1/(101*dx)
!   freq(3)   =  2/(100*dx)               freq(3)   =  2/(101*dx)
!   freq(50)  = 49/(100*dx)               freq(50)  = 49/(101*dx)
!   freq(51)  = +/- 50/(100*dx)           freq(51)  = 50/(101*dx)
!   freq(52)  = -freq(50) = -49/(100*dx)  freq(52)  = -freq(51) = -50/(101*dx)
!   freq(53)  = -freq(49) = -48/(100*dx)  freq(53)  = -freq(50) = -49/(101*dx)
!   freq(99)  = -freq(3)  =  -2/(100*dx)  freq(99)  = -freq(4)  =  -3/(101*dx)
!   freq(100) = -freq(2)  =  -1/(100*dx)  freq(100) = -freq(3)  =  -2/(101*dx)
!                                         freq(101) = -freq(2)  =  -1/(101*dx)
!
!----------------------------------------------------------------------------------

    subroutine wavnum(nx    ,dx    ,freq  )

    implicit none

    integer,             intent(in)  :: nx
    real,                intent(in)  :: dx
    real, dimension(nx), intent(out) :: freq

    integer                          :: i

    do i=1,nx/2+1
      freq(i) = float(i-1)/(nx*dx)
    end do

    do i=nx/2+2,nx
      freq(i) = -freq(nx-i+2)
    end do

    return
    end subroutine wavnum    

!----------------------------------------------------------------------------------
!
!   SUBROUTINE rearr_coef1d
!
!   PURPOSE:
!
!   To rearrange one-dimensional Fourier coefficients
!
!   INPUT:
!
!   nx      : The number of grid points in x direction.
!   vari    : One-dimensional real or imaginary part of Fourier coefficients 
!             obtained using fft1df.
!   mx      : The number of grid points of rearranged Fourier coefficients.
!             mx should be equal to nx+1 when nx is an even number
!             mx should be equal to nx when nx is an odd number
!
!   OUTPUT:
!
!   varo    : Rearranged one-dimensional Fourier coefficients.
!
!   DESCRIPTION:  
!  
!   When nx = 100  (even number)
!   kx   :   1 -  51  ==>  51 - 101
!   kx   :  51 - 100  ==>   1 -  50
!
!   When nx = 101  (odd number)
!   kx   :   1 -  51  ==>  51 - 101
!   kx   :  52 - 101  ==>   1 -  50
!
!----------------------------------------------------------------------------------

    subroutine rearr_coef1d (nx    ,vari  ,mx    ,varo  )

    implicit none

    integer,             intent(in)  :: nx    ,mx
    real, dimension(nx), intent(in)  :: vari
    real, dimension(mx), intent(out) :: varo

    integer                          :: i     ,itmp1 ,itmp2
    logical                          :: flag

    if ( mod(nx,2) == 0 ) then
      flag = .true.      ! nx is an even number.
      if ( mx /= nx+1 ) then
        write (6,*) '(REARR_COEF1D): mx should be equal to nx+1 when nx is even number.'
        stop
      end if
    end if

    if ( mod(nx,2) /= 0 ) then
      flag = .false.     ! nx is an odd number.
      if ( mx /= nx ) then
        write (6,*) '(REARR_COEF1D): mx should be equal to nx when nx is odd number.'
        stop
      end if
    end if

    if ( flag ) then
      itmp1 = 1
      itmp2 = 0
    else
      itmp1 = 0
      itmp2 = 1
    end if

    do i=nx/2+1,nx+itmp1
      varo(i) = vari(i-nx/2)
    end do
    do i=1,nx/2
      varo(i) = vari(i+nx/2+itmp2)
    end do

    return
    end subroutine rearr_coef1d

!----------------------------------------------------------------------------------
!
!   SUBROUTINE rearr_coef2d
!
!   Purpose:
!
!   To rearrange 2-dimensional Fourier coefficients. 
!
!   INPUT:
!
!   nx      : The number of grid points in x direction.
!   ny      : The number of grid points in y direction.
!   vari    : Two-dimensional real or imaginary part of Fourier coefficients 
!             obtained using fft2df.
!   mx      : The number of grid points of rearranged coefficients in x direction.
!             mx should be equal to nx+1 when nx is an even number
!             mx should be equal to nx when nx is an odd number
!   my      : The number of grid points of rearranged coefficients in y direction.
!             my should be equal to ny+1 when ny is an even number
!             my should be equal to ny when ny is an odd number
!
!   OUTPUT:
!
!   varo    : Rearranged two-dimensional Fourier coefficients.
!
!
!   #############        #############
!   #     #     #        #     #     #
!   #  2  #  1  #        #  4  #  3  #
!   #     #     #        #     #     #
!   #############  ===>  #############
!   #     #     #        #     #     #
!   #  3  #  4  #        #  1  #  2  #
!   #     #     #        #     #     #
!   #############        #############
!
!----------------------------------------------------------------------------------

    subroutine rearr_coef2d (nx    ,ny    ,vari  ,mx    ,my    ,varo)

    implicit none

    integer                   , intent(in)    :: nx     ,ny
    integer                   , intent(in)    :: mx     ,my
    real, dimension(nx,ny)    , intent(in)    :: vari
    real, dimension(nx,ny)    , intent(out)   :: varo

    integer                                   :: i      ,j
    integer                                   :: itmp1  ,itmp2
    integer                                   :: jtmp1  ,jtmp2
    logical                                   :: flagx  ,flagy

    if ( mod(nx,2) == 0 ) then
      flagx = .true.      ! nx is an even number.
      if ( mx /= nx+1 ) then
        write (6,*) '(REARR_COEF2D): mx should be equal to nx+1 when nx is even number.'
        stop
      end if
    end if

    if ( mod(nx,2) /= 0 ) then
      flagy = .false.     ! nx is an odd number.
      if ( mx /= nx ) then
        write (6,*) '(REARR_COEF2D): mx should be equal to nx when nx is odd number.'
        stop
      end if
    end if

    if ( mod(ny,2) == 0 ) then
      flagy = .true.      ! ny is an even number.
      if ( my /= ny+1 ) then
        write (6,*) '(REARR_COEF2D): my should be equal to ny+1 when ny is even number.'
        stop
      end if
    end if

    if ( mod(ny,2) /= 0 ) then
      flagy = .false.     ! ny is an odd number.
      if ( my /= ny ) then
        write (6,*) '(REARR_COEF2D): my should be equal to ny when ny is odd number.'
        stop
      end if
    end if

    if ( flagx ) then
      itmp1 = 1
      itmp2 = 0
    else
      itmp1 = 0
      itmp2 = 1
    end if

    if ( flagy ) then
      jtmp1 = 1
      jtmp2 = 0
    else
      jtmp1 = 0
      jtmp2 = 1
    end if

!   The third quadrant -> the first quadrant

    do i=nx/2+1,nx+itmp1
    do j=ny/2+1,ny+jtmp1
      varo(i,j) = vari(i-nx/2,j-ny/2)
    end do
    end do

!   The fourth quadrant -> the second quadrant

    do i=1,nx/2
    do j=ny/2+1,ny+jtmp1
      varo(i,j) = vari(i+nx/2+itmp2,j-ny/2)
    end do
    end do

!   The second quadrant -> the fourth quadrant

    do i=nx/2+1,nx+itmp1
    do j=1,ny/2
      varo(i,j) = vari(i-nx/2,j+ny/2+jtmp2)
    end do
    end do

!   The first quadrant -> the third quadrant

    do i=1,nx/2
    do j=1,ny/2
      varo(i,j) = vari(i+nx/2+itmp2,j+ny/2+jtmp2)
    end do
    end do 

    return
    end subroutine rearr_coef2d

    end module fft

!**************************************************************************
!      FFTPACK: Available at http://www.scd.ucar.edu/softlib/mathlib.html
!
!   ***NOTE (C.T.) *** this is a subset of FFTPACK, which only includes
!            the routines for the complex FFT (forward & inverse).
!
! Modified: November 1999 by Arjan van Dijk to include IMPLICIT NONE and
!           to convert all routines to DOUBLE precision.
!**************************************************************************




!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!      *                                                               *
!      *                           FFTPACK                             *
!      *                                                               *
!      *     A package of Fortran subprograms for calculating          *
!      *     fast Fourier transforms for both complex and real         *
!      *      periodic sequences and certain other symmetric           *
!      *             sequences that are listed below                   *
!      *               (Version 4.1 November 1988)                     *
!      *                             by                                *
!      *                      Paul Swarztrauber                        *
!      *                             of                                *
!      *         The National Center for Atmospheric Research          *
!      *                Boulder, Colorado  (80307)  U.S.A.             *
!      *                   which is sponsored by                       *
!      *              the National Science Foundation                  *
!      *                                                               *
!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!  Any source code available through this distribution interface at NCAR
!  is free of charge, but there is no guarantee.
!
! FFTPACK breaks the FORTRAN 77 ANSI Standard
! by passing REAL arrays to subroutines and using the arrays within
! the subroutines as DOUBLE PRECISION or other types.  This infraction
! may cause data alignment problems when the source code is compiled
! and loaded in an executable.




!****************************************************************************
!     SUBROUTINE CFFTI(N,WSAVE)
!
!     SUBROUTINE CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.
!
      SUBROUTINE CFFTI (N,WSAVE)
      IMPLICIT NONE
      INTEGER N,IW1,IW2
      DOUBLE PRECISION WSAVE
      DIMENSION       WSAVE(*)
!
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

!****************************************************************************
!     SUBROUTINE CFFTF(N,C,WSAVE)
!
!     SUBROUTINE CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER
!     TRANSFORM (THE FOURIER ANALYSIS). EQUIVALENTLY , CFFTF COMPUTES
!     THE FOURIER COEFFICIENTS OF A COMPLEX PERIODIC SEQUENCE.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
!
!     THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED TRANSFORM
!     THE OUTPUT MUST BE DIVIDED BY N. OTHERWISE A CALL OF CFFTF
!     FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE SEQUENCE BY N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTF MUST BE
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES. N
!
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
!             IN THE PROGRAM THAT CALLS CFFTF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
!
!     OUTPUT PARAMETERS
!
!     C      FOR J=1,...,N
!
!                C(J)=THE SUM FROM K=1,...,N OF
!
!                      C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
!
!                            WHERE I=SQRT(-1)
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
!
      SUBROUTINE CFFTF (N,C,WSAVE)
      IMPLICIT NONE
      INTEGER N,IW1,IW2
      DOUBLE PRECISION C,WSAVE
      DIMENSION       C(*)       ,WSAVE(*)
!
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

!****************************************************************************
!     SUBROUTINE CFFTB(N,C,WSAVE)
!
!     SUBROUTINE CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
!     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , CFFTB COMPUTES
!     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
!
!     A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE
!     SEQUENCE BY N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTB MUST BE
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
!
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
!             IN THE PROGRAM THAT CALLS CFFTB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
!
!     OUTPUT PARAMETERS
!
!     C      FOR J=1,...,N
!
!                C(J)=THE SUM FROM K=1,...,N OF
!
!                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
!
!                            WHERE I=SQRT(-1)
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
!
      SUBROUTINE CFFTB (N,C,WSAVE)
      IMPLICIT NONE
      INTEGER N,IW1,IW2
      DOUBLE PRECISION C,WSAVE
      DIMENSION       C(*)       ,WSAVE(*)
!
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

!****************************************************************************
      SUBROUTINE CFFTI1 (N,WA,WIFAC)
      IMPLICIT NONE
      INTEGER N,NTRYH,NL,NF,J,NTRY,NQ,NR,IB,I,L1,K1,IP,LD,L2,  &
        IDO,IDOT,IPM,I1,II
      DOUBLE PRECISION WA,TPI,ARGH,FI,ARGLD,ARG,WIFAC
      DIMENSION       WA(*)      ,WIFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      WIFAC(NF+2) = DBLE(NTRY)
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         WIFAC(IB+2) = WIFAC(IB+1)
  106 CONTINUE
      WIFAC(3) = 2.D0
  107 IF (NL .NE. 1) GO TO 104
      WIFAC(1) = DBLE(N)
      WIFAC(2) = DBLE(NF)
      TPI = 2.D0*(4.D0*ATAN(1.D0))
      ARGH = TPI/DBLE(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = NINT(WIFAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.D0
            WA(I) = 0.D0
            LD = LD+L1
            FI = 0.D0
            ARGLD = DBLE(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE CFFTF1 (N,C,CH,WA,WIFAC)
      IMPLICIT NONE
      INTEGER N,NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,N2,  &
        I,NAC
      DOUBLE PRECISION C,CH,WA,WIFAC
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WIFAC(*)
      NF = NINT(WIFAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = NINT(WIFAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE CFFTB1 (N,C,CH,WA,WIFAC)
      IMPLICIT NONE
      INTEGER N,NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,N2,  &
        I,NAC
      DOUBLE PRECISION C,CH,WA,WIFAC
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WIFAC(*)
      NF = NINT(WIFAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = NINT(WIFAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT NONE
      INTEGER NAC,IDO,IP,L1,IDL1,IDOT,IPP2,IPPH,IDP,J,K,I,IDL,LC,L,  &
        IK,IDLJ,JC,INC,IDIJ,IDJ
      DOUBLE PRECISION CC,C1,C2,CH,CH2,WA,WAR,WAI
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,  &
                      C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),  &
                      CH2(IDL1,IP)
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,TR2,TI2
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,  &
                      WA1(1)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,  &
        DR2,DR3,DI2,DI3
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,   &
                      WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5D0,.866025403784439D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,WA3,TI1,TI2,TI3,TI4,TR1,TR2,TR3,  &
        TR4,CR3,CR2,CI3,CR4,CI4,CI2
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
                      WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,WA3,WA4,TR11,TI11,TR12,TI12,         &
        TI2,TI3,TI4,TI5,TR2,TR3,TR4,TR5,CR2,CI2,CR3,CI3,CR4,CI4,CR5,CI5,  &
        DR3,DI3,DR4,DI4,DR5,DI5,DR2,DI2
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,    &
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,.951056516295154D0,    &
      -.809016994374947D0,.587785252292473D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT NONE
      INTEGER NAC,IDO,IP,L1,IDL1,IDOT,IPP2,IPPH,IDP,J,JC,K,I,IDL,    &
        INC,L,LC,IK,IDLJ,IDIJ,IDJ
      DOUBLE PRECISION CC,C1,C2,CH,CH2,WA,WAR,WAI
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,   &
                      C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),   &
                      CH2(IDL1,IP)
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,TR2,TI2
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,  &
                      WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,  &
        DR2,DR3,DI2,DI3
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,   &
                      WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5D0,-.866025403784439D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,WA3,TI1,TI2,TI3,TI4,TR1,TR2,TR3,    &
        TR4,CR2,CR3,CR4,CI2,CI3,CI4
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,   &
                      WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

!****************************************************************************
      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT NONE
      INTEGER IDO,L1,K,I
      DOUBLE PRECISION CC,CH,WA1,WA2,WA3,WA4,TR11,TI11,TR12,TI12,        &
        TI2,TI3,TI4,TI5,TR2,TR3,TR4,TR5,CI2,CI3,CI4,CI5,CR2,CR3,CR4,     &
        CR5,DI2,DI3,DI4,DI5,DR2,DR3,DR4,DR5
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,   &
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,-.951056516295154D0,  &
      -.809016994374947D0,-.587785252292473D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
