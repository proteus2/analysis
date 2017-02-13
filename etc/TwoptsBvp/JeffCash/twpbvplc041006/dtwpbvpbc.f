c Driver for TWPBVPC with all the test problems
c
c  version of July 2006
c 
c  subroutine  used by TWPBVPL
c     fsub, dfsub, gsub, dgsub, initu
C 



      PROGRAM  DTWPBVPBC
      IMPLICIT NONE
      INTEGER IPROB, NMSH
      DOUBLE PRECISION ES
      LOGICAL BSUCC
C     
C
      WRITE(*,'(/,A,/)') 'PROGRAM DTWPBVPBC '
      
      WRITE(*,'(A)') 'Input the problem number (0 to 35): '
      WRITE(*,'(A)') '(0 = batch process - save to ''batchc.res'')' 
      READ(*,*) IPROB
      
      IF (IPROB .eq. 0) THEN
         CALL RUNBATCH()
         STOP
      ENDIF
      
      IF ((IPROB .gt. 35) .or. (IPROB .lt. 0)) THEN
         WRITE(*,*) ' The problem number is incorrect'
         STOP
      ENDIF   
            
      CALL RUNPROB(IPROB, 1.0D-8, 1.0D-8, BSUCC, NMSH, .FALSE., ES)
      STOP
      END

      SUBROUTINE RUNBATCH()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TOLS(3)
      LOGICAL BSUCC
      
      OPEN(UNIT=26,FILE='batchc.res',STATUS='UNKNOWN')
      
      TOLS(1) = 1.0D-4
      TOLS(2) = 1.0D-6
      TOLS(3) = 1.0D-8

      DO IPROB=1,35
         WRITE(26,*) 'NP'
         DO J=1,3
            EPS = 1.0D0
            DO WHILE (.TRUE.)
              CALL RUNPROB(IPROB, EPS, TOLS(J), BSUCC, NMSH, .TRUE., ES)
                IF (BSUCC) THEN
2609            FORMAT(1x, I2, A, I10, A, D10.3, A, D10.3, A, D10.3)
                WRITE(26,2609) IPROB, ',', NMSH, ',', EPS, ',',
     +                         TOLS(J), ',', ES
                ELSE
                EXIT
                END IF
                EPS = EPS/1.0D1
                IF (EPS .LT. 1.0D-6) EXIT
            END DO
         END DO
      END DO
      END
      
      SUBROUTINE RUNPROB(IPROB, EPSIN, ETOL, BSUCC, NMSH, BATCH, ES)
      IMPLICIT NONE
      INTEGER LWRKFL, LWRKIN, NUDIM, NXXDIM, LISERIES, IPRINT, IND,IWRK,
     +        IPROB, NTOL, IPAR, NCOMP, LTOL, NLBC, NMSH, NFXPNT, NMAX,
     +        IFLBVP, ISERIES, INDNMS, I, NMSHACC, K, J, NMINIT, IDUM
      
      DOUBLE PRECISION PI, EPS, EPSIN, TOL, ETOL, WRK, RPAR, ALEFT,
     +                 ARIGHT, UVAL0, FIXPNT, XX, U, CKAPPA1, GAMMA1,
     +                 CKAPPA, FIXPNTACC, XXACC, UACC, XPT, ES, TNORM,
     +                 ERRMAX, XPTOLD, EXX, EXACT, DABS
      PARAMETER(LWRKFL=2500000)
      PARAMETER(LWRKIN=500000)
      DIMENSION WRK(LWRKFL),IWRK(LWRKIN)
      PARAMETER(NUDIM=6,NXXDIM=100000)
      DIMENSION U(NUDIM,NXXDIM), XX(NXXDIM)
      DIMENSION UACC(NUDIM,NXXDIM), XXACC(NXXDIM)      
      DIMENSION TOL(NUDIM), LTOL(NUDIM), FIXPNT(1), FIXPNTACC(NXXDIM)
      DIMENSION RPAR(2), IPAR(1)
      PARAMETER (LISERIES=100)
      LOGICAL   LINEAR, GIVEU, GIVMSH, BSUCC, BATCH
      EXTERNAL  TWPBVPC, FSUB, DFSUB, GSUB, DGSUB
      
      LOGICAL PDEBUG, use_c, comp_c
c this common need to be defined in order to run TWPBVPL successfully
c give information about some parameters      
      COMMON /ALGPRS/ NMINIT,PDEBUG,IPRINT,IDUM,UVAL0,use_c,comp_c

      BSUCC = .FALSE.
      PI=4.0D0*DATAN(1.0D0)
    
c     If you do not like to use the conditioning in the mesh selection
c     USE_C  = .false.
c     If you do not like to compute the conditioning parameters
c     COMP_C = .false.     
      USE_C  =  .true.
      COMP_C =  .true.
      
      EPS = EPSIN
C     take the reciprocal if we are dealing with problem 34
      IF (IPROB .EQ. 34) EPS=1.D0/EPSIN
      TOL(1) = ETOL      
     
      IPRINT = -1

      IF (.NOT. BATCH) THEN      
        write(*,'(A)') 'Input epsilon '    
        read(*,*) eps
        write(*,'(A)') 'Input tolerance '   
        read(*,*) TOL(1)
      ENDIF

      DO ind=1,LWRKFL
       WRK(ind) = 0d0
      ENDDO
      DO ind=1,LWRKIN
       IWRK(ind) = 0
      ENDDO

 
      write(6,*) 'Problem =', iprob
      write(6,*) 'Tol = ',tol(1),' eps = ',eps



      IF (IPROB.EQ.33) THEN
        NTOL = 6
      ELSE
        NTOL = 2
      END IF 
      
      IPAR(1)=IPROB
      RPAR(1)=EPS
      RPAR(2)=PI
      
      CALL INITIAL(NCOMP,ALEFT,ARIGHT,NTOL,TOL,LTOL,RPAR,IPAR)

      IF (IPROB.LE.30.or.iprob.ge.34) THEN
        NLBC=1
      ELSE IF (IPROB.LE.32) THEN
        NLBC=2
      ELSE IF (IPROB.EQ.33) THEN
        NLBC=3
      ENDIF
      
      NMSH  = 0
      NFXPNT= 0
c The initial approximation is equal to UVAL0      
      UVAL0 = 0.D0     
      IF (IPROB.EQ.24) UVAL0=0.5D0
      IF (IPROB.EQ.25) UVAL0=0.5D0
      IF (IPROB.EQ.27) UVAL0=0.5D0
      IF (IPROB.EQ.28) UVAL0=0.5D0
      IF (IPROB.EQ.29) UVAL0=0.5D0
    

      IF (IPROB.LT.19.OR.IPROB.GT.34) THEN
         LINEAR = .TRUE.
      ELSE
         LINEAR = .FALSE.
      ENDIF

      GIVEU  = .FALSE.
      GIVMSH = .FALSE.
      PDEBUG = .FALSE.
    
C
C  main subroutine

      CALL TWPBVPC(NCOMP,NLBC,ALEFT,ARIGHT,NFXPNT,FIXPNT,NTOL,
     +            LTOL,TOL,LINEAR,GIVMSH,GIVEU,NMSH,NXXDIM,
     +            XX,NUDIM,U,NMAX,LWRKFL,WRK,LWRKIN,IWRK,
     +            FSUB,DFSUB,GSUB,DGSUB,
     +            ckappa1,gamma1,ckappa,rpar,ipar,IFLBVP)

C     SCMODIFIED     
C     CHECK THE NUMBER OF MESHPOINTS
      IF (NMSH .gt. NXXDIM) THEN
         WRITE(*,*) 'Assertion Failure, max meshpoints bypassed!'
         WRITE(*,*) 'Returning to avoid corruption...'
         RETURN
      ENDIF
      
C     NOW WE WANT TO MEASURE THE DIFFERENCE BETWEEN SOLVED MESH
C     AND A "MORE ACCURATE" MESH      
      DO I=2,NMSH-1
        FIXPNTACC(2*(I-1)-1)=(XX(I-1)+XX(I))/2.0D0
        FIXPNTACC(2*(I-1))=XX(I)
      END DO
      FIXPNTACC(2*(NMSH-2)+1)=(XX(NMSH-1)+XX(NMSH))/2.0
      NFXPNT= 2*(NMSH-2)+1
      NMSHACC= 0
      
      TOL(1) = EPS/1.0D2
C      TOL(1) = EPS
      
      DO ind=1,LWRKFL
       WRK(ind) = 0d0
      ENDDO
      DO ind=1,LWRKIN
       IWRK(ind) = 0
      ENDDO
      
      CALL TWPBVPC(NCOMP,NLBC,ALEFT,ARIGHT,NFXPNT,FIXPNTACC,NTOL,
     +            LTOL,TOL,LINEAR,GIVMSH,GIVEU,NMSHACC,NXXDIM,
     +            XXACC,NUDIM,UACC,NMAX,LWRKFL,WRK,LWRKIN,IWRK,
     +            FSUB,DFSUB,GSUB,DGSUB,
     +            ckappa1,gamma1,ckappa,rpar,ipar,IFLBVP)    
          
      IF (NMSHACC .gt. NXXDIM) THEN
         WRITE(*,*) 'Assertion Failure, max meshpoints bypassed!'
         WRITE(*,*) 'Returning to avoid corruption...'
         RETURN
      ENDIF     

      IF (NMSHACC .eq. NMSH) THEN
         WRITE(*,*) 'Assertion Failure, ACC mesh not good enough!'
         RETURN
      ENDIF
      
      IF (NMSHACC .eq. 0) THEN
         WRITE(*,*) 'Computation of ''accurate'' mesh failed.'
         RETURN
      ENDIF

C
C.... CHECK ERRORS 
C
      XPT=ALEFT
      ES = TNORM(U, UACC, 1, 1, NLBC*2, NUDIM)
      ES = ES + TNORM(U, UACC, NMSH, NMSHACC, NLBC*2, NUDIM)
      
C     LOOP OVER THE OTHER MESH POINTS      
      K = 2
      DO I=2,NMSH-1
        DO WHILE(DABS(XX(I)-XXACC(K)) .gt. 1.0D-10)
            K = K + 1
            
            IF (K .ge. NMSHACC) THEN
                WRITE(*, *)
     +          'Assertion failure, could not check all points!'
                RETURN
            END IF            
        END DO
        ES= ES +TNORM(U, UACC, I, K, NLBC*2, NUDIM)
      END DO
      
      ES = ES/NMSH
       
      ERRMAX=1.0D-16
      
      IF (.NOT. BATCH) THEN
       if (iprob.le.21.or.iprob.eq.35) then
        write(*,*) '    X             U(X)       EXACT      ERROR  ' 
      else
        write(*,*) '    X             U(X)         '
       endif 
      end if
      
      
      DO 36 j=1,NMSH
        XPTOLD=XPT
        XPT=XX(J)

        if (iprob.le.21.or.iprob.eq.35) then
         if (j.eq.1.or.j.eq.nmsh) then 
             EXX=u(1,j)
         else  
            EXX = EXACT(XPT,RPAR,IPAR)
         end if   
          ERRMAX=DMAX1(ERRMAX,DABS((U(1,J)-EXX)/(1d0+DABS(EXX))))
          IF (.NOT. BATCH) THEN
          write(6,109) xpt,u(1,j),EXX,DABS((U(1,J)-EXX))/(1d0+DABS(EXX))
         endif    
  109 FORMAT(1x,  4(D12.3))
        else 
  108 FORMAT(1x,  2(D12.3))   
          IF (.NOT. BATCH) THEN
          write(6,108) xpt,u(1,j)
         endif 
          ERRMAX=0.D0
        end if
      
       
        
  36  CONTINUE

      IF (.NOT. BATCH) THEN  
      write(*,'(/,/)')

      if (iprob.le.21.or.iprob.eq.35) then
         if (comp_c) then
        
           write(6,*) '   eps      nmsh    error      kappa1     ',
     *        ' gamma1      kappa   '
        
           write(6,103) eps,  NMSH, ERRMAX, ckappa1, gamma1, ckappa

         else

           write(6,*) '   eps     nmsh      error   '
           write(6,103) eps, NMSH, ERRMAX
         endif

      else
      
              if (comp_c) then
        
           write(6,*) '   eps     nmsh     kappa1     ',
     *        ' gamma1       kappa  '
        
           write(6,103) eps, NMSH,  ckappa1, gamma1, ckappa
 
         else

           write(6,*) '   eps     nmsh '
           write(6,103) eps, NMSH
         endif
      
      
      endif
      write(6,*)'------------------------------------------------'
      ENDIF
        
  100 FORMAT(/,/)
  103 FORMAT(1x, D9.3, I8, 5(D12.3))
       
 1555 continue
 2555 continue 
 
      IF (iflbvp .eq. 0) BSUCC = .TRUE.
      IF (.not. BATCH) WRITE(*,*) 'Actual accuracy = ', es

      RETURN
      END
      
      DOUBLE PRECISION FUNCTION TNORM(A, B, I1, I2, NCMP, nudim)
      IMPLICIT NONE
      INTEGER I, NCMP, I1, I2, nudim
      DOUBLE PRECISION A, B, DEN
      DIMENSION A(nudim,*), B(nudim,*)
      
      TNORM = 0.D0
      
      DO I=1, NCMP
       DEN = DABS(B(I,I2))
       DEN = MAX(DEN, 1.0D0)
       TNORM = TNORM + (A(I, I1)-B(I, I2))*(A(I, I1)-B(I, I2))/(DEN*DEN)
      END DO
         
      TNORM = DSQRT(TNORM)
      RETURN      
      END

      SUBROUTINE INITIAL(NCOMP,ALEFT,ARIGHT,NTOL,TOL,LTOL,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TOL(*),LTOL(*)
      DIMENSION RPAR(*), IPAR(*)
      INTEGER J
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
     
 
        DO 44 J=1,ntol
         LTOL(J)=J
         TOL(J) =TOL(1)
 44     CONTINUE

      IF (IPROB.LE.30.or.IPROB.GE.34) THEN
        NCOMP = 2
      ELSE IF (IPROB.LE.32) THEN
        NCOMP = 4
      ELSE 
        NCOMP = 6
      ENDIF
C
      
      IF (IPROB.EQ.1.OR.IPROB.EQ.2.OR.IPROB.EQ.8.OR.IPROB.EQ.16.OR.
     *    (IPROB.GT.18.AND.IPROB.LE.34)) THEN
        ALEFT  = 0.D0
        ARIGHT = 1.0D0
      ELSE IF (IPROB.EQ.17) THEN
        ALEFT = -0.1D0
        ARIGHT = 0.1D0
      ELSE IF (IPROB.EQ.18) THEN
        ALEFT = 0.D0
        ARIGHT = 0.25D0
      ELSE 
        ALEFT = -1.0D0
        ARIGHT = 1.0D0
      ENDIF
      RETURN
      END

C
      SUBROUTINE FSUB(NCOMP,X,Z,F,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),F(*)
      DIMENSION RPAR(*), IPAR(*)
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
      
      IF (IPROB.EQ.1) THEN
       F(1)=Z(2)
       F(2)=Z(1)/EPS
      ELSE IF (IPROB.EQ.2) THEN
       F(1)=Z(2)
       F(2)=Z(2)/EPS
      ELSE IF (IPROB.EQ.3) THEN
       F(1)=Z(2)
       F(2)=(-(2.0D0+DCOS(PI*X))*Z(2)+Z(1)-(1.0D0+EPS*PI*PI)*DCOS(PI*X)
     +                            -(2.0D0+DCOS(PI*X))*PI*SIN(PI*X))/EPS
      ELSE IF (IPROB.EQ.4) THEN
       F(1)=Z(2)
       F(2)=((1.0D0+EPS)*Z(1)-Z(2))/EPS
      ELSE IF (IPROB.EQ.5) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(Z(1)+X*Z(2)-(1.0D0+EPS*PI*PI)*DCOS(PIX)+X*PI*DSIN(PIX))/EPS
      ELSE IF (IPROB.EQ.6) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(-X*Z(2)-EPS*PI*PI*DCOS(PIX)-PI*X*DSIN(PIX))/EPS
      ELSE IF (IPROB.EQ.7) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(-X*Z(2)+Z(1)-(1.0D0+EPS*PI*PI)*DCOS(PIX)-PIX*DSIN(PIX))/EPS
      ELSE IF (IPROB.EQ.8) THEN
       F(1)=Z(2)
       F(2)=-Z(2)/EPS
      ELSE IF (IPROB.EQ.9) THEN
       F(1)=Z(2)
       F(2)=-(4.0D0*X*Z(2)+2.0D0*Z(1))/(EPS+X*X)
      ELSE IF (IPROB.EQ.10) THEN
       F(1)=Z(2)
       F(2)=-X*Z(2)/EPS
      ELSE IF (IPROB.EQ.11) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(Z(1)-EPS*PI*PI*DCOS(PIX)-DCOS(PIX))/(EPS)
      ELSE IF (IPROB.EQ.12) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(Z(1)-EPS*PI*PI*DCOS(PIX)-DCOS(PIX))/(EPS)
      ELSE IF (IPROB.EQ.13) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(Z(1)-EPS*PI*PI*DCOS(PIX)-DCOS(PIX))/(EPS)
      ELSE IF (IPROB.EQ.14) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=(Z(1)-EPS*PI*PI*DCOS(PIX)-DCOS(PIX))/(EPS)
      ELSE IF (IPROB.EQ.15) THEN
       F(1)=Z(2)
       F(2)=X*Z(1)/EPS
      ELSE IF (IPROB.EQ.16) THEN
       F(1)=Z(2)
       F(2)=-Z(1)*PI*PI/(4.0D0*EPS*EPS)
      ELSE IF (IPROB.EQ.17) THEN
       F(1)=Z(2)
       F(2)=-3.0D0*EPS*Z(1)/(EPS+X*X)**2
      ELSE IF (IPROB.EQ.18) THEN
       F(1)=Z(2)
       F(2)=-Z(2)/EPS
      ELSE IF (IPROB.EQ.19) THEN
       PIX=PI*X
       F(1)=Z(2)
       F(2)=((PI/2.D0)*DSIN(PIX/2.D0)*DEXP(2.D0*Z(1))
     *                                 -DEXP(Z(1))*Z(2))/EPS
      ELSE IF (IPROB.EQ.20) THEN
       F(1)=Z(2)
       F(2)=(1.D0-Z(2)*Z(2))/EPS
      ELSE IF (IPROB.EQ.21) THEN
       F(1)=Z(2)
       F(2)=(Z(1)*(1.D0+Z(1))-DEXP((-2.D0*X)/EPS))/(EPS*EPS)
      ELSE IF (IPROB.EQ.22) THEN
       F(1)=Z(2)
       F(2)=-(Z(2)+Z(1)*Z(1))/EPS
      ELSE IF (IPROB.EQ.23) THEN
       F(1)=Z(2)
       F(2)=DSINH(Z(1)/EPS)/EPS
      ELSE IF (IPROB.EQ.24) THEN
       AX=1.D0+X**2
       APX=2.D0*X
       GA=1.4D0
       F(1)=Z(2)
       F(2)=(((1.D0+GA)/2.D0-EPS*APX)*Z(1)*Z(2)-Z(2)/Z(1)-
     *          (APX/AX)*(1.D0-(GA-1.D0)*Z(1)**2/2.D0))/(EPS*AX*Z(1))
      ELSE IF (IPROB.GE.25.AND.IPROB.LE.30) THEN
       F(1)=Z(2)
       F(2)=(Z(1)*(1.D0-Z(2)))/EPS
      ELSE IF (IPROB.EQ.31) THEN
       F(1)=DSIN(Z(2))
       F(2)=Z(3)
       F(3)=-Z(4)/EPS
       F(4)=((Z(1)-1.D0)*DCOS(Z(2))-Z(3)*(1.D0/DCOS(Z(2))+
     *                                       EPS*Z(4)*DTAN(Z(2))))/EPS
      ELSE IF (IPROB.EQ.32) THEN
       F(1)=Z(2)
       F(2)=Z(3)
       F(3)=Z(4)
       F(4)=(Z(2)*Z(3)-Z(1)*Z(4))/EPS
      ELSE IF (IPROB.EQ.33) THEN
       F(1)=Z(2)
       F(2)=(Z(1)*Z(4) - Z(3)*Z(2))/EPS
       F(3)=Z(4)
       F(4)=Z(5)
       F(5)=Z(6)
       F(6)=(-Z(3)*Z(6) - Z(1)*Z(2))/EPS
      ELSE IF (IPROB.EQ.34) THEN
       F(1)=Z(2)
       F(2)=-EPS*DEXP(Z(1))
      ELSE IF (IPROB.EQ.35) THEN
       F(1)=Z(2)
       F(2)=X*Z(2)/EPS-Z(1)/EPS
      ENDIF
      RETURN
      END 

      SUBROUTINE DFSUB(NCOMP,X,Z,DF,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),DF(NCOMP,*)
      DIMENSION RPAR(*), IPAR(*)
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
    
      IF (IPROB.EQ.1) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.D0/EPS
       DF(2,2)=0.D0
      ELSE IF (IPROB.EQ.2) THEN
       DF(1,1)=0.D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=1.0D0/EPS
      ELSE IF (IPROB.EQ.3) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/EPS
       DF(2,2)=-(2.0D0+DCOS(PI*X))/EPS
      ELSE IF (IPROB.EQ.4) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=(1.0D0+EPS)/EPS
       DF(2,2)=-1.0D0/EPS
      ELSE IF (IPROB.EQ.5) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/EPS
       DF(2,2)=X/EPS
      ELSE IF (IPROB.EQ.6) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=-X/EPS
      ELSE IF (IPROB.EQ.7) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/EPS
       DF(2,2)=-X/EPS
      ELSE IF (IPROB.EQ.8) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=-1.0D0/EPS
      ELSE IF (IPROB.EQ.9) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=-2.0D0/(EPS+X*X)
       DF(2,2)=-4.0D0*X/(EPS+X*X)
      ELSE IF (IPROB.EQ.10) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=-X/EPS
      ELSE IF (IPROB.EQ.11) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/(EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.12) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/(EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.13) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/(EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.14) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=1.0D0/(EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.15) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=X/EPS
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.16) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=-PI*PI/(4.0D0*EPS*EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.17) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=-3.0D0*EPS/(EPS+X*X)**2
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.18) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=-1.D0/EPS
      ELSE IF (IPROB.EQ.19) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=(PI*DSIN(PI*X/2.D0)*DEXP(2.D0*Z(1))-DEXP(Z(1))*Z(2))/EPS
       DF(2,2)=-DEXP(Z(1))/EPS
      ELSE IF (IPROB.EQ.20) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=0.0D0
       DF(2,2)=-2.D0*Z(2)/EPS
      ELSE IF (IPROB.EQ.21) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=(1.D0+2.D0*Z(1))/(EPS*EPS)
       DF(2,2)=0.0D0
      ELSE IF (IPROB.EQ.22) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=-2.D0*Z(1)/EPS
       DF(2,2)=-1.D0/EPS
      ELSE IF (IPROB.EQ.23) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=DCOSH(Z(1)/EPS)/(EPS*EPS)
       DF(2,2)=0.D0
      ELSE IF (IPROB.EQ.24) THEN
       AX=1.D0+X**2
       APX=2.D0*X
       GA=1.4D0
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=((2.D0*Z(2))/Z(1)**3+APX/(AX*Z(1)**2)+
     *            (APX*(GA-1.D0))/(2.D0*AX))/(EPS*AX)
       DF(2,2)=(((1.D0+GA)/2.D0-EPS*APX)*Z(1)-1.D0/Z(1))/(EPS*AX*Z(1))
      ELSE IF (IPROB.GE.25.AND.IPROB.LE.30) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.0D0
       DF(2,1)=(1.D0-Z(2))/EPS
       DF(2,2)=-Z(1)/EPS
      ELSE IF (IPROB.EQ.31) THEN
       Z2CS=DCOS(Z(2))
       DF(1,1)=0.0D0
       DF(1,2)=Z2CS
       DF(1,3)=0.0D0
       DF(1,4)=0.0D0
       DF(2,1)=0.0D0
       DF(2,2)=0.0D0
       DF(2,3)=1.0D0
       DF(2,4)=0.0D0
       DF(3,1)=0.0D0
       DF(3,2)=0.0D0
       DF(3,3)=0.0D0
       DF(3,4)=-1.0D0/EPS
       DF(4,1)=Z2CS/EPS
       Z2SC=1.D0/Z2CS
       DF(4,2)=(-(Z(1)-1.D0)*DSIN(Z(2))-Z(3)*Z2SC*(DTAN(Z(2))+
     *                                       EPS*Z(4)*Z2SC))/EPS
       DF(4,3)=-(Z2SC+EPS*Z(4)*DTAN(Z(2)))/EPS
       DF(4,4)=(-Z(3)*EPS*DTAN(Z(2)))/EPS
      ELSE IF (IPROB.EQ.32) THEN
       DF(1,1)=0.0D0
       DF(1,2)=1.D0
       DF(1,3)=0.0D0
       DF(1,4)=0.0D0
       DF(2,1)=0.0D0
       DF(2,2)=0.0D0
       DF(2,3)=1.0D0
       DF(2,4)=0.0D0
       DF(3,1)=0.0D0
       DF(3,2)=0.0D0
       DF(3,3)=0.0D0
       DF(3,4)=1.0D0
       DF(4,1)=-Z(4)/EPS
       DF(4,2)=Z(3)/EPS
       DF(4,3)=Z(2)/EPS
       DF(4,4)=-Z(1)/EPS
      ELSE IF (IPROB.EQ.33) THEN
       DF(1,2)=1.0D0
       DF(2,1)=Z(4)/EPS
       DF(2,2)=-Z(3)/EPS
       DF(2,3)=-Z(2)/EPS
       DF(2,4)=Z(1)/EPS
       DF(3,4)=1.0D0
       DF(4,5)=1.0D0
       DF(5,6)=1.0D0
       DF(6,1)=-Z(2)/EPS
       DF(6,2)=-Z(1)/EPS
       DF(6,3)=-Z(6)/EPS
       DF(6,6)=-Z(3)/EPS
      ELSE IF(IPROB.EQ.34) THEN
       DF(1,1)=0
       DF(1,2)=1
       DF(2,1)=-EPS*DEXP(Z(1))
       DF(2,2)=0
      ELSE IF (IPROB.EQ.35) THEN
       DF(1,1) = 0.0D0
       DF(1,2) = 1.0D0
       DF(2,1) = -1.0D0/eps
       DF(2,2) = x/eps         
      ENDIF
      RETURN
      END 
C 
      SUBROUTINE GSUB(I,NCOMP,Z,G,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*)
      DIMENSION RPAR(*), IPAR(*)
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
      IF (IPROB.EQ.1) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)
      ELSE IF (IPROB.EQ.2) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)
      ELSE IF (IPROB.EQ.3) THEN
       G=Z(1)+1.0D0
      ELSE IF (IPROB.EQ.4) THEN
       IF (I.EQ.1) G=Z(1)-1.0D0-DEXP(-2.0D0)
       IF (I.EQ.2) G=Z(1)-1.0D0-DEXP(-2.0D0*(1.0D0+EPS)/EPS)
      ELSE IF (IPROB.EQ.5) THEN
       G=Z(1)+1.0D0
      ELSE IF (IPROB.EQ.6) THEN
       IF (I.EQ.1) G=Z(1)+2.0D0
       IF (I.EQ.2) G=Z(1)
      ELSE IF (IPROB.EQ.7) THEN
       IF (I.EQ.1) G=Z(1)+1.0D0
       IF (I.EQ.2) G=Z(1)-1.0D0
      ELSE IF (IPROB.EQ.8) THEN
       IF (I.EQ.1) G=Z(1)-1.0D0
       IF (I.EQ.2) G=Z(1)-2.0D0
      ELSE IF (IPROB.EQ.9) THEN
       G=Z(1)-1.0D0/(1.0D0+EPS)
      ELSE IF (IPROB.EQ.10) THEN
       IF (I.EQ.1) G=Z(1)-0.0D0
       IF (I.EQ.2) G=Z(1)-2.0D0
      ELSE IF (IPROB.EQ.11) THEN
       G=Z(1)+1.0D0
      ELSE IF (IPROB.EQ.12) THEN
       IF (I.EQ.1) G=Z(1)+1.0D0
       IF (I.EQ.2) G=Z(1)
      ELSE IF (IPROB.EQ.13) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)+1.0D0
      ELSE IF (IPROB.EQ.14) THEN
       G=Z(1)
      ELSE IF (IPROB.EQ.15) THEN
       G=Z(1)-1.0D0
      ELSE IF (IPROB.EQ.16) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)-DSIN(PI/(2.0D0*EPS))
      ELSE IF (IPROB.EQ.17) THEN
       IF (I.EQ.1) G=Z(1)+0.1D0/DSQRT(EPS+0.01D0)
       IF (I.EQ.2) G=Z(1)-0.1D0/DSQRT(EPS+0.01D0)
      ELSE IF (IPROB.EQ.18) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)-DEXP(-1.D0/(4.D0*EPS))
      ELSE IF (IPROB.EQ.19) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)
      ELSE IF (IPROB.EQ.20) THEN
        IF (I.EQ.1) THEN
           X = -0.745D0/EPS
           G=Z(1)-1.D0-EPS*(-X+DLOG((DEXP(2.D0*X)+1.D0)/2.D0))
        ELSE
           X = 0.255D0/EPS
           G=Z(1)-1.D0-EPS*(X+DLOG((DEXP(-2.D0*X)+1.D0)/2.D0))
        ENDIF
      ELSE IF (IPROB.EQ.21) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)-DEXP(-1.D0/EPS)
      ELSE IF (IPROB.EQ.22) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)-0.5D0
      ELSE IF (IPROB.EQ.23) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)-1.D0
      ELSE IF (IPROB.EQ.24) THEN
       IF (I.EQ.1) G=Z(1)-0.9129D0
       IF (I.EQ.2) G=Z(1)-0.375D0
      ELSE IF (IPROB.EQ.25) THEN
       IF (I.EQ.1) G=Z(1)+1.D0/3.D0
       IF (I.EQ.2) G=Z(1)-1.D0/3.D0
      ELSE IF (IPROB.EQ.26) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)+1.D0/3.D0
      ELSE IF (IPROB.EQ.27) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)-1.D0/3.D0
      ELSE IF (IPROB.EQ.28) THEN
       IF (I.EQ.1) G=Z(1)-1.D0
       IF (I.EQ.2) G=Z(1)-3.D0/2.D0
      ELSE IF (IPROB.EQ.29) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(1)-3.D0/2.D0
      ELSE IF (IPROB.EQ.30) THEN
       IF (I.EQ.1) G=Z(1)+7.D0/6.D0
       IF (I.EQ.2) G=Z(1)-3.D0/2.D0
      ELSE IF (IPROB.EQ.31) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(3)
       IF (I.EQ.3) G=Z(1)
       IF (I.EQ.4) G=Z(3)
      ELSE IF (IPROB.EQ.32) THEN
       IF (I.EQ.1) G=Z(1)
       IF (I.EQ.2) G=Z(2)
       IF (I.EQ.3) G=Z(1)-1.D0
       IF (I.EQ.4) G=Z(2)
      ELSE IF (IPROB.EQ.33) THEN
       IF (I.EQ.1) G=Z(1)+1.D0
       IF (I.EQ.2) G=Z(3)
       IF (I.EQ.3) G=z(4)
       IF (I.EQ.4) G=Z(1)-1.D0
       IF (I.EQ.5) G=Z(3)
       IF (I.EQ.6) G=Z(4)
      ELSE IF (IPROB.EQ.34) THEN
       G=Z(1)
      ELSE IF (IPROB.EQ.35) THEN
        if (i.eq.1) G = Z(1)-1.0d0
        if (i.eq.2) G = Z(1)-2.0d0  
      ENDIF
      RETURN 
      END 
C 
C      
      SUBROUTINE DGSUB(I,NCOMP,Z,DG,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),DG(*)
      DIMENSION RPAR(*), IPAR(*)
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
     
      DG(1)=1.D0
      DG(2)=0.D0
      IF (IPROB.EQ.31) THEN
        DG(1)=0.D0
        DG(2)=0.D0
        DG(3)=0.D0
        DG(4)=0.D0  
        IF (I.EQ.1.OR.I.EQ.3) DG(1)=1.D0
        IF (I.EQ.2.OR.I.EQ.4) DG(3)=1.D0
      ELSE IF (IPROB.EQ.32) THEN
        DG(1)=0.D0
        DG(2)=0.D0
        DG(3)=0.D0
        DG(4)=0.D0  
        IF (I.EQ.1.OR.I.EQ.3) DG(1)=1.D0
        IF (I.EQ.2.OR.I.EQ.4) DG(2)=1.D0
      ELSE IF (IPROB.EQ.33) THEN
        DG(1)=0.D0
        DG(2)=0.D0
        DG(3)=0.D0
        DG(4)=0.D0
        DG(5)=0.D0
        DG(6)=0.D0
        IF (I.EQ.1.OR.I.EQ.4) DG(1)=1.D0
        IF (I.EQ.2.OR.I.EQ.5) DG(3)=1.D0
        IF (I.EQ.3.OR.I.EQ.6) DG(4)=1.D0
      ENDIF
      RETURN
      END 
C
      
      DOUBLE PRECISION FUNCTION EXACT(X,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RPAR(*), IPAR(*)
      IPROB = IPAR(1)
      EPS   = RPAR(1)
      PI    = RPAR(2)
      
      IF (IPROB.EQ.1) THEN
       EXACT=(DEXP(-X/DSQRT(EPS))-DEXP(-(2.D0-X)/DSQRT(EPS)))
     +                             /(1.D0-DEXP(-2.D0/DSQRT(EPS)))
      ELSE IF (IPROB.EQ.2) THEN
       EXACT=(1.D0-DEXP((X-1.D0)/EPS))/(1.D0-DEXP(-1.D0/EPS))
      ELSE IF (IPROB.EQ.3) THEN
       EXACT=DCOS(PI*X)
      ELSE IF (IPROB.EQ.4) THEN
       AA=0.D0
       BB=(1.0D0+EPS)*(1.0D0+X)/EPS
       IF (BB.LT.100.0D0) AA=DEXP(-BB)
       EXACT=DEXP(X-1.0D0)+AA
      ELSE IF (IPROB.EQ.5) THEN
       EXACT=DCOS(PI*X)
      ELSE IF (IPROB.EQ.6) THEN
       PIX = PI*X
       SQEP = DSQRT(2.D0*EPS)
       EXACT = DCOS(PIX)+DERF(X/SQEP)/DERF(1.D0/SQEP)
      ELSE IF (IPROB.EQ.7) THEN
       CC=0.0D0
       DD=0.0D0
       AA=X*X/(2.0D0*EPS)
       BB=1.0D0/(2.0D0*EPS)
       IF (AA.LT.100.0D0) CC=DEXP(-AA)
       IF (BB.LT.100.0D0) DD=DEXP(-BB)
       EXACT=DCOS(PI*X)+X+(X*DERF(X/SQRT(2.0D0*EPS))
     +                                 +SQRT(2.0D0*EPS/PI)*CC)/
     +            (DERF(1.0D0/SQRT(2.0D0*EPS))+SQRT(2.0D0*EPS/PI)*DD)
      ELSE IF (IPROB.EQ.8) THEN
       AA=0.0D0
       CC=0.0D0
       BB=X/EPS
       DD=1.0D0/EPS
       IF (BB.LT.100.0D0) AA=DEXP(-BB)
       IF (DD.LT.100.0D0) CC=DEXP(-DD)
       EXACT=(2.0D0-CC-AA)/(1.0D0-CC)
      ELSE IF (IPROB.EQ.9) THEN
       EXACT=1.0D0/(X*X+EPS)
      ELSE IF (IPROB.EQ.10) THEN
       EXACT=1.0D0+DERF(X/SQRT(2.0D0*EPS))/DERF(1.0D0/SQRT(2.0D0*EPS))
      ELSE IF (IPROB.EQ.11) THEN
       EXACT=DCOS(PI*X)
      ELSE IF (IPROB.EQ.12) THEN
       EXACT=DCOS(PI*X)+DEXP(-(1.0D0-X)/DSQRT(EPS))
      ELSE IF (IPROB.EQ.13) THEN
       EXACT=DCOS(PI*X)+DEXP(-(1.0D0+X)/DSQRT(EPS))
      ELSE IF (IPROB.EQ.14) THEN
       EXACT=DCOS(PI*X)+DEXP(-(1.0D0+X)/DSQRT(EPS))+
     *                                  DEXP(-(1.0D0-X)/DSQRT(EPS))
      ELSE IF (IPROB.EQ.15) THEN
       AA=0.0D0
       BB=0.0D0
       IF ((1.0D0-X)/EPS.LT.100.0D0) BB=DEXP((X-1.0D0)/EPS)
       IF ((1.0D0+X)/EPS.LT.100.0D0) AA=DEXP(-(X+1.0D0)/EPS)
       ZA=1.0D0/(EPS**(1.0D0/3.0D0))*X
       Z1=1.0D0/(EPS**(1.0D0/3.0D0))
C       Z2=S17AGF(-Z1,IFAIL)
C       Z3=S17AHF(-Z1,IFAIL)
C       Z4=S17AGF(Z1,IFAIL)
C       Z5=S17AHF(Z1,IFAIL)
C       BC=(Z4-Z2)/(Z4*Z3-Z2*Z5)
C       AC=(1.0D0-BC*Z3)/Z2
C       EXACT=AC*S17AGF(ZA,IFAIL)+BC*S17AHF(ZA,IFAIL)
      ELSE IF (IPROB.EQ.16) THEN
       EXACT=DSIN(PI*X/(2.0D0*EPS))
      ELSE IF (IPROB.EQ.17) THEN
       EXACT=X/SQRT(EPS+X*X)
      ELSE IF (IPROB.EQ.18) THEN
       EXACT=DEXP(-X/EPS)
      ELSE IF (IPROB.EQ.19) THEN
       EXACT=-DLOG((1.D0+DCOS(PI*X/2.D0))*
     *                            (1.D0-0.5D0*DEXP(-X/(2.D0*EPS))))
      ELSE IF (IPROB.EQ.20) THEN
       XX = (X-0.745D0)/EPS
       IF (XX.GT.0.D0) THEN
         EXACT=1.D0+EPS*(XX+DLOG((1.D0+DEXP(-2.D0*XX))/2.D0))
       ELSE
         EXACT=1.D0+EPS*(-XX+DLOG((1.D0+DEXP(2.D0*XX))/2.D0))
       ENDIF
      ELSE IF (IPROB.EQ.21) THEN
       EXACT=DEXP(-X/EPS)
      ELSE IF (IPROB.EQ.35) THEN
       IF (X.EQ.-1.0D+0) THEN
          EXACT=1.0D+0
       ELSE IF (X.EQ.1.D+0) THEN
          EXACT=2.0D+0
       ELSE IF (X.LE.-0.9999D+0) THEN
          EXACT=-0.5D+0 + 3.0D+0/2.0D+0*EXP(-(X+1)/EPS)
       ELSE IF (X.GE.0.9999D+0) THEN
          EXACT=0.5D+0 + 3.0D+0/2.0D+0*EXP(-(1-X)/EPS)
       ELSE
          EXACT=X/2.0D+0
       ENDIF 
      ENDIF
      RETURN
      END
      
      
      
      subroutine initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension xx(*), u(nudim, *)

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      iprob = ipar(1)

*  This routine must be provided to reset u after re-meshing 
*  for linear problems or for nonlinear problems
*  when interpolation of the old solution is not used.
 
      if (iprint .ne. -1) write(6,99) uval0
   99 format(1h ,'initu, uval0',1pd15.5)
      
      if (Iprob.eq.33) then 
       call mtload33(ncomp, nmsh, xx, nudim, u)
      else 
*  This version sets all elements of u to the constant uval0.
       call mtload(ncomp, nmsh, uval0, nudim, u)
      end if
      return
      end

      subroutine mtload33( nrow, ncol, xx, nrowx, xmat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrowx, ncol), xx(*)

*  for problem 33 the initial condition is set to 
*  xmat(1,:) = 2 x  - 1
*  xmat(2,:) = 2
*  xmat(3:end,:) 0
* 

      if (nrow .le. 0 .or. ncol .le. 0)  return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
C     SCMODIFIED: from
C     xmat(i,j) = 2.D0 * xx(i-1)-1.D0
      if (i.eq.1) then
        xmat(i,j) = 2.D0 * xx(j)-1.D0
      else if (i.eq.2) then
        xmat(i,j) = 2.D0
      else
        xmat(i,j) = 0.D0
      end if
  100 continue
      return
      end

        

