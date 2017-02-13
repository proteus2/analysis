MODULE LMExtern
!
   USE Base,     ONLY: i4, r4, r8
   USE CommonData
!
   IMPLICIT NONE
!
   PRIVATE
!
   INTERFACE FCN
      MODULE PROCEDURE FCNJAC
      MODULE PROCEDURE FCNNOJAC
   END INTERFACE
!
   PUBLIC :: FCN
   PUBLIC :: FCN_ORG
!
CONTAINS
!!!
!!!
   SUBROUTINE FCNJAC(m, n, x, fvec, fjac, iflag)
!
      USE Base, ONLY: i4, r8
      USE CommonData, ONLY: NDAT, NUMP, NUMQ, NCOF, SMAX, DS,  &
                            nfev, njev, s, fs
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: x(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      REAL(r8),    INTENT(OUT)   :: fjac(:,:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
      REAL(r8), DIMENSION(NUMP) :: pc
      REAL(r8), DIMENSION(NUMQ) :: qc
      REAL(r8), DIMENSION(NDAT) :: pp, qp
      INTEGER(i4) :: i, j, jj
!!!
!
      DO i = 1, NDAT
         s(i) = -SMAX + (i-1)*DS
      END DO
!
      pc(1:NUMP) = x(1:NUMP)
      qc(1:NUMQ) = x(NUMP+1:NCOF)
!
      DO i = 1, NDAT
         pp(i) = 0._r8
         qp(i) = 0._r8
         DO jj = 1, NUMP
            pp(i) = pp(i) + pc(jj) * (s(i)**(jj-1))
         END DO
         DO jj = 1, NUMQ
            qp(i) = qp(i) + qc(jj) * (s(i)**(jj-1))
         END DO
      END DO
! 
      IF (iflag == 1) THEN
         nfev = nfev + 1
         DO i = 1, NDAT
            fvec(i) = fs(i) - pp(i) / qp(i)
         END DO
      END IF
!
      IF (iflag == 2) THEN
         njev = njev + 1
         DO i = 1, NDAT
            DO j = 1, NUMP
               jj = j
               fjac(i,j) = -(s(i)**(jj-1)) / qp(i)
            END DO
            DO j = NUMP+1, NCOF
               jj = j - NUMP
               fjac(i,j) =  (s(i)**(jj-1)) / qp(i) * (pp(i) / qp(i))
            END DO
         END DO
      END IF
!
      RETURN
   END SUBROUTINE FCNJAC
!
   SUBROUTINE FCNNOJAC(m, n, x, fvec, iflag)
!
      USE Base, ONLY: i4, r4, r8
!
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: x(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
!     Local variables
!
      INTEGER(i4)         :: i
!!!
!
      RETURN
   END SUBROUTINE FCNNOJAC
!!!
!!!
SUBROUTINE FCN_ORG(m, fvec, nprob)
!
   USE Base, ONLY: i4, r8
   USE CommonData, ONLY: NDAT, NUMP, NUMQ, NCOF, SMAX, DS,  &
                         nfev, njev, s
!
   IMPLICIT NONE
!
   INTEGER(i4), INTENT(IN)  :: m
   REAL(r8),    INTENT(OUT) :: fvec(:)
   INTEGER(i4), INTENT(IN)  :: nprob
!
   REAL(r8), PARAMETER :: PI = 3.141592653589793_r8
   INTEGER(i4) :: i
!!!
!
   DO i = 1, NDAT
      s(i) = -SMAX + (i-1)*DS
   END DO
!
   SELECT CASE (nprob)
      CASE ( 1)   !   1/(s + 0.5)
         DO i = 1, m
            IF (s(i) /= -0.5_r8) THEN
              fvec(i) = 1._r8 / (s(i) + 0.5_r8)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == -0.5_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 2)   !  (s+0.2)/((s+0.2)^2+1)
         DO i = 1, m
            IF (s(i) /= 0.2_r8) THEN
              fvec(i) = (s(i)+0.2_r8) / ((s(i)+0.2_r8)**2 + 1._r8)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.2_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 3)   !  1/s
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
              fvec(i) = 1._r8 / s(i)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 4)   !  1/(s*s)
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
              fvec(i) = 1._r8 / (s(i) * s(i))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 5)   !  1/(s*s*s)
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
              fvec(i) = 1._r8 / (s(i) * s(i) * s(i))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 6)   !  1/sqrt(abs(s))
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
              fvec(i) = 1._r8 / SQRT(ABS(s(i)))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE ( 7)   !  1/(s*s + 1)
         DO i = 1, m
            fvec(i) = 1._r8 / (s(i) * s(i) + 1._r8)
         END DO
      CASE ( 8)   !  1/sqrt(s*s + 1)
         DO i = 1, m
            fvec(i) = 1._r8 / SQRT(s(i) * s(i) + 1._r8)
         END DO
      CASE ( 9)   !  exp(-1/s)/s
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = EXP(- 1._r8 / s(i)) / s(i)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (10)   !  1/(s*sqrt(abs(s)))
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = 1._r8 / (s(i) * SQRT(ABS(s(i))))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (11)   !  log((s*s+1)/(s*s + 4))
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = LOG((s(i)*s(i)+1._r8)/(s(i)*s(i)+4._r8))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (12)   !  atan(1/s)
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = ATAN(1._r8/s(i))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (13)   !  exp(-s)/s
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = EXP(-s(i)) / s(i)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (14)   !  exp(-s)/(s(s+1))
         DO i = 1, m
            IF (s(i) /= 0.0_r8 .AND. s(i) /= -1.0_r8) THEN
               fvec(i) = EXP(-s(i)) / (s(i) * (s(i) + 1._r8))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
            IF (s(i) == -1.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (15)   !  1/(s^2+s+1)
         DO i = 1, m
            fvec(i) = 1._r8 / (s(i)**2 + s(i) + 1._r8)
         END DO
      CASE (16)   ! Log(s)/s
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = LOG(s(i)) / s(i)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (17)   ! 1/(s - 1)
         DO i = 1, m
            IF (s(i) /= 1.0_r8) THEN
               fvec(i) = 1._r8 / (s(i) - 1._r8)
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 1.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO
      CASE (18)   ! exp(-4*sqrt(abs(s)))
         DO i = 1, m
            fvec(i) = EXP(-4._r8 * SQRT(ABS(s(i))))
         END DO
      CASE (19)   ! 1/(s^3 - 8)
         DO i = 1, m
            IF (s(i)**3 /= 8.0_r8) THEN
               fvec(i) = 1._r8 / (s(i)**3 - 8._r8)
            END IF
         END DO
         DO i = 1, m
            IF (s(i)**3 == 8.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO 
      CASE (20)   ! 1/(s*(1+exp(s)))
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = 1._r8 / (s(i)*(1._r8 + EXP(s(i))))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
              fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO 
      CASE (21)   ! 1/(s*(s+1))*(1+exp(-pi*s))/(1-exp(-pi*s))
         DO i = 1, m
            IF (s(i) /= 0.0_r8 .AND. s(i) /= -1.0_r8) THEN
               fvec(i) = 1._r8 / (s(i)*(s(i)+1._r8)) *  &
                         (1._r8 + EXP(-PI*s(i))) / (1._r8 - EXP(-PI*s(i)))
            END IF
         END DO
         DO i = 1, m
            IF (s(i) == 0.0_r8 .OR. s(i) == -1.0_r8) THEN
               fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO 
      CASE (22)   ! 1/((s-1)*sqrt(abs(s)))
         DO i = 1, m
            IF (s(i) /= 0.0_r8 .AND. s(i) /= 1.0_r8) THEN
               fvec(i) = 1._r8/((s(i)-1._r8)*SQRT(ABS(s(i))))
            END IF
         END DO 
         DO i = 1, m
            IF (s(i) == 0.0_r8 .OR. s(i) == 1.0_r8) THEN
               fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO 
      CASE (23)   ! exp(-5*sqrt(abs(s)))/s
         DO i = 1, m
            IF (s(i) /= 0.0_r8) THEN
               fvec(i) = EXP(-5._r8*SQRT(ABS(s(i)))) / s(i)
            END IF
         END DO 
         DO i = 1, m
            IF (s(i) == 0.0_r8) THEN
               fvec(i) = (fvec(i-1) + fvec(i+1))*0.5_r8
            END IF
         END DO 
      CASE (24)   ! 24 (5s^4+10s^2+1) / (s^2+1)^5
         DO i = 1, m
            fvec(i) = 24._r8 * (5._r8 * s(i)**4 + 10._r8 * s(i)**2 + 1._r8) /  &
                               (s(i)**2 + 1._r8)**5
         END DO 
   END SELECT
!
   RETURN
END SUBROUTINE FCN_ORG
!
END
