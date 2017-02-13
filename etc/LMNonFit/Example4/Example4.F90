PROGRAM Example4
!
!-------------------------------------------------------------------------------
!
!  Tests for fitting rational functions to tablated results of 
!  inverse Laplace transform
!
!-------------------------------------------------------------------------------
!
   USE Base,       ONLY: i4, r4, r8
   USE LMNonFit,   ONLY: LMDER1, LMDIF1, ENORM
   USE CommonData, ONLY: NDAT, NUMP, NUMQ, NCOF, nfev, njev, s, fs
   USE LMExtern,   ONLY: FCN, FCN_ORG
!
   IMPLICIT NONE
!
   REAL(r8),    DIMENSION(NDAT)      :: fvec
   REAL(r8),    DIMENSION(NDAT,NCOF) :: fjac
   REAL(r8),    DIMENSION(NCOF)      :: coef
   INTEGER(i4), DIMENSION(NCOF)      :: iwa
   INTEGER(i4), DIMENSION(NDAT)      :: ma, na, nf, nj, np, nx
   REAL(r8),    DIMENSION(NCOF,24)   :: coefall
   REAL(r8)    :: fnorm1, fnorm2, tol
   INTEGER(i4) :: i, j, k, info, iflag
!
!  LOGICAL INPUT  UNIT IS ASSUMED TO BE NUMBER 5.
!  LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
!
   INTEGER(i4), PARAMETER :: NREAD = 5, NWRITE = 6
!!!
!
   tol = SQRT(EPSILON(1._r8))
!
   DO k = 1, 24
      nfev = 0
      njev = 0
      DO j = 1, NCOF
         coef(j) = 1._r8
      END DO
      CALL FCN_ORG(NDAT, fs, k)
      CALL LMDER1(NDAT, NCOF, coef, fvec, fjac, tol, info, iwa)
      iflag = 1
      CALL FCN(NDAT, NCOF, coef, fvec, fjac, iflag)
      fvec = fs - fvec
      fnorm1 = ENORM(NDAT, fs)
      fnorm2 = ENORM(NDAT, fvec)
      WRITE(NWRITE,10) k, NCOF, NDAT
      WRITE(NWRITE,20) fnorm1, fnorm2, nfev, njev, info, coef(1:ncof)
      coefall(1:ncof,k) = coef(1:ncof)
   END DO
   DO k = 1, 24
      WRITE(6,*) coefall(1:ncof,k)
   END DO
!
5  FORMAT(I5,1X,E15.7,1X,E15.7)
10 FORMAT(' PROBLEM: ', I5, ' DIMENSIONS: ', 2I5)
20 FORMAT(' INITIAL L2 NORM OF THE RESIDUALS', G15.7/  &
          ' FINAL L2 NORM OF THE RESIDUALS  ', G15.7/  &
          ' NUMBER OF FUNCTION EVALUATIONS  ', I10/  &
          ' NUMBER OF JACOBIAN EVALUATIONS  ', I10/  &
          ' EXIT PARAMETER', I10/  &
          ' FINAL APPROXIMATE SOLUTION'/ (17G15.7)/)
30 FORMAT (' SUMMARY OF ', i3, ' CALLS TO LMDER1'/)
40 FORMAT (' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'/)
50 FORMAT (3I5, 3I6, ' ', g15.7)
!
END PROGRAM Example4
