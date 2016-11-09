
 
C   IMSL ROUTINE NAME   - FFTSC                                         FFTS0010
C                                                                       FFTS0020
C-----------------------------------------------------------------------FFTS0030
C                                                                       FFTS0040
C   COMPUTER            - IBM/SINGLE                                    FFTS0050
C                                                                       FFTS0060
C   LATEST REVISION     - JUNE 1, 1980                                  FFTS0070
C                                                                       FFTS0080
C   PURPOSE             - COMPUTE THE SINE AND COSINE TRANSFORMS OF     FFTS0090
C                           A REAL VALUED SEQUENCE                      FFTS0100
C                                                                       FFTS0110
C   USAGE               - CALL FFTSC (A,N,ST,CT,IWK,WK,CWK)             FFTS0120
C                                                                       FFTS0130
C   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH N WHICH           FFTS0140
C                           CONTAINS THE DATA TO BE TRANSFORMED.        FFTS0150
C                N      - INPUT NUMBER OF DATA POINTS TO BE TRANSFORMED.FFTS0160
C                           N MUST BE A POSITIVE EVEN INTEGER.          FFTS0170
C                ST     - OUTPUT REAL VECTOR OF LENGTH N/2+1            FFTS0180
C                           CONTAINING THE COEFFICIENTS OF THE          FFTS0190
C                           SINE TRANSFORM.                             FFTS0200
C                CT     - OUTPUT REAL VECTOR OF LENGTH N/2+1            FFTS0210
C                           CONTAINING THE COEFFICIENTS OF THE          FFTS0220
C                           COSINE TRANSFORM.                           FFTS0230
C                IWK    - INTEGER WORK VECTOR.                          FFTS0240
C                           IF N IS A POWER OF 2, THEN IWK SHOULD BE OF FFTS0250
C                           LENGTH M WHERE N=2**M.                      FFTS0260
C                           OTHERWISE, IWK SHOULD BE OF LENGTH          FFTS0270
C                           6*(N/2)+150.                                FFTS0280
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTS0290
C                WK     - REAL WORK VECTOR OF LENGTH 6*(N/2)+150.       FFTS0300
C                           WK IS NOT USED IF N IS A POWER OF 2.        FFTS0310
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTS0320
C                CWK    - COMPLEX WORK VECTOR OF LENGTH N/2+1.          FFTS0330
C                                                                       FFTS0340
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         FFTS0350
C                       - SINGLE/H36,H48,H60                            FFTS0360
C                                                                       FFTS0370
C   REQD. IMSL ROUTINES - FFTCC,FFTRC,FFT2C                             FFTS0380
C                                                                       FFTS0390
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           FFTS0400
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      FFTS0410
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  FFTS0420
C                                                                       FFTS0430
C   REMARKS  1.  FFTSC COMPUTES THE SINE TRANSFORM, ST, ACCORDING       FFTS0440
C                TO THE FOLLOWING FORMULA;                              FFTS0450
C                                                                       FFTS0460
C                  ST(K+1) = 2.0 * SUM FROM J = 0 TO N-1 OF             FFTS0470
C                            A(J+1)*SIN(2.0*PI*J*K/N)                   FFTS0480
C                  FOR K=0,1,...,N/2 AND PI=3.1415...                   FFTS0490
C                                                                       FFTS0500
C                FFTSC COMPUTES THE COSINE TRANSFORM, CT, ACCORDING     FFTS0510
C                TO THE FOLLOWING FORMULA;                              FFTS0520
C                                                                       FFTS0530
C                  CT(K+1) = 2.0 * SUM FROM J = 0 TO N-1 OF             FFTS0540
C                            A(J+1)*COS(2.0*PI*J*K/N)                   FFTS0550
C                  FOR K=0,1,...,N/2 AND PI=3.1415...                   FFTS0560
C            2.  THE FOLLOWING RELATIONSHIP EXISTS BETWEEN THE DATA     FFTS0570
C                AND THE COEFFICIENTS OF THE SINE AND COSINE TRANSFORM  FFTS0580
C                                                                       FFTS0590
C                  A(J+1) = CT(1)/(2*N) + CT(N/2+1)/(2*N)*(-1)**J +     FFTS0600
C                           SUM FROM K = 1 TO N/2-1 OF                  FFTS0610
C                            (CT(K+1)/N*COS((2.0*PI*J*K)/N) +           FFTS0620
C                             ST(K+1)/N*SIN((2.0*PI*J*K)/N))            FFTS0630
C                  FOR J=0,1,...,N-1 AND PI=3.1415...                   FFTS0640
C                                                                       FFTS0650
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       FFTS0660
C                                                                       FFTS0670
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN FFTS0680
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    FFTS0690
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        FFTS0700
C                                                                       FFTS0710
C-----------------------------------------------------------------------FFTS0720
C                                                                       FFTS0730
      SUBROUTINE FFTSC  (A,N,ST,CT,IWK,WK,CWK)                          FFTS0740
C                                  SPECIFICATIONS FOR ARGUMENTS         FFTS0750
      INTEGER            N,IWK(1)                                       FFTS0760
      REAL               A(N),ST(1),CT(1),WK(1)                         FFTS0770
      COMPLEX            CWK(1)                                         FFTS0780
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   FFTS0790
      INTEGER            I,ND2M1                                        FFTS0800
C                                  FIRST EXECUTABLE STATEMENT           FFTS0810
      ND2M1 = N/2+1                                                     FFTS0820
C                                  CALL FFTRC TO COMPUTE THE FOURIER    FFTS0830
C                                    COEFFICIENTS                       FFTS0840
      CALL FFTRC (A,N,CWK,IWK,WK)                                       FFTS0850
C                                  THE COSINE TRANSFORM COEFFICIENT IS  FFTS0860
C                                    THE REAL PART OF THE FOURIER       FFTS0870
C                                    COEFFICIENT AND THE SINE TRANSFORM FFTS0880
C                                    COEFFICIENT IS THE IMAGINARY PART  FFTS0890
      DO 5 I=1,ND2M1                                                    FFTS0900
         CT(I) = 2.0*REAL(CWK(I))                                       FFTS0910
         ST(I) = 2.0*AIMAG(CWK(I))                                      FFTS0920
    5 CONTINUE                                                          FFTS0930
      RETURN                                                            FFTS0940
      END                                                               FFTS0950
 
R; T=0.03/0.28 23:30:07
       FFTS0940
      END                                                               FFTS0950
 
R; T=0.03/