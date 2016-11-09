C   IMSL ROUTINE NAME   - FTFPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - FAST FOURIER TRANSFORM ESTIMATES OF POWER
C                           SPECTRA AND CROSS SPECTRA OF TIME SERIES
C
C   USAGE               - CALL FTFPS (X,Y,N,L,IND,PSX,PSY,XPS,IWK,WK,
C                           CWK,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           FIRST TIME SERIES. SEE REMARKS.
C                           X IS DESTROYED ON OUTPUT.
C                Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SECOND TIME SERIES. DEFINED ONLY WHEN
C                           IND IS GREATER THAN ZERO. SEE REMARKS.
C                           Y IS DESTROYED ON OUTPUT.
C                N      - INPUT LENGTH OF THE TIME SERIES.  N MUST BE
C                           EVENLY DIVISIBLE BY L. SEE REMARKS.
C                L      - INPUT PARAMETER USED TO SEGMENT THE TIME
C                           SERIES. L MUST BE A POWER OF TWO. SPECTRAL
C                           COMPUTATIONS ARE AT (L/2)+1 FREQUENCIES.
C                           SUGGESTED VALUES FOR L ARE 16, 32, 64, OR
C                           128.
C                IND    - INPUT CONTROL PARAMETER.
C                         IF IND IS GREATER THAN ZERO, THEN TWO TIME
C                           SERIES ARE INPUT AND PSX,PSY,AND XPS
C                           ARE OUTPUT.
C                         IF IND IS LESS THAN OR EQUAL TO ZERO, THEN
C                           ONLY ONE TIME SERIES IS INPUT AND PSX IS
C                           THE OUTPUT.
C                PSX    - OUTPUT VECTOR OF LENGTH (L/2)+1 CONTAINING
C                           THE SPECTRAL ESTIMATES OF X.
C                PSY    - OUTPUT VECTOR OF LENGTH (L/2)+1 CONTAINING
C                           THE SPECTRAL ESTIMATES OF Y, DEFINED FOR
C                           POSITIVE IND.
C                         NOTE THAT THE SPECTRAL ESTIMATES ARE TAKEN
C                           AT FREQUENCES (I-1)/L FOR I=1,2,...,(L/2)+1.
C                XPS    - OUTPUT VECTOR OF LENGTH L+2, DEFINED FOR IND
C                           POSITIVE. THE FIRST (L/2)+1 LOCATIONS CON-
C                           TAIN THE MAGNITUDE OF THE CROSS-SPECTRUM,
C                           AND THE LAST (L/2)+1 LOCATIONS CONTAIN THE
C                           PHASE OF THE CROSS-SPECTRUM, GIVEN IN
C                           FRACTIONS OF A CIRCLE (THAT IS, ON THE
C                           CLOSED INTERVAL (0,1)).
C                IWK    - INTEGER WORK VECTOR OF LENGTH M WHERE L=2**M.
C                WK     - REAL WORK VECTOR OF LENGTH L/2.
C                CWK    - COMPLEX WORK VECTOR OF LENGTH L+2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES L DOES NOT DIVIDE N EVENLY
C
C   REQD. IMSL ROUTINES - FFTCC,FFTRC,FFT2C,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  PRIOR TO CALLING FTFPS, THE MEANS OF TIME SERIES X
C                AND Y SHOULD BE REMOVED FROM THEIR RESPECTIVE SERIES
C                AND THE RESULTANT SERIES PADDED WITH ZEROS AS FOLLOWS;
C                A. FOR EACH TIME SERIES, COMPUTE THE MEAN OF THAT
C                   TIME SERIES AND CENTER THE TIME SERIES BY SUB-
C                   TRACTING THE MEAN FROM EACH ELEMENT OF THAT SERIES.
C                B. GIVEN INPUT L (AN INTEGER POWER OF TWO), THE
C                   LENGTH OF THE TIME SERIES (CALL IT LPS) IS USUALLY
C                   NOT EVENLY DIVISIBLE BY L. THE CENTERED TIME
C                   SERIES SHOULD BE PADDED WITH ZEROS ON THE RIGHT
C                   TO A LENGTH OF N, WHERE N IS THE FIRST INTEGER
C                   GREATER THAN OR EQUAL TO LPS THAT IS EVENLY
C                   DIVISIBLE BY L. N MAY BE COMPUTED ACCORDING TO
C                   THE FORMULA
C                   N = LPS + L - MOD(LPS,L)
C            2.  THE STABILITY OF THE ESTIMATES DEPENDS ON THE
C                AVERAGING PROCESS. THAT IS, THE GREATER THE NUMBER
C                OF SEGMENTS N/L, THE MORE STABLE THE RESULTING
C                ESTIMATES. N/L GREATER THAN 15 IS APPROPRIATE.
C                OTHERWISE, IMSL ROUTINE FTFREQ IS RECOMMENDED.
C            3.  THE OUTPUT IS RETURNED IN UNITS WHICH ARE THE
C                SQUARE OF THE DATA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTFPS   (X,Y,N,L,IND,PSX,PSY,XPS,IWK,WK,CWK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,L,IND,IWK(1),IER
      REAL               X(N),Y(1),PSX(1),PSY(1),XPS(1),WK(1)
      COMPLEX            CWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IPNF,J,K,KM1SL,LD2,LPK,LP1,M,NF
      REAL               C1,C2,PHASE,PI,PII,PI2,XIM1,XM,XSAVE,YSAVE
      DATA               PI/3.141593/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF ((N/L)*L .EQ. N) GO TO 5
      IER = 129
      GO TO 9000
    5 M = N / L
      LP1 = L + 1
      LD2 = L / 2
      NF = LD2 + 1
C                                  ZERO OUTPUT FOR SUMMING
      DO 10 I = 1,NF
         PSX(I) = 0.0
   10 CONTINUE
      IF (IND .LE. 0) GO TO 20
      DO 15 I = 1,NF
         PSY(I) = 0.0
         XPS(I) = 0.0
         XPS(I+NF) = 0.0
   15 CONTINUE
C                                  CALCULATE DATA WINDOW
   20 C1 = 0.5 * (L-1)
      C2 = 2.0 / (L+1)
      DO 25 I = 1,LD2
         XIM1 = I - 1
         WK(I) = 1.0 + C2*(XIM1-C1)
   25 CONTINUE
C                                  SEGMENT TIME SERIES
      DO 50 K = 1,M
         KM1SL = (K-1) * L
         LPK = LP1 + KM1SL
C                                  APPLY DATA WINDOW
         DO 30 I = 1,LD2
            J = I + KM1SL
            X(J) = X(J) * WK(I)
            J = LPK - I
            X(J) = X(J) * WK(I)
   30    CONTINUE
         CALL FFTRC (X(1+KM1SL),L,CWK,IWK,WK)
C                                  POWER SPECTRUM OF X
         DO 35 I = 1,LD2
            PSX(I) = PSX(I) + REAL(CWK(I))**2 + AIMAG(CWK(I))**2
   35    CONTINUE
         PSX(NF) = PSX(NF) + REAL(CWK(NF))**2
         IF (IND .LE. 0) GO TO 50
         DO 40 I = 1,LD2
            J = I + KM1SL
            Y(J) = Y(J) * WK(I)
            J = LPK - I
            Y(J) = Y(J) * WK(I)
   40    CONTINUE
         CALL FFTRC (Y(1+KM1SL),L,CWK(NF+1),IWK,WK)
C                                  SPECTRA OF Y AND X CROSS Y
         DO 45 I = 1,LD2
            IPNF = I + NF
            PSY(I) = PSY(I) + REAL(CWK(IPNF))**2 + AIMAG(CWK(IPNF))**2
            XPS(I) = XPS(I) + REAL(CWK(I))*REAL(CWK(IPNF)) +
     1               AIMAG(CWK(I))*AIMAG(CWK(IPNF))
            XPS(IPNF) = XPS(IPNF) + AIMAG(CWK(I))*REAL(CWK(IPNF)) -
     1                  REAL(CWK(I))*AIMAG(CWK(IPNF))
   45    CONTINUE
         IPNF = NF+NF
         PSY(NF) = PSY(NF) + REAL(CWK(IPNF))**2
         XPS(NF) = XPS(NF) + REAL(CWK(NF))*REAL(CWK(IPNF))
   50 CONTINUE
C                                  BOTH CWK(NF) AND CWK(NF+NF) ARE
C                                  USED AS REAL AVERAGING AND SCALING
C                                  OUTPUT
      XM = 3.0 / N
      DO 55 I = 1,NF
         PSX(I) = PSX(I) * XM
   55 CONTINUE
      IF (IND .LE. 0) GO TO 9005
      XPS(NF+NF) = 0.0
      PI2 = PI + PI
      PII = 1.0 / PI2
      DO 65 I = 1,NF
         IPNF = I + NF
         PSY(I) = PSY(I) * XM
         XSAVE = XPS(I) * XM
         YSAVE = XPS(IPNF) * XM
         XPS(I) = XSAVE*XSAVE + YSAVE*YSAVE
         IF (XSAVE.EQ.0.0 .AND. YSAVE.EQ.0.0) GO TO 60
         PHASE = ATAN2(YSAVE,XSAVE)
         IF (PHASE .LT. 0.0) PHASE = PHASE + PI2
         XPS(IPNF) = PHASE * PII
         GO TO 65
   60    XPS(IPNF) = 0.0
   65 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'FTFPS ')
 9005 RETURN
      END
