C   IMSL ROUTINE NAME   - VSAR
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SORTING OF MATRICES (WITH OPTIONS)
C
C   USAGE               - CALL VSAR (A,IA,NR,NC,IOP,KPOS,NK,IR,WK,IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT NR BY NC ARRAY
C                           ON INPUT, A CONTAINS THE ARRAY TO BE SORTED.
C                           ON OUTPUT, A CONTAINS THE SORTED ARRAY.
C                IA     - ROW DIMENSION OF ARRAY A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NR     - INPUT NUMBER OF ROWS OF A.
C                NC     - INPUT NUMBER OF COLUMNS OF A.
C                IOP    - INPUT VECTOR OF LENGTH 4 SPECIFYING OPTIONS.
C                           IOP(1) INDICATES WHETHER ROWS OR COLUMNS ARE
C                           THE RECORDS TO BE SORTED. WHEN
C                             IOP(1)=0, THE COLUMNS ARE THE RECORDS.
C                             IOP(1)=1, THE ROWS ARE THE RECORDS.
C                           IOP(2) INDICATES THE METHOD OF COMPARISON OF
C                           THE RECORDS. IOP(2)=I IMPLIES THE FOLLOWING
C                           METHOD, WHEN
C                             I=0, COMPONENTWISE, BY ALGEBRAIC VALUES.
C                             I=1, COMPONENTWISE, BY ABSOLUTE VALUES.
C                             I=2, BY THEIR L-2 NORM (EUCLIDEAN NORM.)
C                             I=3, BY THEIR L-1 NORM (SUM NORM.)
C                             I=4, BY THEIR L-INFINITY NORM (MAX NORM.)
C                           IOP(3)=I IMPLIES THE FOLLOWING ORDER, WHEN
C                             I=0, ASCENDING.
C                             I=1, DESCENDING.
C                           IOP(4)=I  IMPLIES THE FOLLOWING RESULTING
C                           ARRAY A AND RETURN INDEX VECTOR IR, WHEN
C                             I=0, A IS SORTED AND IR IS USED AS A WORK
C                               VECTOR.
C                             I=1, A IS SORTED, IR CONTAINS A RECORD OF
C                               THE PERMUTATIONS MADE ON THE COLUMNS OR
C                               ROWS OF A.
C                             I=2, A IS UNCHANGED AND IR IS THE SAME AS
C                               IN OPTION I=1 ABOVE (DETACHED KEY SORT.)
C                             I=3, A IS SORTED AND IR CONTAINS THE
C                               INDICES WITH RESPECT TO THE SORTED ORDER
C                               OF THE DISTINCT COLUMNS OR ROWS OF A.
C                             I=4, IMPLEMENTS OPTIONS I=1 AND I=3 ABOVE.
C                               IR CONTAINS THE PERMUTATIONS FOLLOWED BY
C                               THE INDICES OF THE DISTINCT RECORDS.
C                             I=5, A IS UNCHANGED AND IR IS THE SAME AS
C                               IN OPTION I=4 ABOVE (DETACHED KEY SORT.)
C                KPOS   - INPUT VECTOR OF LENGTH NK OF INDEX KEYS.
C                           KPOS CONTAINS INDICES OF SELECTED COLUMNS IF
C                             IOP(1)=1, OR SELECTED ROWS IF IOP(1)=0.
C                NK     - INPUT NUMBER OF KEYS.
C                IR     - OUTPUT INDEX VECTOR OF LENGTH IIR, WHERE
C                           IF IOP(4)=0, IIR=(1-IOP(1))*NR + IOP(1)*NC,
C                           IF IOP(4)=1,2,OR 3 IIR = NR+NC,
C                           IF IOP(4)=4,OR 5 IIR =
C                             (1-IOP(1))*MAX(2*NC+1,NC+NR) +
C                             IOP(1)*MAX(2*NR+1,NC+NR)
C                           IF IOP(4)=1,2,4, OR 5,
C                             ON OUTPUT, IR CONTAINS, IN THE FIRST NC
C                               OR NR POSITIONS, A RECORD OF THE
C                               PERMUTATIONS MADE ON THE COLUMNS OR
C                               ROWS OF ARRAY A. SEE REMARKS.
C                           IF IOP(4)=3, IR CONTAINS THE INDICES OF THE
C                             DISTINCT COLUMNS OR ROWS, WITH RESPECT TO
C                             SORTED ORDER. THE NUMBER OF THESE IS IN
C                             IR(NC+1) IF IOP(1)=0 OR IN IR(NR+1) IF
C                             IOP(1)=1.
C                           IF IOP(4)=5 IR CONTAINS THE PERMUTATION
C                             VECTOR FOLLOWED BY THE DISTINCT RECORDS
C                             INDICES. THE NUMBER OF DISTINCT RECORDS IS
C                             IN IR(2*NC+1) FOR COLUMNS, OR
C                             IN IR(2*NR+1) FOR ROWS. SEE REMARKS.
C                WK     - WORK VECTOR OF LENGTH IIW, WHERE IIW =
C                           (1-IOP(1))*MAX(NC,2*NR)+IOP(1)*MAX(NR,2*NC).
C                IER    - OUTPUT ERROR PARAMETER.
C                           TERMINAL ERROR
C                             IER=129 INDICATES THAT IA, NR, NC OR NK
C                               IS .LE. 0, OR IOP IS NOT SET PROPERLY,
C                               OR IER=65 OCCURS FOR ALL KEYS.
C                               NO SORTING IS PERFORMED.
C                           WARNING WITH FIX ERROR
C                             IER=65 INDICATES THAT, FOR SOME J
C                               KPOS(J).LE.0 OR KPOS(J).GT.NC (OR NR),
C                               IN WHICH CASE KPOS(J) IS IGNORED, AND/OR
C                               NR IS .GT. IA, IN WHICH CASE NR IS
C                               LIMITED TO IA.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,VBLA=ISAMAX,VBLA=SASUM,
C                           VBLA=SCOPY,VBLA=SNRM2,VBLA=SSWAP,VNABSX,
C                           VNAXPB,VNDIST,VNINI,VNPINV,VNSCPY,VNSWAP,
C                           VSARER,VSARX,VSCMP,VSCMPA,VSCMPE,VSENA,
C                           VSIRA,VSIRAQ,VSIRE,VSIREQ,VSIRU,VSIRUQ,
C                           VSORA,VSORAQ,VSORE,VSOREQ,VSORU,VSORUQ,
C                           VSRTP,VSSWAM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      1.  WHEN A IS SORTED BY COLUMNS AND BY ALGEBRAIC VALUES
C                    (IOP(2)=0), IN ASCENDING ORDER, THE RESULTING
C                    ARRAY A IS SUCH THAT, FOR K=1,2,...,NK-1
C                    1. A(1,KPOS(K)).LE.A(1,KPOS(K+1)), AND
C                    2. FOR J=2,..,NK, IF A(I,KPOS(K)).EQ.A(I,KPOS(K+1))
C                       FOR I=1,2,...,J-1, THEN
C                              A(J,KPOS(K)).LE.A(J,KPOS(K+1)).
C                    WHEN IOP(2)=1, THE ABSOLUTE VALUES ARE COMPARED
C                    INSTEAD.
C                2.  WHEN IOP(4)=1, 2, 4, OR 5, THE ACTION OF VECTOR IR
C                    IS AS FOLLOWS. LET AU REPRESENT THE ARRAY A IN ITS
C                    UNSORTED INPUT ORDER, AND LET AO BE THE ARRAY A
C                    SORTED. THEN, ON OUTPUT IR CONTAINS A PERMUTATION
C                    OF THE FIRST NR OR NC INTEGERS SUCH THAT,
C                      AO(*,I)=AU(*,IR(I)) IF IOP(1)=0 AND,
C                      AO(I,*)=AU(IR(I),*) IF IOP(1)=1.
C                3.  WHEN THE COLUMNS ARE SORTED (IOP(1) = 0), THE
C                    RETURN VECTOR IR IS ORGANIZED SO THAT, WHEN
C                      IOP(4)=1,2,4 OR 5, IR(1),...,IR(NC) CONTAIN
C                        A RECORD OF THE PERMUTATIONS ON THE COLUMNS.
C                      IOP(4)=3, IR(1),...,IR(IR(NC+1)) CONTAIN THE
C                        INDICES OF THE DISTINCT COLUMNS, AND IR(NC+1)
C                        CONTAIN THE NUMBER OF DISTINCT COLUMNS.
C                      IOP(4)=4 OR 5, IR(NC+1),...,IR(NC + IR(2*NC+1))
C                        CONTAIN THE INDICES OF THE DISTINCT COLUMNS,
C                        AND IR(2*NC+1) CONTAINS THE NUMBER OF DISTINCT
C                        COLUMNS. (THE ROW CASE IS ANALOGOUS.)
C                4.  WHEN THE COLUMNS ARE SORTED (IOP(1)=1) AND IOP(4)=5
C                    (A IS UNCHANGED), IR(NC+J) CONTAINS THE INDEX
C                    OF THE J-TH DISTINCT COLUMN WITH RESPECT TO THE
C                    SORTED ORDER, AND THIS COLUMN CAN BE ACCESSED AT
C                    AU(*,IR(IR(NC+J))). (THE ROW CASE IS ANALOGOUS.)
C                5.  NO DUPLICATIONS ARE PERMISSIBLE IN VECTOR KPOS.
C                    IF THIS CONDITION IS NOT MET, THE ROUTINE MAY NOT
C                    PERFORM CORRECTLY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSAR (A,IA,NR,NC,IOP,KPOS,NK,IR,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,NK,IER,IOP(1),KPOS(1),IR(1)
      REAL               A(IA,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ID,IE,IK,LC,LH,LHM1,LIR,LJR,LK,LR,N,ND
C                                  FIRST EXECUTABLE STATEMENT
C                                  ERROR CHECKING AND INITIALIZATION
      CALL VSARER(IA,NR,NC,IOP,KPOS,NK,LR,LC,LK,IR,LIR,LJR,LH,IK,ID,IER)
      IF (IER.EQ.129) GO TO 9000
      CALL VSARX(A,IA,LR,LC,1-IOP(1),IR(IK),WK)
      N = 5*IOP(1)+IOP(2)+1
      GO TO (5, 10, 15, 15, 15, 20, 25, 30, 30, 30), N
C                                  SORT COLUMNS
    5 IF (IOP(4).EQ.0) CALL VSORA(A,IA,LR,LC,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSORAQ(A,IA,LR,LC,LK,IR,WK,IE)
      GO TO 35
   10 IF (IOP(4).EQ.0) CALL VSORU(A,IA,LR,LC,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSORUQ(A,IA,LR,LC,LK,IR,WK,IE)
      GO TO 35
   15 IF (IOP(4).EQ.0) CALL VSORE(A,IA,LR,LC,IOP(2)-1,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSOREQ(A,IA,LR,LC,IOP(2)-1,LK,IR,WK,IE)
      GO TO 35
C                                  SORT ROWS
   20 IF (IOP(4).EQ.0) CALL VSIRA(A,IA,LR,LC,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSIRAQ(A,IA,LR,LC,LK,IR,WK,IE)
      GO TO 35
   25 IF (IOP(4).EQ.0) CALL VSIRU(A,IA,LR,LC,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSIRUQ(A,IA,LR,LC,LK,IR,WK,IE)
      GO TO 35
   30 IF (IOP(4).EQ.0) CALL VSIRE(A,IA,LR,LC,IOP(2)-1,LK,WK,IE)
      IF (IOP(4).GT.0) CALL VSIREQ(A,IA,LR,LC,IOP(2)-1,LK,IR,WK,IE)
C                                  DESCENDING/DETACHED/DISTINCT CASES
   35 CALL VNPINV(IR(IK),LJR,WK)
      CALL VSARX(A,IA,LR,LC,1-IOP(1),IR(IK),WK)
      IF (IOP(4).GE.3) CALL VNDIST(IR,LIR,ID,ND)
      IF (IOP(4).EQ.1 .OR. IOP(4).EQ.2) CALL VNABSX(LIR,IR,1)
      IF (IOP(4)*IOP(3).NE.2 .AND. IOP(4).NE.5 .AND. IOP(3).NE.0) CALL
     *VSSWAM(A,IA,LR,LC,IOP(1),1,1)
      IE = -1
      IF (IOP(4)*IOP(3).GE.3) CALL VNAXPB(ND,IE,IR(ID),1,LIR+1,0)
      IF (IOP(4).EQ.2 .OR. IOP(4).EQ.5) CALL VNPINV(IR,LIR,WK)
      IF (IOP(4).EQ.2 .OR. IOP(4).EQ.5) CALL VSARX(A,IA,LR,LC,IOP(1),IR,
     *WK)
      IF (IOP(4).EQ.2 .OR. IOP(4).EQ.5) CALL VNPINV(IR,LIR,WK)
      IF (IOP(4)*IOP(3).LT.1 .OR. IOP(4).EQ.3) GO TO 9000
      IF (MOD(NC,2) .NE. 0 .AND. IOP(1) .EQ. 0) GO TO 40
      IF (MOD(NR,2) .NE. 0 .AND. IOP(1) .EQ. 1) GO TO 40
      LHM1 = LH + 1
      GO TO 45
   40 LHM1 = LH
   45 CALL VNSWAP(LH,IR,1,IR(LHM1),-1)
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HVSAR  )
 9005 RETURN
      END
