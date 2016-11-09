C   IMSL ROUTINE NAME   - UHELP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - DISPLAY METHODS OF OBTAINING INFORMATION ON
C                           IMSL CONVENTIONS REGARDING VARIOUS SUBJECTS
C                           AND PROVIDE A MEANS FOR INDIVIDUAL SITES TO
C                           SUPPLY USERS WITH SITE SPECIFIC INFORMATION
C
C   USAGE               - CALL UHELP
C
C   REQD. IMSL ROUTINES - UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF THE USER DESIRES THAT OUTPUT BE WRITTEN TO A DEVICE
C                OTHER THAN THE STANDARD OUTPUT DEVICE, HE MUST FIRST
C                CALL IMSL ROUTINE UGETIO TO RESET THE OUTPUT DEVICE.
C                SEE THE UGETIO DOCUMENT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UHELP
C
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IOPT,NIN,NOUT
C                                  FIRST EXECUTABLE STATEMENT
      IOPT = 1
      CALL UGETIO (IOPT,NIN,NOUT)
      WRITE(NOUT,5)
      WRITE(NOUT,10)
    5 FORMAT(1X,                                             65HDETAILED
     * INFORMATION CONCERNING IMSL CONVENTIONS ON THE FOLLOWING,/,1X,63H
     *SUBJECTS IS AVAILABLE BY CALLING THE INDICATED IMSL SUBROUTINE.,/,
     */,1X,57HSUBJECT                                      SUBROUTINE  ,
     */,1X,57H-------                                      ----------  ,
     */,1X,57HDOCUMENTATION CONVENTIONS AND NOTATION       CALL UHELP1 ,
     */,1X,57HINPUT/OUTPUT CONVENTIONS                     CALL UHELP2 ,
     */,1X,57HERROR DETECTING FACILITIES                   CALL UHELP3 ,
     */,1X,57HMATRIX/VECTOR STORAGE MODES                  CALL UHELP4 )
   10 FORMAT(
     *  1X,57H                                                         ,
     */,1X,57HEXAMPLE                                                  ,
     */,1X,57H                                                         ,
     */,1X,57HTO OBTAIN INFORMATION REGARDING IMSL INPUT AND OUTPUT    ,
     */,1X,57HCONVENTIONS, ONE WOULD RUN THE FOLLOWING PROGRAM.        ,
     */,1X,57H                                                         ,
     */,1X,57H      CALL UHELP2                                        ,
     */,1X,57H      STOP                                               ,
     */,1X,13H      END    )
      RETURN
      END

