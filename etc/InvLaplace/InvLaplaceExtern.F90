MODULE InvLaplaceExtern
!
   USE Base, ONLY: i4, r4, r8, r16
!
   IMPLICIT NONE
!
   PRIVATE
!
   PUBLIC :: FS
!
CONTAINS
!
   FUNCTION FS(s) RESULT(ftval)
!
      IMPLICIT NONE
!
      COMPLEX(r8), INTENT(IN) :: s
!
      COMPLEX(r8) :: ftval, ai
!!!
!
      ai = (0._r8, 1._r8)
      ftval = 1.0_r8 / (s + 0.5_r8)
!
      RETURN
   END FUNCTION FS
!
END
