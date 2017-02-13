Module InvLaplace
!
   USE Base, ONLY: i4, r4, r8, r16
   USE InvLaplaceExtern, ONLY: FS
!
   IMPLICIT NONE
!
   REAL(r8),    PARAMETER :: PI = 3.14159265358979323846_r8
   INTEGER(i4), PARAMETER :: NGAUSS = 64
!
   PUBLIC  :: LINVSIDI
   PRIVATE :: GAUSS, ITGAUSS, EXTRAP
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE LINVSIDI
!
!  Computes inverse f(t) of Laplace transformed function F(s)
!
!        F(s) = \int_0^{\infty} \exp(-st) f(t) dt, where Re(s) > c
!
!  based on Sidi's modification of the algorithm by Dubner and Abate (1968).
!
!  References
!
!     [1] Dubner, H., and J. Abate, 1968: Numerical Inversion of Laplace
!         Transforms and the Finite Fourier Transform, JACM 15, 115-123.
!     [2] Cohen, A. M., 2007: Numerical Methods for Laplace Transform
!         Inversion, Springer-Verlag, New York (see Ch. 4.3.2: The Sidi mW-
!         Transformation for the Bromwich integral).
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE LINVSIDI(p, c, mm, t, y, ierr)
!
      IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!     p    : Accuracy of iterative Gauss integration
!     c    : A real number larger than real parts of complex numbers at which
!            the function "func" becomes singular
!     mm   : The number of dimension for results
!     t    : Array for t values
!     y    : Array for f(t)
!     ierr : Error code
!            = 1 for mm <= 0
!            = 2 for t(i) <= 0
!            = 3 for p <= 0
!-------------------------------------------------------------------------------
!
      INTEGER(i4), INTENT(IN) :: p
      REAL(r8),    INTENT(IN) :: c
      INTEGER(i4), INTENT(IN) :: mm
      REAL(r8),    DIMENSION(1:mm), INTENT(IN)  :: t
      REAL(r8),    DIMENSION(1:mm), INTENT(OUT) :: y
      INTEGER(i4), INTENT(OUT) :: ierr
!
      INTEGER(i4), PARAMETER :: M = 40
!
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x, w
      REAL(r8), DIMENSION(0:M)   :: psia, fa
      REAL(r8), DIMENSION(0:M+1) :: za
      REAL(r8), DIMENSION(0:M)   :: psib, fb
      REAL(r8), DIMENSION(0:M+1) :: zb
      REAL(r8) :: a, b, tt, xa, xb
      REAL(r8) :: valr, vali
      COMPLEX(r8) :: ai, val
      INTEGER(i4) :: n
      INTEGER(i4) :: i, k, l
      INTEGER(i4) :: ierr1, ierr2
!
      y(:) = 0._r8
!
      ai = (0._r8, 1._r8)
!
      ierr = 0
      IF (p <= 0) THEN
         ierr = 3; RETURN
      END IF
      IF (mm <= 0) THEN
         ierr = 1; RETURN
      END IF
!
      n = NGAUSS * p
      ALLOCATE(x(0:n))
      ALLOCATE(w(0:n))
!
      ierr = 0
!
      CALL GAUSS(NGAUSS, n, x, w, ierr1)
      IF (p > 0) CALL ITGAUSS(NGAUSS, n, p, x, w, ierr2)
!
      DO l = 1, mm
!
         tt = t(l)
!
         IF (tt <= 0._r8) THEN
            ierr = 2
            IF (ALLOCATED(x)) DEALLOCATE(x)
            IF (ALLOCATED(w)) DEALLOCATE(w)
            RETURN
         END IF
!
         DO i = 1, m+1
            za(i-1) = (i*PI/tt)
            zb(i-1) = za(i-1)
         END DO
!
         za(0) = za(0) / 2._r8
         k = 0
!
         fa(0) = 0._r8
         psia(k) = 0._r8 
         DO i = 1, n
            xa = k*za(0) + za(0)*x(i)
            val = FS(c + ai*xa)
            valr = REAL(val)
            valr = valr * COS(tt*xa)
            psia(k) = psia(k) + w(i)*valr
         END DO
         psia(k) = (PI*psia(k))/2._r8
!
         fb(0) = 0._r8
         psib(k) = 0._r8
         DO i = 1, n
            xb = k*zb(0) + zb(0)*x(i)
            val = FS(c + ai*xb)
            vali = AIMAG(val)
            vali = -vali*SIN(tt*xb)
            psib(k) = psib(k) + w(i)*vali
         END DO
         psib(k) = PI*psib(k)
!
         DO k = 1, m
            psia(k) = 0._r8
            psib(k) = 0._r8
            DO i = 1, n
               xa = (2*k-1)*za(0) + 2*za(0)*x(i)
               xb = k*zb(0) + zb(0)*x(i)
               val = FS(c + ai*xa) 
               valr = REAL(val)
               valr = valr*COS(tt*xa)
               psia(k) = psia(k) + w(i)*valr
               val = FS(c + ai*xb) 
               vali = AIMAG(val)
               vali = -vali*SIN(tt*xb)
               psib(k) = psib(k) + w(i)*vali
            END DO
            psia(k) = psia(k)*PI
            fa(k) = fa(k-1) + psia(k-1)
            psib(k) = psib(k)*PI
            fb(k) = fb(k-1) + psib(k-1)
         END DO
!
         CALL EXTRAP(m, psia, fa, za, a)
         a = 2*a/PI
         a = (a*EXP(c*tt))/tt
!
         CALL EXTRAP(m, psib, fb, zb, b)
         b = 2*b/PI
         b = (b*EXP(c*tt))/tt
!
         y(l) = (a+b)/2._r8 
!
      END DO
!
      IF (ALLOCATED(x)) DEALLOCATE(x)
      IF (ALLOCATED(w)) DEALLOCATE(w)
!
      RETURN
   END SUBROUTINE LINVSIDI
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE GAUSS
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE GAUSS(n, nn, x, w, ierr)
!
      IMPLICIT NONE
!
      INTEGER(i4),                   INTENT(IN)    :: n, nn
      REAL(r8),     DIMENSION(0:nn), INTENT(INOUT) :: x, w
      INTEGER(i4),                   INTENT(OUT)   :: ierr
!
      INTEGER(i4) :: i, j
!!!
!
      x( 1) = 0.999652520867886069728452812172818_r16
      w( 1) = 0.0008916403608482164736480395724898344_r16
      x( 2) = 0.998170058385977639673462250338200_r16
      w( 2) = 0.002073516630281233817643767864276530_r16
      x( 3) = 0.995506685738372160369691191721652_r16
      w( 3) = 0.003252228984489181428058680199989531_r16
      x( 4) = 0.991668126942312958465649651078416_r16
      w( 4) = 0.004423379913181973861515457329864311_r16
      x( 5) = 0.986663413894955481870926753676136_r16
      w( 5) = 0.005584069730065564409295246509604289_r16
      x( 6) = 0.980504399826026859459307060948579_r16
      w( 6) = 0.006731523948359321299030383342978218_r16
      x( 7) = 0.973205687429201408031240745673632_r16
      w( 7) = 0.007863015238012359660982997648769673_r16
      x( 8) = 0.964784586065969787910745077279613_r16
      w( 8) = 0.008975857887848671542522651000559224_r16
      x( 9) = 0.955261068539251402878190334004165_r16
      w( 9) = 0.01006741157676510468617015836427180_r16
      x(10) = 0.944657722997557052926702019136426_r16
      w(10) = 0.01113508690419162707964916519207727_r16
      x(11) = 0.932999699077046409880391692535079_r16
      w(11) = 0.01217635128435543666908877520453430_r16
      x(12) = 0.920314648126290181375845772347937_r16
      w(12) = 0.01318873485752732933584589631261280_r16
      x(13) = 0.906632657561398779870961669043152_r16
      w(13) = 0.01416983630712974161375565260011866_r16
      x(14) = 0.891986179471670703805110262606884_r16
      w(14) = 0.01511732853620123943398702990977421_r16
      x(15) = 0.876409953630265948305931887442847_r16
      w(15) = 0.01602896417742577679273375217394922_r16
      x(16) = 0.859940925085805413424470108915974_r16
      w(16) = 0.01690258091857080469578274105536263_r16
      x(17) = 0.842618156527116621281779185515688_r16
      w(17) = 0.01773610662844119190534657335762311_r16
      x(18) = 0.824482735627328669928880615996702_r16
      w(18) = 0.01852756427012002302020755090479165_r16
      x(19) = 0.805577677586196625124426485509274_r16
      w(19) = 0.01927507658930781456448124847340455_r16
      x(20) = 0.785947823101317017141939058329594_r16
      w(20) = 0.01997687056636017069332846306416800_r16
      x(21) = 0.765639732009947272829006951772228_r16
      w(21) = 0.02063128162131176430507814873681898_r16
      x(22) = 0.744701572853526478739263153510961_r16
      w(22) = 0.02123675756182679450366988395440869_r16
      x(23) = 0.723183008626732043992473857379458_r16
      w(23) = 0.02179186226466172668841393048686875_r16
      x(24) = 0.701135078981995801847883385630079_r16
      w(24) = 0.02229527908187828153006735501547241_r16
      x(25) = 0.678610079168834057975221307523101_r16
      w(25) = 0.02274581396370907223988549848563456_r16
      x(26) = 0.655661435995105478078756349280078_r16
      w(26) = 0.02314239829065720864797662461613062_r16
      x(27) = 0.632343581104383708186982086255010_r16
      w(27) = 0.02348409140810500866266314287729056_r16
      x(28) = 0.608711821870003542074824374494411_r16
      w(28) = 0.02377008285741515433114110347211160_r16
      x(29) = 0.584822210211996409018656814874135_r16
      w(29) = 0.02399969429822915386406308993567301_r16
      x(30) = 0.560731409648060277235188231746124_r16
      w(30) = 0.02417238111740147858488476357900890_r16
      x(31) = 0.536496560893899519724771470970169_r16
      w(31) = 0.02428773372075171346739953339198905_r16
      x(32) = 0.512175146331712216254477921426858_r16
      w(32) = 0.02434547850456986019168269536737497_r16
!
      j = n/2
      DO i = 1, j
         x(j+i) = 1.0_r8 - x(j-i+1)
         w(j+i) = w(j-i+1)
      END DO
!
      ierr = 0
!
      RETURN
   END SUBROUTINE GAUSS
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE ITGAUSS
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE ITGAUSS(n, nn, p, x, w, ierr)
!
      IMPLICIT NONE
!
      INTEGER(i4),                  INTENT(IN)    :: n, nn, p
      REAL(r8),    DIMENSION(0:nn), INTENT(INOUT) :: x, w
      INTEGER(i4),                  INTENT(OUT)   :: ierr
!
      REAL(r8) :: dp, h, hh
      INTEGER(i4) :: i, k, kk
!!!
!
      dp = DBLE(p)
!
      DO i = 1, n
         x(i) = x(i) / dp
         w(i) = w(i) / dp
      END DO
!
      h = 1._r8 / dp
!
      DO k = 2, p
         kk = k -1
         hh = kk * h
         DO i = 1, n
            x(kk*n+i) = hh + x(i)
            w(kk*n+i) = w(i)
         END DO 
      END DO
!
      RETURN
   END SUBROUTINE ITGAUSS
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE EXTRAP
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE EXTRAP(m, psi, f, z, a)
!
      IMPLICIT NONE
!
      INTEGER(i4),                 INTENT(IN)  :: m
      REAL(r8),    DIMENSION(0:m), INTENT(IN)  :: psi, f, z
      REAL(r8),                    INTENT(OUT) :: a
!
      REAL(r8), DIMENSION(0:m) :: p, q
      REAL(r8) :: b
      INTEGER(i4) :: j, k
!!!
!
      DO j = 0, m
         p(j) = f(j) /psi(j)
         q(j) = 1._r8/psi(j)
      END DO
!
      DO k = 1, m
         DO j = 0, m-k
            p(j) = p(j+1) - p(j)
            q(j) = q(j+1) - q(j)
               b = (1.0_r8/z(j+k)) - (1.0_r8/z(j))
            p(j) = p(j)/b
            q(j) = q(j)/b
         END DO
         a = p(0)/q(0)
      END DO
!
      RETURN
   END SUBROUTINE EXTRAP
!
END
