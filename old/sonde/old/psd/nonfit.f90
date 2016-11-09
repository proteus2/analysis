!2345
!------------------------------------------------------------------------------------
!
!   This module NONFIT is based on MINPACK for minimization of the
!   sum of the squares of m nonlinear functions in n variables by a modification of    
!   Leveberg-Marquardt algorithm.
!   Documentation for double precision version MINPACK programs can be found in 
!   /export1/f90lib/NONFIT.
!
!------------------------------------------------------------------------------------

    module nonfit

    contains

!------------------------------------------------------------------------------------
!
!   SUBROUTINE LMDERDRV
!
!   LMDERDRV is LMDER1 in original MINPACK.
!
!------------------------------------------------------------------------------------

    subroutine lmderdrv(fcn,m,ux,uy,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)

    integer m,n,ldfjac,info,lwa
    integer ipvt(n)
    double precision tol
    double precision ux(m),uy(m)
    double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)

    external fcn

    integer maxfev,mode,nfev,njev,nprint
    double precision factor,ftol,gtol,xtol,zero
    data factor,zero /1.0d2,0.0d0/
    info = 0
!
!   check the input parameters for errors.
!
    if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m .or. tol .lt. zero  &
        .or. lwa .lt. 5*n + m) go to 10
!
!   call lmder.
!
    maxfev = 100*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    mode = 1
    nprint = 0
    call lmder(fcn,m,ux,uy,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev,  &
               wa(1),mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1),  &
               wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
    if (info .eq. 8) info = 4
10  continue
    return
!
!   last card of subroutine lmder1.
!
    end subroutine lmderdrv

!------------------------------------------------------------------------------------
!
!   SUBROUTINE LMDIFDRV
!
!   LMDIFDRV is LMDIF1 in original MINPACK.
!
!------------------------------------------------------------------------------------

    subroutine lmdifdrv(fcn,m,ux,uy,n,x,fvec,tol,info,iwa,wa,lwa)

    integer m,n,info,lwa
    integer iwa(n)
    double precision ux(m), uy(m)
    double precision tol
    double precision x(n),fvec(m),wa(lwa)
    external fcn

    integer maxfev,mode,mp5n,nfev,nprint
    double precision epsfcn,factor,ftol,gtol,xtol,zero
    data factor,zero /1.0d2,0.0d0/
    info = 0
!
!   check the input parameters for errors.
!
    if (n .le. 0 .or. m .lt. n .or. tol .lt. zero .or. lwa .lt. m*n + 5*n + m) go to 10
!
!   call lmdif.
!
    maxfev = 200*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    epsfcn = zero
    mode = 1
    nprint = 0
    mp5n = m + 5*n
    call lmdif(fcn   ,m     ,ux    ,uy    ,n     ,x     ,fvec   ,ftol     ,  &
               xtol  ,gtol  ,maxfev,epsfcn,wa(1) ,mode  ,factor ,nprint   ,  &
               info  ,nfev  ,wa(mp5n+1)   ,m     ,iwa   ,wa(n+1),wa(2*n+1),  &
               wa(3*n+1),wa(4*n+1),wa(5*n+1))
    if (info .eq. 8) info = 4
10  continue
    return
!
!   last card of subroutine lmdifdrv.
!
    end subroutine lmdifdrv

!------------------------------------------------------------------------------------
!
!   SUBROUTINE LMDER
! 
!------------------------------------------------------------------------------------

    subroutine lmder(fcn,m,ux,uy,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,  &
                     maxfev,diag,mode,factor,nprint,info,nfev,njev,  &
                     ipvt,qtf,wa1,wa2,wa3,wa4)
!
    integer m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
    integer ipvt(n)
    double precision ftol,xtol,gtol,factor
    double precision ux(m),uy(m)
    double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n)
    double precision wa1(n),wa2(n),wa3(n),wa4(m)
!
    integer i,iflag,iter,j,l
    double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm
    double precision one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio
    double precision sum,temp,temp1,temp2,xnorm,zero
!   double precision dpmpar,enorm
    data one,p1,p5,p25,p75,p0001,zero /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!   epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    info = 0
    iflag = 0
    nfev = 0
    njev = 0
!
!   check the input parameters for errors.
!
    if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m   &
        .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero  &
        .or. maxfev .le. 0 .or. factor .le. zero) go to 300
    if (mode .ne. 2) go to 20
    do 10 j = 1, n
       if (diag(j) .le. zero) go to 300
10  continue
20  continue
!
!   evaluate the function at the starting point
!   and calculate its norm.
!
    iflag = 1
    call fcn(m,ux,uy,n,x,fvec,fjac,ldfjac,iflag)
    nfev = 1
    if (iflag .lt. 0) go to 300
    fnorm = enorm(m,fvec)
!
!   initialize levenberg-marquardt parameter and iteration counter.
!
    par = zero
    iter = 1
!
!   beginning of the outer loop.
!
30  continue
!
!   calculate the jacobian matrix.
!
    iflag = 2
    call fcn(m,ux,uy,n,x,fvec,fjac,ldfjac,iflag)
    njev = njev + 1
    if (iflag .lt. 0) go to 300
!
!   if requested, call fcn to enable printing of iterates.
!
    if (nprint .le. 0) go to 40
    iflag = 0

    if (mod(iter-1,nprint) .eq. 0) call fcn(m,ux,uy,n,x,fvec,fjac,ldfjac,iflag)
    if (iflag .lt. 0) go to 300
40  continue
!
!   compute the qr factorization of the jacobian.
!
    call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!   on the first iteration and if mode is 1, scale according
!   to the norms of the columns of the initial jacobian.
!
    if (iter .ne. 1) go to 80
    if (mode .eq. 2) go to 60
    do 50 j = 1, n
      diag(j) = wa2(j)
      if (wa2(j) .eq. zero) diag(j) = one
50  continue
60  continue
!
!   on the first iteration, calculate the norm of the scaled x
!   and initialize the step bound delta.
!
    do 70 j = 1, n
      wa3(j) = diag(j)*x(j)
70  continue
    xnorm = enorm(n,wa3)
    delta = factor*xnorm
    if (delta .eq. zero) delta = factor
80  continue
!
!   form (q transpose)*fvec and store the first n components in qtf.
!
    do 90 i = 1, m
      wa4(i) = fvec(i)
90  continue
    do 130 j = 1, n
      if (fjac(j,j) .eq. zero) go to 120
      sum = zero
      do 100 i = j, m
        sum = sum + fjac(i,j)*wa4(i)
100   continue
      temp = -sum/fjac(j,j)
      do 110 i = j, m
        wa4(i) = wa4(i) + fjac(i,j)*temp
110   continue
120   continue
      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)
130 continue
!
!   compute the norm of the scaled gradient.
!
    gnorm = zero
    if (fnorm .eq. zero) go to 170
    do 160 j = 1, n
      l = ipvt(j)
      if (wa2(l) .eq. zero) go to 150
      sum = zero
      do 140 i = 1, j
        sum = sum + fjac(i,j)*(qtf(i)/fnorm)
140   continue
      gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
150   continue
160 continue
170 continue
!
!   test for convergence of the gradient norm.
!
    if (gnorm .le. gtol) info = 4
    if (info .ne. 0) go to 300
!
!   rescale if necessary.
!
    if (mode .eq. 2) go to 190
    do 180 j = 1, n
      diag(j) = dmax1(diag(j),wa2(j))
180 continue
190 continue
!
!   beginning of the inner loop.
!
200 continue
!
!   determine the levenberg-marquardt parameter.
!
    call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4)
!
!   store the direction p and x + p. calculate the norm of p.
!
    do 210 j = 1, n
      wa1(j) = -wa1(j)
      wa2(j) = x(j) + wa1(j)
      wa3(j) = diag(j)*wa1(j)
210 continue
    pnorm = enorm(n,wa3)
!
!   on the first iteration, adjust the initial step bound.
!
    if (iter .eq. 1) delta = dmin1(delta,pnorm)
!
!   evaluate the function at x + p and calculate its norm.
!
    iflag = 1
    call fcn(m,ux,uy,n,wa2,wa4,fjac,ldfjac,iflag)
    nfev = nfev + 1
    if (iflag .lt. 0) go to 300
    fnorm1 = enorm(m,wa4)
!
!   compute the scaled actual reduction.
!
    actred = -one
    if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!   compute the scaled predicted reduction and
!   the scaled directional derivative.
!
    do 230 j = 1, n
      wa3(j) = zero
      l = ipvt(j)
      temp = wa1(l)
      do 220 i = 1, j
        wa3(i) = wa3(i) + fjac(i,j)*temp
220   continue
230 continue
    temp1 = enorm(n,wa3)/fnorm
    temp2 = (dsqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)
!
!   compute the ratio of the actual to the predicted
!   reduction.
!
    ratio = zero
    if (prered .ne. zero) ratio = actred/prered
!
!   update the step bound.
!
    if (ratio .gt. p25) go to 240
    if (actred .ge. zero) temp = p5
    if (actred .lt. zero) temp = p5*dirder/(dirder + p5*actred)
    if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
    delta = temp*dmin1(delta,pnorm/p1)
    par = par/temp
    go to 260
240 continue
    if (par .ne. zero .and. ratio .lt. p75) go to 250
    delta = pnorm/p5
    par = p5*par
250 continue
260 continue
!
!   test for successful iteration.
!
    if (ratio .lt. p0001) go to 290
!
!   successful iteration. update x, fvec, and their norms.
!
    do 270 j = 1, n
      x(j) = wa2(j)
      wa2(j) = diag(j)*x(j)
270 continue
    do 280 i = 1, m
      fvec(i) = wa4(i)
280 continue
    xnorm = enorm(n,wa2)
    fnorm = fnorm1
    iter = iter + 1
290 continue
!
!   tests for convergence.
!
    if (dabs(actred) .le. ftol .and. prered .le. ftol .and. p5*ratio .le. one) info = 1
    if (delta .le. xtol*xnorm) info = 2
    if (dabs(actred) .le. ftol .and. prered .le. ftol  &
              .and. p5*ratio .le. one .and. info .eq. 2) info = 3
    if (info .ne. 0) go to 300
!
!   tests for termination and stringent tolerances.
!
    if (nfev .ge. maxfev) info = 5
    if (dabs(actred) .le. epsmch .and. prered .le. epsmch .and. p5*ratio .le. one) info = 6
    if (delta .le. epsmch*xnorm) info = 7
    if (gnorm .le. epsmch) info = 8
    if (info .ne. 0) go to 300
!
!   end of the inner loop. repeat if iteration unsuccessful.
!
    if (ratio .lt. p0001) go to 200
!
!   end of the outer loop.
!
    go to 30
300 continue
!
!   termination, either normal or user imposed.
!
    if (iflag .lt. 0) info = iflag
    iflag = 0
    if (nprint .gt. 0) call fcn(m,ux,uy,n,x,fvec,fjac,ldfjac,iflag)
    return
!
!   last card of subroutine lmder.
!
    end subroutine lmder

!------------------------------------------------------------------------------------
!
!   SUBROUTINE LMDIF
!
!------------------------------------------------------------------------------------

    subroutine lmdif(fcn   ,m     ,ux    ,uy    ,n     ,x     ,fvec  ,ftol  ,  &
                     xtol  ,gtol  ,maxfev,epsfcn,diag  ,mode  ,factor,nprint,  &
                     info  ,nfev  ,fjac  ,ldfjac,ipvt  ,qtf   ,wa1   ,wa2   ,  &
                     wa3   ,wa4   )

    integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
    integer ipvt(n)
    double precision ux(m), uy(m)
    double precision ftol,xtol,gtol,epsfcn,factor
    double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n)
    double precision wa1(n),wa2(n),wa3(n),wa4(m)

    external fcn

    integer i,iflag,iter,j,l
    double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm
    double precision one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio
    double precision sum,temp,temp1,temp2,xnorm,zero
!   double precision dpmpar,enorm
    data one,p1,p5,p25,p75,p0001,zero /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/ 
!
!   epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    info = 0
    iflag = 0
    nfev = 0
!
!   check the input parameters for errors.
!
    if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m   &
        .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero  &
        .or. maxfev .le. 0 .or. factor .le. zero) go to 300
    if (mode .ne. 2) go to 20
    do 10 j = 1, n
      if (diag(j) .le. zero) go to 300
10  continue
20  continue
!
!   evaluate the function at the starting point
!   and calculate its norm.
!
    iflag = 1
    call fcn(m,ux,uy,n,x,fvec,iflag)
    nfev = 1
    if (iflag .lt. 0) go to 300
    fnorm = enorm(m,fvec)
!
!   initialize levenberg-marquardt parameter and iteration counter.
!
    par = zero
    iter = 1
!
!   beginning of the outer loop.
!
30  continue
!
!   calculate the jacobian matrix.
!
    iflag = 2
    call fdjac2(fcn,m,ux,uy,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
    nfev = nfev + n
    if (iflag .lt. 0) go to 300
!
!   if requested, call fcn to enable printing of iterates.
!
    if (nprint .le. 0) go to 40
    iflag = 0
    if (mod(iter-1,nprint) .eq. 0) call fcn(m,ux,uy,n,x,fvec,iflag)
    if (iflag .lt. 0) go to 300
40  continue
!
!   compute the qr factorization of the jacobian.
!
    call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!   on the first iteration and if mode is 1, scale according
!   to the norms of the columns of the initial jacobian.
!
    if (iter .ne. 1) go to 80
    if (mode .eq. 2) go to 60
    do 50 j = 1, n
      diag(j) = wa2(j)
      if (wa2(j) .eq. zero) diag(j) = one
50  continue
60  continue
!
!   on the first iteration, calculate the norm of the scaled x
!   and initialize the step bound delta.
!
    do 70 j = 1, n
      wa3(j) = diag(j)*x(j)
70  continue
    xnorm = enorm(n,wa3)
    delta = factor*xnorm
    if (delta .eq. zero) delta = factor
80  continue
!
!   form (q transpose)*fvec and store the first n components in qtf.
!
    do 90 i = 1, m
      wa4(i) = fvec(i)
90  continue
    do 130 j = 1, n
      if (fjac(j,j) .eq. zero) go to 120
      sum = zero
      do 100 i = j, m
        sum = sum + fjac(i,j)*wa4(i)
100   continue
      temp = -sum/fjac(j,j)
      do 110 i = j, m
        wa4(i) = wa4(i) + fjac(i,j)*temp
110   continue
120   continue
      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)
130 continue
!
!   compute the norm of the scaled gradient.
!
    gnorm = zero
    if (fnorm .eq. zero) go to 170
    do 160 j = 1, n
      l = ipvt(j)
      if (wa2(l) .eq. zero) go to 150
      sum = zero
      do 140 i = 1, j
        sum = sum + fjac(i,j)*(qtf(i)/fnorm)
140   continue
      gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
150   continue
160   continue
170 continue
!
!   test for convergence of the gradient norm.
!
    if (gnorm .le. gtol) info = 4
    if (info .ne. 0) go to 300
!
!   rescale if necessary.
!
    if (mode .eq. 2) go to 190
    do 180 j = 1, n
      diag(j) = dmax1(diag(j),wa2(j))
180 continue
190 continue
!
!   beginning of the inner loop.
!
200 continue
!
!   determine the levenberg-marquardt parameter.
!
    call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4)
!
!   store the direction p and x + p. calculate the norm of p.
!
    do 210 j = 1, n
      wa1(j) = -wa1(j)
      wa2(j) = x(j) + wa1(j)
      wa3(j) = diag(j)*wa1(j)
210 continue
    pnorm = enorm(n,wa3)
!
!   on the first iteration, adjust the initial step bound.
!
    if (iter .eq. 1) delta = dmin1(delta,pnorm)
!
!   evaluate the function at x + p and calculate its norm.
!
    iflag = 1
    call fcn(m,ux,uy,n,wa2,wa4,iflag)
    nfev = nfev + 1
    if (iflag .lt. 0) go to 300
    fnorm1 = enorm(m,wa4)
!
!   compute the scaled actual reduction.
!
    actred = -one
    if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!   compute the scaled predicted reduction and
!   the scaled directional derivative.
!
    do 230 j = 1, n
      wa3(j) = zero
      l = ipvt(j)
      temp = wa1(l)
      do 220 i = 1, j
        wa3(i) = wa3(i) + fjac(i,j)*temp
220   continue
230 continue
    temp1 = enorm(n,wa3)/fnorm
    temp2 = (dsqrt(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)
!
!   compute the ratio of the actual to the predicted
!
    ratio = zero
    if (prered .ne. zero) ratio = actred/prered
!
!   update the step bound.
!
    if (ratio .gt. p25) go to 240
    if (actred .ge. zero) temp = p5
    if (actred .lt. zero) temp = p5*dirder/(dirder + p5*actred)
    if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
    delta = temp*dmin1(delta,pnorm/p1)
    par = par/temp
    go to 260
240 continue
    if (par .ne. zero .and. ratio .lt. p75) go to 250
    delta = pnorm/p5
    par = p5*par
250 continue
260 continue
!
!   test for successful iteration.
!
    if (ratio .lt. p0001) go to 290
!
!   successful iteration. update x, fvec, and their norms.
!
    do 270 j = 1, n
      x(j) = wa2(j)
      wa2(j) = diag(j)*x(j)
270 continue
    do 280 i = 1, m
      fvec(i) = wa4(i)
280 continue
    xnorm = enorm(n,wa2)
    fnorm = fnorm1
    iter = iter + 1
290 continue
!
!   tests for convergence.
!
    if (dabs(actred) .le. ftol .and. prered .le. ftol   &
        .and. p5*ratio .le. one) info = 1
    if (delta .le. xtol*xnorm) info = 2
    if (dabs(actred) .le. ftol .and. prered .le. ftol  &
        .and. p5*ratio .le. one .and. info .eq. 2) info = 3
    if (info .ne. 0) go to 300
!
!   tests for termination and stringent tolerances.
!
    if (nfev .ge. maxfev) info = 5
    if (dabs(actred) .le. epsmch .and. prered .le. epsmch .and. p5*ratio .le. one) info = 6
    if (delta .le. epsmch*xnorm) info = 7
    if (gnorm .le. epsmch) info = 8
    if (info .ne. 0) go to 300
!
!   end of the inner loop. repeat if iteration unsuccessful.
!
    if (ratio .lt. p0001) go to 200
!
!   end of the outer loop.
!
    go to 30
300 continue
!
!   termination, either normal or user imposed.
!
    if (iflag .lt. 0) info = iflag
    iflag = 0
    if (nprint .gt. 0) call fcn(m,ux,uy,n,x,fvec,iflag)
    return
!
!   last card of subroutine lmdif.
!
    end subroutine lmdif

!------------------------------------------------------------------------------------
!
!   SUBROUTINE LMPAR
!
!------------------------------------------------------------------------------------

    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)

    integer n,ldr
    integer ipvt(n)
    double precision delta,par
    double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),wa2(n)

    integer i,iter,j,jm1,jp1,k,l,nsing
    double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001
    double precision sum,temp,zero
!   double precision dpmpar,enorm
    data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
!
!   dwarf is the smallest positive magnitude.
!
    dwarf = dpmpar(2)
!
!   compute and store in x the gauss-newton direction. if the
!   jacobian is rank-deficient, obtain a least squares solution.
!
    nsing = n
    do 10 j = 1, n
       wa1(j) = qtb(j)
       if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
       if (nsing .lt. n) wa1(j) = zero
 10    continue
    if (nsing .lt. 1) go to 50
    do 40 k = 1, nsing
       j = nsing - k + 1
       wa1(j) = wa1(j)/r(j,j)
       temp = wa1(j)
       jm1 = j - 1
       if (jm1 .lt. 1) go to 30
       do 20 i = 1, jm1
          wa1(i) = wa1(i) - r(i,j)*temp
 20       continue
 30    continue
 40    continue
 50 continue
    do 60 j = 1, n
       l = ipvt(j)
       x(l) = wa1(j)
 60    continue
!
!   initialize the iteration counter.
!   evaluate the function at the origin, and test
!   for acceptance of the gauss-newton direction.
!
    iter = 0
    do 70 j = 1, n
       wa2(j) = diag(j)*x(j)
 70    continue
    dxnorm = enorm(n,wa2)
    fp = dxnorm - delta
    if (fp .le. p1*delta) go to 220
!
!   if the jacobian is not rank deficient, the newton
!   step provides a lower bound, parl, for the zero of
!   the function. otherwise set this bound to zero.
!
    parl = zero
    if (nsing .lt. n) go to 120
    do 80 j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
 80    continue
    do 110 j = 1, n
       sum = zero
       jm1 = j - 1
       if (jm1 .lt. 1) go to 100
       do 90 i = 1, jm1
          sum = sum + r(i,j)*wa1(i)
 90       continue
100    continue
       wa1(j) = (wa1(j) - sum)/r(j,j)
110    continue
    temp = enorm(n,wa1)
    parl = ((fp/delta)/temp)/temp
120 continue
!
!   calculate an upper bound, paru, for the zero of the function.
!
    do 140 j = 1, n
       sum = zero
       do 130 i = 1, j
          sum = sum + r(i,j)*qtb(i)
130    continue
       l = ipvt(j)
       wa1(j) = sum/diag(l)
140 continue
    gnorm = enorm(n,wa1)
    paru = gnorm/delta
    if (paru .eq. zero) paru = dwarf/dmin1(delta,p1)
!
!   if the input par lies outside of the interval (parl,paru),
!   set par to the closer endpoint.
!
    par = dmax1(par,parl)
    par = dmin1(par,paru)
    if (par .eq. zero) par = gnorm/dxnorm
!
!   beginning of an iteration.
!
150 continue
    iter = iter + 1
!
!   evaluate the function at the current value of par.
!
    if (par .eq. zero) par = dmax1(dwarf,p001*paru)
    temp = dsqrt(par)
    do 160 j = 1, n
       wa1(j) = temp*diag(j)
160 continue
    call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
    do 170 j = 1, n
       wa2(j) = diag(j)*x(j)
170 continue
    dxnorm = enorm(n,wa2)
    temp = fp
    fp = dxnorm - delta
!
!   if the function is small enough, accept the current value
!   of par. also test for the exceptional cases where parl
!   is zero or the number of iterations has reached 10.
!
    if (dabs(fp) .le. p1*delta                       &
        .or. parl .eq. zero .and. fp .le. temp       &
          .and. temp .lt. zero .or. iter .eq. 10) go to 220
!
!   compute the newton correction.
!
    do 180 j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
180 continue
    do 210 j = 1, n
       wa1(j) = wa1(j)/sdiag(j)
       temp = wa1(j)
       jp1 = j + 1
       if (n .lt. jp1) go to 200
       do 190 i = jp1, n
          wa1(i) = wa1(i) - r(i,j)*temp
190    continue
200 continue
210 continue
    temp = enorm(n,wa1)
    parc = ((fp/delta)/temp)/temp
!
!   depending on the sign of the function, update parl or paru.
!
    if (fp .gt. zero) parl = dmax1(parl,par)
    if (fp .lt. zero) paru = dmin1(paru,par)
!
!   compute an improved estimate for par.
!
    par = dmax1(parl,par+parc)
!
!   end of an iteration.
!
    go to 150
220 continue
!
!   termination.
!
    if (iter .eq. 0) par = zero
    return
!
!   last card of subroutine lmpar.
!
    end subroutine lmpar

!------------------------------------------------------------------------------------
!
!   SUBROUTINE QRSOLV
!
!------------------------------------------------------------------------------------

    subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

    integer n,ldr
    integer ipvt(n)
    double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)

    integer i,j,jp1,k,kp1,l,nsing
    double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
    data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
!
!   copy r and (q transpose)*b to preserve input and initialize s.
!   in particular, save the diagonal elements of r in x.
!
    do 20 j = 1, n
      do 10 i = j, n
        r(i,j) = r(j,i)
10    continue
      x(j) = r(j,j)
      wa(j) = qtb(j)
20  continue
!
!   eliminate the diagonal matrix d using a givens rotation.
!
    do 100 j = 1, n
!
!   prepare the row of d to be eliminated, locating the
!   diagonal element using p from the qr factorization.
!
      l = ipvt(j)
      if (diag(l) .eq. zero) go to 90
      do 30 k = j, n
        sdiag(k) = zero
30    continue
      sdiag(j) = diag(l)
!
!     the transformations to eliminate the row of d
!     modify only a single element of (q transpose)*b
!     beyond the first n, which is initially zero.
!
      qtbpj = zero
!
      do 80 k = j, n
!
!       determine a givens rotation which eliminates the
!       appropriate element in the current row of d.
!
        if (sdiag(k) .eq. zero) go to 70
        if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40
        cotan = r(k,k)/sdiag(k)
        sin = p5/dsqrt(p25+p25*cotan**2)
        cos = sin*cotan
        go to 50
40      continue
        tan = sdiag(k)/r(k,k)
        cos = p5/dsqrt(p25+p25*tan**2)
        sin = cos*tan
50      continue
!
!       compute the modified diagonal element of r and
!       the modified element of ((q transpose)*b,0).
!
        r(k,k) = cos*r(k,k) + sin*sdiag(k)
        temp = cos*wa(k) + sin*qtbpj
        qtbpj = -sin*wa(k) + cos*qtbpj
        wa(k) = temp
!
!       accumulate the tranformation in the row of s.
!
        kp1 = k + 1
        if (n .lt. kp1) go to 70
        do 60 i = kp1, n
          temp = cos*r(i,k) + sin*sdiag(i)
          sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
          r(i,k) = temp
60      continue
70     continue
80    continue
90    continue
!
!     store the diagonal element of s and restore
!     the corresponding diagonal element of r.
!
      sdiag(j) = r(j,j)
      r(j,j) = x(j)
100 continue 
!
!   solve the triangular system for z. if the system is
!   singular, then obtain a least squares solution.
!
    nsing = n
    do 110 j = 1, n
      if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
      if (nsing .lt. n) wa(j) = zero
110 continue
    if (nsing .lt. 1) go to 150
    do 140 k = 1, nsing
      j = nsing - k + 1
      sum = zero
      jp1 = j + 1
      if (nsing .lt. jp1) go to 130
      do 120 i = jp1, nsing
        sum = sum + r(i,j)*wa(i)
120   continue
130   continue
      wa(j) = (wa(j) - sum)/sdiag(j)
140 continue
150 continue
!
!   permute the components of z back to components of x.
!
    do 160 j = 1, n
      l = ipvt(j)
      x(l) = wa(j)
160 continue
    return
!
!   last card of subroutine qrsolv.
!
    end subroutine qrsolv

!------------------------------------------------------------------------------------
!
!   SUBROUTINE QRFAC
!
!------------------------------------------------------------------------------------

    subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

    integer m,n,lda,lipvt
    integer ipvt(lipvt)
    logical pivot
    double precision a(lda,n),rdiag(n),acnorm(n),wa(n)

    integer i,j,jp1,k,kmax,minmn
    double precision ajnorm,epsmch,one,p05,sum,temp,zero
!   double precision dpmpar,enorm
    data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!   epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
!   compute the initial column norms and initialize several arrays.
!
    do 10 j = 1, n
      acnorm(j) = enorm(m,a(1,j))
      rdiag(j) = acnorm(j)
      wa(j) = rdiag(j)
      if (pivot) ipvt(j) = j
10  continue
!
!   reduce a to r with householder transformations.
!
    minmn = min0(m,n)
    do 110 j = 1, minmn
      if (.not.pivot) go to 40
!
!     bring the column of largest norm into the pivot position.
!
      kmax = j
      do 20 k = j, n
        if (rdiag(k) .gt. rdiag(kmax)) kmax = k
20    continue
      if (kmax .eq. j) go to 40
      do 30 i = 1, m
        temp = a(i,j)
        a(i,j) = a(i,kmax)
        a(i,kmax) = temp
30    continue
      rdiag(kmax) = rdiag(j)
      wa(kmax) = wa(j)
      k = ipvt(j)
      ipvt(j) = ipvt(kmax)
      ipvt(kmax) = k
40    continue
!
!     compute the householder transformation to reduce the
!     j-th column of a to a multiple of the j-th unit vector.
!
      ajnorm = enorm(m-j+1,a(j,j))
      if (ajnorm .eq. zero) go to 100
      if (a(j,j) .lt. zero) ajnorm = -ajnorm
      do 50 i = j, m
        a(i,j) = a(i,j)/ajnorm
50    continue
      a(j,j) = a(j,j) + one
!
!     apply the transformation to the remaining columns
!     and update the norms.
!
      jp1 = j + 1
      if (n .lt. jp1) go to 100
      do 90 k = jp1, n
        sum = zero
        do 60 i = j, m
          sum = sum + a(i,j)*a(i,k)
60      continue
        temp = sum/a(j,j)
        do 70 i = j, m
          a(i,k) = a(i,k) - temp*a(i,j)
70      continue
        if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
        temp = a(j,k)/rdiag(k)
        rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
        if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
        rdiag(k) = enorm(m-j,a(jp1,k))
        wa(k) = rdiag(k)
80      continue
90    continue
100   continue
      rdiag(j) = -ajnorm
110 continue

    return

    end subroutine qrfac

!------------------------------------------------------------------------------------
!
!   SUBROUTINE FDJAC2
!
!------------------------------------------------------------------------------------

    subroutine fdjac2(fcn,m,ux,uy,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

    integer m,n,ldfjac,iflag
    double precision epsfcn
    double precision ux(m),uy(m)
    double precision x(n),fvec(m),fjac(ldfjac,n),wa(m)

    integer i,j
    double precision eps,epsmch,h,temp,zero
!   double precision dpmpar
    data zero /0.0d0/
!
!   epsmch is the machine precision.
!
    epsmch = dpmpar(1)
!
    eps = dsqrt(dmax1(epsfcn,epsmch))
    do 20 j = 1, n
      temp = x(j)
      h = eps*dabs(temp)
      if (h .eq. zero) h = eps
      x(j) = temp + h
      call fcn(m,ux,uy,n,x,wa,iflag)
      if (iflag .lt. 0) go to 30
      x(j) = temp
      do 10 i = 1, m
        fjac(i,j) = (wa(i) - fvec(i))/h
10    continue
20  continue
30  continue

    return
!
!   last card of subroutine fdjac2.
!
    end subroutine fdjac2

!------------------------------------------------------------------------------------
!
!   FUNCTION ENORM
!
!------------------------------------------------------------------------------------

    double precision function enorm(n,x)

    integer n
    double precision x(n)

    integer i
    double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs
    double precision x1max,x3max,zero

    data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = n
    agiant = rgiant/floatn
    do 90 i = 1, n
      xabs = dabs(x(i))
      if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
      if (xabs .le. rdwarf) go to 30
!
!     sum for large components.
!
      if (xabs .le. x1max) go to 10
      s1 = one + s1*(x1max/xabs)**2
      x1max = xabs
      go to 20
10    continue
      s1 = s1 + (xabs/x1max)**2
20    continue
      go to 60
30    continue
!
!     sum for small components.
!
      if (xabs .le. x3max) go to 40
      s3 = one + s3*(x3max/xabs)**2
      x3max = xabs
      go to 50
40    continue
      if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
50    continue
60    continue
      go to 80
70    continue
!
!     sum for intermediate components.
!
      s2 = s2 + xabs**2
80    continue
90  continue
!
!   calculation of norm.
!
    if (s1 .eq. zero) go to 100
    enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
    go to 130
100 continue
    if (s2 .eq. zero) go to 110
    if (s2 .ge. x3max) enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
    if (s2 .lt. x3max) enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
    go to 120
110 continue
    enorm = x3max*dsqrt(s3)
120 continue
130 continue
!
    return
!
!   last card of function enorm.
!
    end function enorm

!------------------------------------------------------------------------------------
!
!   FUNCTION DPMPAR
!
!------------------------------------------------------------------------------------

    double precision function dpmpar(i)

    integer, intent(in) :: i

    integer mcheps(4)
    integer minmag(4)
    integer maxmag(4)
    double precision dmach(3)
!   equivalence (dmach(1),mcheps(1))
!   equivalence (dmach(2),minmag(1))
!   equivalence (dmach(3),maxmag(1))
!
!   Machine constants for IEEE machines.
! 
    data dmach(1) /2.22044604926d-16/
    data dmach(2) /2.22507385852d-308/
    data dmach(3) /1.79769313485d+308/
!
    dpmpar = dmach(i)
    return
!
!   Last card of function dpmpar.
!
    end function dpmpar

    end module nonfit
