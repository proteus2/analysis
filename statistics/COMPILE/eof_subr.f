      subroutine pcf(tt,mx,w,z,lxend,lxch,nch,std,avg,ccov,
     +               datamiss,lout,lout2,pc,rc,ts)
      dimension tt(mx,nch),ts(mx,nch),rc(mx,lout,nch),pc(mx,lout)
      dimension w(lxch),z(lxch,lxch),std(nch),avg(nch)
      character*3 ccov
      print *,'@@ welcome to subroutine pcf   @@'

      do 20 k=1,nch
      do 20 i= 1,mx
      do 20 j= 1,lout
      rc(i,j,k)=0.0
      pc(i,j)=0.0
   20 continue
      lxchd2=lxch/2
      do 50 i=1,lxchd2
      ii=lxch-i+1
      temp=w(i)
      w(i)=w(ii)
      w(ii)=temp
   50 continue
      do 30 i=1,lxch
      do 30 j=1,lxchd2
      jj=lxch-j+1
      temp=z(i,j)
      z(i,j)=z(i,jj)
      z(i,jj)=temp
   30 continue
      do 60 i=1,mx-lxend+1
      do 60 k=1,lout
      nc=0
      do 65 l=1,nch
      nn=(l-1)*lxend
      do 65 m=1,lxend
      if (tt(i+m-1,l).ne.datamiss) then
      nc=nc+1
      pc(i,k)=pc(i,k)+tt(i+m-1,l)*z(nn+m,k)
      endif
   65 continue
      if (nc.eq.0) pc(i,k)=datamiss
   60 continue
      do 70 l=1,nch
      nn=(l-1)*lxend
      do 75 i=lxend,mx-lxend+1
      do 75 k=1,lout
      nc=0
      do 751 m=1,lxend
      if (pc(i-m+1,k).ne.datamiss) then
      nc=nc+1
      rc(i,k,l)=rc(i,k,l)+pc(i-m+1,k)*z(nn+m,k)
      endif
  751 continue
      if (nc.eq.0) then
      rc(i,k,l)=datamiss
      else
        rc(i,k,l)=rc(i,k,l)/float(nc)
      endif
   75 continue
      do 76 i=1,lxend-1
      do 76 k=1,lout
      nc=0
      do 761 m=1,i
      if (pc(i-m+1,k).ne.datamiss) then
      nc=nc+1
      rc(i,k,l)=rc(i,k,l)+pc(i-m+1,k)*z(nn+m,k)
      endif
  761 continue
      if (nc.eq.0) then
      rc(i,k,l)=datamiss
      else
      rc(i,k,l)=rc(i,k,l)/float(nc)
      endif
   76 continue
      do 77 i=mx-lxend+2,mx
      do 77 k=1,lout
      nc=0
      do 771 m=i-mx+lxend,lxend
      if (pc(i-m+1,k).ne.datamiss) then
      nc=nc+1
      rc(i,k,l)=rc(i,k,l)+pc(i-m+1,k)*z(nn+m,k)
      endif
  771 continue
      if (nc.eq.0) then
      rc(i,k,l)=datamiss
      else
      rc(i,k,l)=rc(i,k,l)/float(nc)
      endif
   77 continue
   70 continue

c      do 95 i=1,mx
c      if (i.ge.lxend.and.i.le.mx-lxend+1) de=float(lxend)
c      if (i.ge. 1.and.i.le.lxend-1) de=float(i)
c      if (i.ge.mx-lxend+2.and.i.le.mx) de=float(mx-i+1)
c      do 95 k=1,lout
c      do 95 l=1,nch
c      rc(i,k,l)=rc(i,k,l)/de
c   95 continue

      do 898 iii=1,nch
      do 898 jjj=1,mx
      td=0.0
      ts(jjj,iii)=0.
      nc=0
      do 897 kkk=1,lout
      if (rc(jjj,kkk,iii).ne.datamiss) then
      nc=nc+1
      td=td+rc(jjj,kkk,iii)
      ts(jjj,iii)=ts(jjj,iii)+rc(jjj,kkk,iii)
      endif    
  897 continue
      if (nc.ne.0) then
      if (ccov.eq.'cor') ts(jjj,iii)=ts(jjj,iii)*std(iii)
      ts(jjj,iii)=ts(jjj,iii)+avg(iii)
      endif
      ttss=ts(jjj,iii)
      if (tt(jjj,iii).ne.datamiss) then
      if (ccov.eq.'cor') tt(jjj,iii)=tt(jjj,iii)*std(iii)
      tttt=tt(jjj,iii)+avg(iii)
      else
      tttt=datamiss
      endif
  898 continue
      return
      end

cc***
      subroutine mean(tt,mx,nch,datamiss,ccov,std,avg)
c     mean and variance
      dimension tt(mx,nch),std(nch),avg(nch)
      character*3 ccov
      print *,'@@ welcome to subroutine mean  @@'
      do 200 j=1,nch
      sum=0.
      nc=0
      do  100 i=1,mx
      if (tt(i,j).ne.datamiss) then
      nc=nc+1
      sum=sum+tt(i,j)
      endif
  100 continue
      if (nc.eq.0) print*,j,'  err in mean'
      sum=sum/float(nc)
      if (nc.eq.0) then
      avg(j)=datamiss
      else
      avg(j)=sum
      endif
      var=0.
      nc=0
      do 300 i=1,mx
      if (tt(i,j).ne.datamiss) then
      nc=nc+1
      tt(i,j)=tt(i,j)-sum
      var=var+tt(i,j)*tt(i,j)
      endif
  300 continue
      if (ccov.eq.'cov') goto 200
      var=var/float(nc)
      if (nc.eq.0) then
      std(j)=datamiss
      else
      std(j)=sqrt(var)
      endif
        if (var.le.0.) print*,'error in subroutine mean'
        if (nc.eq.0) print*,j,'  err in mean 2'
      do 400 i=1,mx
      if (tt(i,j).ne.datamiss) tt(i,j)=tt(i,j)/sqrt(var)
  400 continue
  200 continue
      return
      end
      subroutine covf(i1,j1,i2,j2,tt,nmax,cov,lxend,nch,datamiss)
c     covariance function
      dimension tt(nmax,nch)
c     print *,'@@ welcome to subroutine covf @@'
      sum=0.
      nc=0
      do 10 i=1,nmax-lxend+1
      a=tt(i+j1-1,i1)
      b=tt(i+j2-1,i2)
      if (a.ne.datamiss.and.b.ne.datamiss) then
      nc=nc+1
      sum=sum+a*b
      endif
   10 continue
      if (nc.eq.0) print *,'err in subroutine covf'
      cov=sum/float(nc)
      return
      end

cc***
      subroutine symtrx(rowij,m,root,eigv,ni,wk,ier)
c$name  symtrx
c$link  cpcon
c solves eigenfunction problem for symmetric matrix
c--- history
c 89.12. 4 modified with cncpu
c 90. 1. 6 cncpu is replaced by prec.
c--- input
c m       i      order of original symmetric matrix
c ni      i      initial dimension of -root-, -eigv- and -wk-
c rowij  r(*)    symmetric storage mode of order m*(m+1)/2
c--- output
c eigv   r(ni,m) eigenvectors of original symmetric matrix
c ier     i      index for root(j) failed to converge (j=ier-128)
c root   r(m)    eigenvalues of original symmetric matrix
c rowij          storage of householder reduction elements
c wk             work area
c$endi
      real rowij(*),root(*),wk(*),eigv(ni,*),ccp(3)
c+++ add epscp
c     data  rdelp/ 1.1921e-07 /
c
      print *,'@@ welcome to subroutine symtrx @@'
      call cpcon(ccp)
      rdelp=ccp(1)*10
c+++
      ier = 0
      mp1 = m + 1
      mm = (m*mp1)/2 - 1
      mbeg = mm + 1- m
c
c+---------------------------------------------------------------------+
c|          loop-100 reduce -rowij- (symmetric storage mode) to a      |
c|          symmetric tridiagonal form by householder method           |
c|                      cf. wilkinson, j.h., 1968,                     |
c|              the algebraic eigenvalue problem, pp 290-293.          |
c|          loop-30&40 and 50 form element of a*u and element p        |
c+---------------------------------------------------------------------+
      do 100 ii=1,m
      i = mp1 - ii
      l = i - 1
      h = 0.0
      scale = 0.0
      if (l.lt.1) then
c|          scale row (algol tol then not needed)
      wk(i) = 0.0
      go to 90
      end if
      mk = mm
      do 10 k=1,l
      scale = scale + abs(rowij(mk))
      mk = mk - 1
   10 continue
      if (scale.eq.0.0) then
      wk(i) = 0.0
      go to 90
      end if
c
      mk = mm
      do 20 k = 1,l
      rowij(mk) = rowij(mk)/scale
      h = h + rowij(mk)*rowij(mk)
      mk = mk - 1
   20 continue
      wk(i) = scale*scale*h
      f = rowij(mm)
      g = - sign(sqrt(h),f)
      wk(i) = scale*g
      h = h - f*g
      rowij(mm) = f - g
      if (l.gt.1) then
      f = 0.0
      jk1 = 1
      do 50 j=1,l
      g = 0.0
      ik = mbeg + 1
      jk = jk1
      do 30 k=1,j
      g = g + rowij(jk)*rowij(ik)
      jk = jk + 1
      ik = ik + 1
   30 continue
      jp1 = j + 1
      if (l.ge.jp1) then
      jk = jk + j - 1
      do 40 k=jp1,l
      g = g + rowij(jk)*rowij(ik)
      jk = jk + k
      ik = ik + 1
   40 continue
      end if
      wk(j) = g/h
      f = f + wk(j)*rowij(mbeg+j)
      jk1 = jk1 + j
   50 continue
      hh = f/(h+h)
c
      jk = 1
      do 70 j=1,l
      f = rowij(mbeg+j)
      g = wk(j) - hh*f
      wk(j) = g
      do 60 k=1,j
      rowij(jk) = rowij(jk) - f*wk(k) - g*rowij(mbeg+k)
      jk = jk + 1
   60 continue
   70 continue
      end if
c
      do 80 k=1,l
      rowij(mbeg+k) = scale*rowij(mbeg+k)
   80 continue
   90 root(i) = rowij(mbeg+i)
      rowij(mbeg+i) = h*scale*scale
      mbeg = mbeg - i + 1
      mm = mm - i
  100 continue
c
c+---------------------------------------------------------------------+
c|          loop-210 compute eigenvalues and eigenvectors              |
c|          setup work area location eigv to the identity matrix       |
c|          loop-140 for finding small sub-diagonal element            |
c|          loop-160 for convergence of eigenvalue j (max. 30 times)   |
c|          loop-190 for ql transformation and loop-180 form vectors   |
c+---------------------------------------------------------------------+
      do 110 i=1,m-1
  110 wk(i) = wk(i+1)
      wk(m) = 0.0
      b = 0.0
      f = 0.0
      do 130 i=1,m
      do 120 j=1,m
  120 eigv(i,j) = 0.0
      eigv(i,i) = 1.0
  130 continue
c
      do 210 l=1,m
      j = 0
      h = rdelp*(abs(root(l))+abs(wk(l)))
      if (b.lt.h) b = h
      do 140 n=l,m
      k = n
      if (abs(wk(k)).le.b) go to 150
  140 continue
  150 n = k
      if (n.eq.l) go to 200
c
  160 continue
      if (j.eq.30) then
      ier = 128 + l
      return
      end if
c
      j = j + 1
      l1 = l + 1
      g = root(l)
      p = (root(l1)-g)/(wk(l)+wk(l))
      r = abs(p)
      if (rdelp*abs(p).lt.1.0) r = sqrt(p*p+1.0)
      root(l) = wk(l)/(p+sign(r,p))
      h = g - root(l)
      do 170 i=l1,m
      root(i) = root(i) - h
  170 continue
      f = f + h
c
      p = root(n)
      c = 1.0
      s = 0.0
      nn1 = n - 1
      nn1pl = nn1 + l
      if (l.le.nn1) then
      do 190 ii=l,nn1
      i = nn1pl - ii
      g = c*wk(i)
      h = c*p
      if (abs(p).lt.abs(wk(i))) then
      c = p/wk(i)
      r = sqrt(c*c+1.0)
      wk(i+1) = s*wk(i)*r
      s = 1.0/r
      c = c*s
      else
      c = wk(i)/p
      r = sqrt(c*c+1.0)
      wk(i+1) = s*p*r
      s = c/r
      c = 1.0/r
      end if
      p = c*root(i) - s*g
      root(i+1) = h + s*(c*g+s*root(i))
      if (ni.ge.m) then
      do 180 k=1,m
      h = eigv(k,i+1)
      eigv(k,i+1) = s*eigv(k,i) + c*h
      eigv(k,i) = c*eigv(k,i) - s*h
  180 continue
      end if
  190 continue
      end if
      wk(l) = s*p
      root(l) = c*p
      if (abs(wk(l)).gt.b) go to 160
  200 root(l) = root(l) + f
  210 continue
c
c+---------------------------------------------------------------------+
c|          back transform eigenvectors of the original matrix from    |
c|          eigenvectors 1 to m of the symmetric tridiagonal matrix    |
c+---------------------------------------------------------------------+
      do 250 i=2,m
      l = i - 1
      ia = (i*l)/2
      if (rowij(ia+i).ne.0.0) then
      do 240 j=1,m
      sum = 0.0
      do 220 k=1,l
      sum = sum + rowij(ia+k)*eigv(k,j)
  220 continue
      sum = sum/rowij(ia+i)
      do 230 k=1,l
      eigv(k,j) = eigv(k,j) - sum*rowij(ia+k)
  230 continue
  240 continue
      end if
  250 continue
c
      return
      end

cc***
      subroutine cpcon(c)
c$name  cpcon
c machine constants of computer
c--- output
c c     r(3)      (1)  minimum positive x for  1+x      .ne. 1
c                 (2)  minimum exponent y for  10.0**y  .ne. 0
c                 (3)  maximum exponent z for  10.0**z  is max. value
c                  if init=1 (data statement) then  set as z=y
c                  if init=2 then this routine gets actual value of z.
c                  - see note -
c--- history
c 90. 1.20  created
c
c--- note
c this program will generate -underflow error- and -overflow error-
c  messages.  on some computer -overflow error- message may be
c  fatal error.  in that case, please set init = 1 in the data
c  satement for suppressing the procedure of getting c(3).
      dimension c(3)
c resolution of computation
      save init,x,y,z
      data init/1/
      if(init.le.0) then
      c(1)=x
      c(2)=y
      c(3)=z
      return
      endif
      n=500
      x=1
      do 1 i=1,n
      x=x/2
      x1=1+x
      if(abs(x1-1) .le. 0) then
      x=2*x
      goto 2
      endif
    1 continue
    2 c(1)=x
c exponent for minimum positive value
c  this procedure will generate -underflow error massage-
      y2=1
      n=500
      do 3 i=1,n
      y1=y2
      y2=y1/10
      if(abs(10*y2/y1-1) .gt. 20*x) goto 4
    3 continue
      i=n+1
    4 y=1-i
      c(2)=y
c exponent for maximum positive value
c this procedure will generate -overflow message-
      if(init.le.1) then
      z=-y
       else
      z2=1
      n=500
      do 5 i=1,n
      z1=z2
      z2=z1*10
      if(abs(z2/z1/10-1) .gt. 20*x) goto 6
    5   continue
      i=n+1
    6   z=i-1
      endif
      c(3)=z
c
      init=0
      return
      end

