      module specanal

      public :: fps, crosss

      contains


!===============================================================
      subroutine fps(a,dx,ldata,m,fr,ps1,whitening)
!---------------------------------------------------------------
!     m         :  maximum lag. usually use ldata/2
!     whitening :  lag-1 autocorrelation for pre-whitening.
!                  0. for no pre-whitening
!     ps1       :  [ (a unit)**2 / (rad/time unit) ]
!     fr        :  [ cycle / time unit ]
!     dx        :  [ time unit ]
!---------------------------------------------------------------

      integer ldata
      integer m, m2, mp1, mp2           ! m : maximum lag
      integer ind(6)
      dimension a(ldata),tp(m+1),ps1(m+1),rn(m+1),ts(ldata)
	real*4 x(ldata*2),xind(2),xy(6), dx, whitening
	real*4 acv(m*2),fr(m+1),ps((m+1)*2),xc(m*2+1)
	real*4 xsp((m+1)*2),aph((m+1)*2),xf((m+1)*2),ch(m+1)

      m2 = m*2
      mp1 = m+1
      mp2 = mp1*2

      ind(1)=0
      ind(2)=ldata
      ind(3)=0
      ind(4)=m
      ind(5)=1
      ind(6)=0 
      xind(1)=dx          ! interval of data
      xind(2)= whitening
      if (whitening .eq. 0.)  ind(5) = 0


      do 100 k=1,ldata
100   x(k)=a(k)

      call ftfreq(x,ind,xind,xy,acv,fr,ps,xc,xsp,aph,xf,ch,ier)

      do 200 i=2,mp1
200   tp(i)=1./fr(i)
      tp(1)=99999.
      do 210 i=1,mp1
210   ps1(i)=ps(i)

!      do 500 i=1,ldata
!500   ts(i)=a(i)
!      call rednos(ts,ldata,ps1,rn,mp1)

!21    format(/,3x,'frequency,   period,     power sp.,    red noise',/)
!22    format(2x,'i=',i4,3f14.3)
!      open(6,file='../pow.res',status='unknown')
!      write(6,21)
!      do 111 i=1,mp1
!111   write(6,22) i,tp(i),ps(i),rn(i)

      return
      end subroutine

!==============================================================
      subroutine crosss(a,b,dx,ldata,m,corr,fr,ps1,ps2,csamp,csphs,ch)

      integer ldata
      integer m, m2, mp1, mp2           ! m : maximum lag
      integer ind(6)
      dimension a(ldata),b(ldata)
      real x(ldata*2),xind(2),xy(6)
      real acv(m*2),fr(m+1),ps((m+1)*2),xc(m*2+1),corr(m*2+1)
      real xsp((m+1)*2),aph((m+1)*2),xf((m+1)*2),ch(m+1)
      real tp1(m+1),ps1(m+1),ps2(m+1),csamp(m+1),csphs(m+1)
      real auto1(m+1),auto2(m+1)

      m2 = m*2
      mp1 = m+1
      mp2 = mp1*2

      ind(1)=0
      ind(2)=ldata
      ind(3)=1
      ind(4)=m
      ind(5)=0
      ind(6)=0
      xind(1)=dx
      xind(2)=1.5

      do 100 k=1,ldata
      x(k)=a(k)
100   x(k+ldata)=b(k)

      call ftfreq(x,ind,xind,xy,acv,fr,ps,xc,xsp,aph,xf,ch,ier)
      ps1(1)=ps(1)
      tp1(1)=99999.
      do 112 i=2,mp1
      ps1(i)=ps(i)
112   tp1(i)=1./fr(i)

      ind(1)=1
      ind(2)=ldata
      ind(3)=1
      ind(4)=m
      ind(5)=0
      ind(6)=0
      xind(1)=dx
      xind(2)=1.5

      call ftfreq(x,ind,xind,xy,acv,fr,ps,xc,xsp,aph,xf,ch,ier)
      do i=1, mp1
        ps2(i) = ps(i+mp1)
        csamp(i) = aph(i)*360.
        csphs(i) = aph(i+mp1)*360.
        if (csphs(i) .gt. 180.)  csphs(i) = csphs(i)-360.
        ch(i) = sqrt(ch(i))
      enddo
      do i=1, m*2+1
        corr(i) = xc(i)/(sqrt(xy(2))*sqrt(xy(4)))
      enddo
      auto1(1)=1.
      auto2(1)=1.
      do 227 i=1,m
        auto1(i+1)=acv(i)/xy(2)
        auto2(i+1)=acv(i+m)/xy(4)
227   continue

!      open(9,file='c-pow',status='unknown')
!21    format(/,11x,'preiod,   power x,  power y,
!     &    amp,     phase,  coherence'/)
!22    format(2x,'i=',i4,8f10.3)
!      write(9,21)
!      do 111 i=1,mp1
!111   write(9,22) i,tp1(i),ps1(i),ps2(i),csamp(i),csphs(i),ch(i)

!      do 220 i=1,m*2+1
!220   write(22,*) corr(i)

!      call dof(auto1,auto2,mp1,ndof)
!      print*, ndof

      return
      end subroutine

!==============================================================
      end module specanal
!==============================================================


      subroutine rednos(ts,l,ps,rn,mp1)

      dimension ts(l),ps(mp1),rn(mp1)
      x=0.
      do 20 i=1,l
20    x=x+ts(i)
      x=x/float(l)
      do 22 i=1,l
22    ts(i)=ts(i)-x
      x=0.
      do 24 i=1,l
24    x=x+ts(i)*ts(i)
      y=ts(1)*ts(l)
      do 26 i=2,l
26    y=y+ts(i-1)*ts(i)
      r1=y/x
      x=0.
      do 28 k=1,mp1
28    x=x+ps(k)
      psb=x/float(mp1)
      upp=psb*(1.-r1*r1)
      de1=     1.+r1*r1
      pm=3.1415926/float(mp1-1)
      do 30 k=1,mp1
      ang=pm*float(k-1)
      de2=2.*r1*cos(ang)
30    rn(k)=upp/(de1-de2)
      return
      end subroutine


      subroutine dof(r1,r2,n,ndof)

      dimension r1(n),r2(n)
      c=0.
      do 10 i=1,n
      c=c+r1(i)*r2(i)
10    continue
      if(c.lt.0) print*,'c less then 0   ',c
      d=1.+2.*c
      ndof=n/d
      return
      end subroutine

