MODULE regrid


  CONTAINS

!=======================================================================

SUBROUTINE regrid1d(nx,x,f,mx,xx,ord_intp,v_fill,g)

  implicit none

  integer,                intent(in)  ::  nx, mx, ord_intp
  real,                   intent(in)  ::  v_fill
  real,    dimension(nx), intent(in)  ::  x
  real,    dimension(nx), intent(in)  ::  f
  real,    dimension(mx), intent(in)  ::  xx
  real,    dimension(mx), intent(out) ::  g

  integer                ::  i, i1, i2, mmx, lw, liw, ier
  logical                ::  l_incx, l_incxx
  real,    dimension(nx) ::  x0
  real,    dimension(nx) ::  f0
  real,    dimension(mx) ::  xx0

  integer, dimension(:), allocatable  ::  iw
  real,    dimension(:), allocatable  ::  w
  real,    dimension(:), allocatable  ::  g0


  ! increasing order
  if ( (x(nx)-x(1)) < 0 ) then
    l_incx = .False.
    do i=1, nx
      x0(i) = x(nx+1-i)
    enddo
  else
    l_incx = .True.
    x0(:) = x(:)
  end if
  if ( (xx(mx)-xx(1)) < 0 ) then
    l_incxx = .False.
    do i=1, mx
      xx0(i) = xx(mx+1-i)
    enddo
  else
    l_incxx = .True.
    xx0(:) = xx(:) 
  end if
  if ( l_incx ) then
    f0(:) = f(:)
  else
    do i=1, nx
      f0(i) = f(nx+1-i)
    enddo
  end if

  call setbdy_regrid(nx,x0,mx,xx0,i1,i2)

  mmx = i2-i1+1

  g(:) = v_fill

  if (mmx > 0) then

    allocate( g0(mmx) )

    lw = mmx+1
    if (ord_intp == 3)  lw = 4*mmx+1
    liw = mmx+1
    allocate( w (lw ) )
    allocate( iw(liw) )

    call rgrd1(nx,x0,f0,mmx,xx0(i1:i2),g0,ord_intp,w,lw,iw,liw,ier)

    if (ier /= 0)  print*, ' ERROR in arguments of rgrd1 : ', ier

    if ( l_incxx ) then
      g(i1:i2) = g0(:)
    else
      do i=i1, i2
        g(i) = g0(mmx+i1-i)
      enddo
    end if

    deallocate( iw, w )
    deallocate( g0 )

  end if

  RETURN

END subroutine regrid1d

!=======================================================================

SUBROUTINE regrid2d(nx,ny,x,y,f,mx,my,xx,yy,ord_intp,v_fill,g)

  implicit none

  integer,                   intent(in)  ::  nx, ny, mx, my, ord_intp(2)
  real,                      intent(in)  ::  v_fill
  real,    dimension(nx),    intent(in)  ::  x
  real,    dimension(ny),    intent(in)  ::  y
  real,    dimension(nx,ny), intent(in)  ::  f
  real,    dimension(mx),    intent(in)  ::  xx
  real,    dimension(my),    intent(in)  ::  yy
  real,    dimension(mx,my), intent(out) ::  g

  integer                   ::  i,j, i1, i2, j1, j2, mmx, mmy, lw, liw, ier
  logical                   ::  l_incx, l_incxx, l_incy, l_incyy
  real,    dimension(nx)    ::  x0
  real,    dimension(ny)    ::  y0
  real,    dimension(nx,ny) ::  f0
  real,    dimension(mx)    ::  xx0
  real,    dimension(my)    ::  yy0

  integer, dimension(:),   allocatable  ::  iw
  real,    dimension(:),   allocatable  ::  w
  real,    dimension(:,:), allocatable  ::  g0


  ! increasing order
  if ( (x(nx)-x(1)) < 0 ) then
    l_incx = .False.
    do i=1, nx
      x0(i) = x(nx+1-i)
    enddo
  else
    l_incx = .True.
    x0(:) = x(:)
  end if
  if ( (xx(mx)-xx(1)) < 0 ) then
    l_incxx = .False.
    do i=1, mx
      xx0(i) = xx(mx+1-i)
    enddo
  else
    l_incxx = .True.
    xx0(:) = xx(:)
  end if
  if ( (y(ny)-y(1)) < 0 ) then
    l_incy = .False.
    do j=1, ny
      y0(j) = y(ny+1-j)
    enddo
  else
    l_incy = .True.
    y0(:) = y(:)
  end if
  if ( (yy(my)-yy(1)) < 0 ) then
    l_incyy = .False.
    do j=1, my
      yy0(j) = yy(my+1-j)
    enddo
  else
    l_incyy = .True.
    yy0(:) = yy(:)
  end if
  if ( l_incx .and. l_incy ) then
    f0(:,:) = f(:,:)
  else if ( l_incx ) then
    do j=1, ny
      f0(:,j) = f(:,ny+1-j)
    enddo
  else if ( l_incy ) then
    do i=1, nx
      f0(i,:) = f(nx+1-i,:)
    enddo
  else
    do j=1, ny
    do i=1, nx
      f0(i,j) = f(nx+1-i,ny+1-j)
    enddo
    enddo
  end if

  call setbdy_regrid(nx,x0,mx,xx0,i1,i2)
  call setbdy_regrid(ny,y0,my,yy0,j1,j2)

  mmx = i2-i1+1
  mmy = j2-j1+1

  g(:,:) = v_fill

  if ( mmx > 0 .and. mmy > 0 ) then

    allocate( g0(mmx,mmy) )

    if (ord_intp(1) == 1) then
      lw = mmx
    else
      lw = 4*mmx+1
    end if
    if (ord_intp(2) == 1) then
      lw = lw + mmy+2*mmx
    else
      lw = lw + 4*(mmx+mmy)
    end if
    lw = lw + 1
    liw = mmx+mmy+1
    allocate( w (lw ) )
    allocate( iw(liw) )

    call rgrd2(nx,ny,x0,y0,f0,mmx,mmy,xx0(i1:i2),yy0(j1:j2),g0,ord_intp, &
               w,lw,iw,liw,ier)

    if (ier /= 0)  print*, ' ERROR in arguments of rgrd2 : ', ier

    if ( l_incxx .and. l_incyy ) then
      g(i1:i2,j1:j2) = g0(:,:)
    else if ( l_incxx ) then
      do j=j1, j2
        g(:,j) = g0(:,mmy+j1-j)
      enddo
    else if ( l_incyy ) then
      do i=i1, i2
        g(i,:) = g0(mmx+i1-i,:)
      enddo
    else
      do j=j1, j2
      do i=i1, i2
        g(i,j) = g0(mmx+i1-i,mmy+j1-j)
      enddo
      enddo
    end if

    deallocate( iw, w )
    deallocate( g0 )

  end if

  RETURN

END subroutine regrid2d

!=======================================================================

SUBROUTINE setbdy_regrid(nx,x,mx,xx,i1,i2)

  implicit none

  integer,                intent(in)    ::  nx, mx
  real,    dimension(mx), intent(in)    ::  xx
  real,    dimension(nx), intent(inout) ::  x
  integer,                intent(out)   ::  i1, i2

  real, parameter ::  v_small = 1.e-3

  integer ::  i
  real    ::  tmp


  i1 = 1  ;  i2 = mx

  ! boundary setting
  tmp = xx(1) - x(1)
  if (tmp < 0) then
    if (abs(tmp) < v_small*abs(x(2)-x(1))) then
      x(1) = xx(1)
    else
      do i=1, mx
        if (xx(i) >= x(1)) then
          i1 = i
          EXIT
        end if
      enddo
    end if
  end if
  tmp = xx(mx) - x(nx)
  if (tmp > 0) then
    if (abs(tmp) < v_small*abs(x(nx)-x(nx-1))) then
      x(nx) = xx(mx)
    else
      do i=mx, 1, -1
        if (xx(i) <= x(nx)) then
          i2 = i
          EXIT
        end if
      enddo
    end if
  end if

  RETURN

END subroutine setbdy_regrid

!=======================================================================

!-----------------------------------------------------------------------
!   RGRD1 and RGRD2 subroutines in REGRID PACK
!   downloaded in http://www.scd.ucar.edu/softlib/mathlib.html
!-----------------------------------------------------------------------                

    subroutine rgrd1(nx,x,p,mx,xx,q,intpol,w,lw,iw,liw,ier)

      implicit none

      real x(*),p(*),xx(*),q(*),w(*)
      integer iw(*)
      integer nx,mx,ier,intpol,lw,liw,i,ii,i1,i2,i3,i4
!
!     check arguments for errors
!
      ier = 1
!
!     check xx grid resolution
!
      if (mx .lt. 1) return
!
!     check intpol
!
      ier = 6
      if (intpol.ne.1 .and. intpol.ne.3) return
!
!     check x grid resolution
!
      ier = 2
      if (intpol.eq.1 .and. nx.lt.2) return
      if (intpol.eq.3 .and. nx.lt.4) return
!
!     check xx grid contained in x grid
!
      ier = 3
      if (xx(1).lt.x(1) .or. xx(mx).gt.x(nx)) return
!
!     check montonicity of grids
!
      do i=2,nx
      if (x(i-1).ge.x(i)) then
        ier = 4
        return
      end if
      end do
      do ii=2,mx
      if (xx(ii-1).gt.xx(ii)) then
        ier = 4
        return
      end if
      end do
!
!     check minimum work space lengths
!
      ier = 5
      if (intpol.eq.1) then
      if (lw .lt. mx) return
      else
      if (lw .lt. 4*mx) return
      end if
      if (liw .lt. mx) return
!
!     arguments o.k.
!
      ier = 0

      if (intpol.eq.1) then
!
!     linear interpolation in x
!
      call linmx(nx,x,mx,xx,iw,w)
      call lint1(nx,p,mx,q,iw,w)
      return
      else
!
!     cubic interpolation in x
!
      i1 = 1
      i2 = i1+mx
      i3 = i2+mx
      i4 = i3+mx
      call cubnmx(nx,x,mx,xx,iw,w(i1),w(i2),w(i3),w(i4))
      call cubt1(nx,p,mx,q,iw,w(i1),w(i2),w(i3),w(i4))
      return
      end if
    end subroutine rgrd1

    subroutine lint1(nx,p,mx,q,ix,dx)
      implicit none
      integer mx,ix(mx),nx,ii,i
      real p(nx),q(mx),dx(mx)
!
!     linearly interpolate p on x onto q on xx
!
      do ii=1,mx
      i = ix(ii)
      q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
      end do
      return
    end subroutine lint1

    subroutine cubt1(nx,p,mx,q,ix,dxm,dx,dxp,dxpp)
      implicit none
      integer mx,ix(mx),nx,i,ii
      real p(nx),q(mx),dxm(mx),dx(mx),dxp(mx),dxpp(mx)
!
!     cubically interpolate p on x to q on xx
!
      do ii=1,mx
      i = ix(ii)
      q(ii) = dxm(ii)*p(i-1)+dx(ii)*p(i)+dxp(ii)*p(i+1)+dxpp(ii)*p(i+2)
      end do
      return
    end subroutine cubt1

    subroutine linmx(nx,x,mx,xx,ix,dx)
!
!     Let x grid pointers for xx grid and interpolation scale terms
!
      implicit none
      real x(*),xx(*),dx(*)
      integer ix(*),isrt,ii,i,nx,mx
      isrt = 1
      do ii=1,mx
!
!     find x(i) s.t. x(i) < xx(ii) <= x(i+1)
!
      do i=isrt,nx-1
        if (x(i+1) .ge. xx(ii)) then
          isrt = i
          ix(ii) = i
          go to 3
        end if
      end do
    3   continue
      end do
!
!     set linear scale term
!
      do ii=1,mx
      i = ix(ii)
      dx(ii) = (xx(ii)-x(i))/(x(i+1)-x(i))
      end do
      return
    end subroutine linmx

    subroutine cubnmx(nx,x,mx,xx,ix,dxm,dx,dxp,dxpp)
      implicit none
      real x(*),xx(*),dxm(*),dx(*),dxp(*),dxpp(*)
      integer ix(*),mx,nx,i,ii,isrt

      isrt = 1
      do ii=1,mx
!
!     set i in [2,nx-2] closest s.t.
!     x(i-1),x(i),x(i+1),x(i+2) can interpolate xx(ii)
!
      do i=isrt,nx-1
        if (x(i+1) .ge. xx(ii)) then
          ix(ii) = min0(nx-2,max0(2,i))
          isrt = ix(ii)
          go to 3
        end if
      end do
    3   continue
      end do
!
!     set cubic scale terms
!
      do ii=1,mx
      i = ix(ii)
      dxm(ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/  &
                ((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
      dx(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/  &
                ((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
      dxp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2))/  &
                ((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
      dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/  &
                ((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
      end do
      return
    end subroutine cubnmx

    subroutine rgrd2(nx,ny,x,y,p,mx,my,xx,yy,q,intpol,w,lw,iw,liw,ier)
      implicit none
      integer nx,ny,mx,my,lw,liw,ier
      integer intpol(2),iw(liw)
      real x(nx),y(ny),p(nx,ny),xx(mx),yy(my),q(mx,my),w(lw)
      integer i,ii,j,jj,j2,j3,j4,j5,j6,j7,j8,j9,i2,i3,i4,i5
      integer jy,lwx,lwy
!
!     check input arguments
!
      ier = 1
!
!     check (xx,yy) grid resolution
!
      if (min0(mx,my) .lt. 1) return
!
!     check intpol
!
      ier = 6
      if (intpol(1).ne.1 .and. intpol(1).ne.3) return
      if (intpol(2).ne.1 .and. intpol(2).ne.3) return
!
!     check (x,y) grid resolution
!
      ier = 2
      if (intpol(1).eq.1 .and. nx.lt.2) return
      if (intpol(1).eq.3 .and. nx.lt.4) return
      if (intpol(2).eq.1 .and. ny.lt.2) return
      if (intpol(2).eq.3 .and. ny.lt.4) return
!
!     check work space lengths
!
      ier = 5
      if (intpol(1).eq.1) then
      lwx = mx
      else
      lwx = 4*mx
      end if
      if (intpol(2).eq.1) then
      lwy = my+2*mx
      else
      lwy = 4*(mx+my)
      end if
      if (lw .lt. lwx+lwy) return
      if (liw .lt. mx+my) return
!
!     check (xx,yy) grid contained in (x,y) grid
!
      ier = 3
      if (xx(1).lt.x(1) .or. xx(mx).gt.x(nx)) return
      if (yy(1).lt.y(1) .or. yy(my).gt.y(ny)) return
!
!     check montonicity of grids
!
      ier = 4
      do i=2,nx
      if (x(i-1).ge.x(i)) return
      end do
      do j=2,ny
      if (y(j-1).ge.y(j)) return
      end do
      do ii=2,mx
      if (xx(ii-1).gt.xx(ii)) return
      end do
      do jj=2,my
      if (yy(jj-1).gt.yy(jj)) return
      end do
!
!     arguments o.k.
!
      ier = 0
!
!     set pointer in integer work space
!
      jy = mx+1
      if (intpol(2) .eq.1) then
!
!     linearly interpolate in y
!
      j2 = 1
      j3 = j2
      j4 = j3+my
      j5 = j4
      j6 = j5
      j7 = j6
      j8 = j7+mx
      j9 = j8+mx
!
!     set y interpolation indices and scales and linearly interpolate
!
      call linmx(ny,y,my,yy,iw(jy),w(j3))
      i2 = j9
!
!     set work space portion and indices which depend on x interpolation
!
      if (intpol(1) .eq. 1) then
      i3 = i2
      i4 = i3
      i5 = i4
      call linmx(nx,x,mx,xx,iw,w(i3))
      else
      i3 = i2+mx
      i4 = i3+mx
      i5 = i4+mx
      call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
      end if
      call lint2(nx,ny,p,mx,my,q,intpol,iw(jy),w(j3),   &
                  w(j7),w(j8),iw,w(i2),w(i3),w(i4),w(i5))
      return

      else
!
!     cubically interpolate in y, set indice pointers
!
      j2 = 1
      j3 = j2+my
      j4 = j3+my
      j5 = j4+my
      j6 = j5+my
      j7 = j6+mx
      j8 = j7+mx
      j9 = j8+mx
      call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
      i2 =  j9+mx
!
!     set work space portion and indices which depend on x interpolation
!
      if (intpol(1) .eq. 1) then
      i3 = i2
      i4 = i3
      i5 = i4
      call linmx(nx,x,mx,xx,iw,w(i3))
      else
      i3 = i2+mx
      i4 = i3+mx
      i5 = i4+mx
      call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
      end if
      call cubt2(nx,ny,p,mx,my,q,intpol,iw(jy),w(j2),w(j3),   &
       w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),iw,w(i2),w(i3),w(i4),w(i5))
      return
      end if
    end subroutine rgrd2

    subroutine lint2(nx,ny,p,mx,my,q,intpol,jy,dy,pj,pjp,   &
                       ix,dxm,dx,dxp,dxpp)
      implicit none
      integer nx,ny,mx,my,intpol(2),jy(my),ix(mx)
      integer jsave,j,jj,ii
      real p(nx,ny),q(mx,my)
      real pj(mx),pjp(mx),dy(my)
      real dxm(mx),dx(mx),dxp(mx),dxpp(mx)
!
!     linearly interpolate in y
!
      if (intpol(1).eq.1) then
!
!     linear in x
!
      jsave = -1
      do jj=1,my
      j = jy(jj)
      if (j.eq.jsave) then
!
!       j pointer has not moved since last pass (no updates or interpolation)
!
      else if (j.eq.jsave+1) then
!
!       update j and interpolate j+1
!
        do ii=1,mx
          pj(ii) = pjp(ii)
        end do
        call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
      else
!
!       interpolate j,j+1in pj,pjp on xx mesh
!
      call lint1(nx,p(1,j),mx,pj,ix,dx)
      call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
      end if
!
!       save j pointer for next pass
!
      jsave = j
!
!       linearly interpolate q(ii,jj) from pjp,pj in y direction
!
      do ii=1,mx
        q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
      end do
      end do

      else
!
!     cubic in x
!
      jsave = -1
      do jj=1,my
      j = jy(jj)
      if (j.eq.jsave) then
!
!       j pointer has not moved since last pass (no updates or interpolation)
!
      else if (j.eq.jsave+1) then
!
!       update j and interpolate j+1
!
        do ii=1,mx
          pj(ii) = pjp(ii)
        end do
        call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
      else
!
!       interpolate j,j+1 in pj,pjp on xx mesh
!
        call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
        call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
      end if
!
!       save j pointer for next pass
!
      jsave = j
!
!       linearly interpolate q(ii,jj) from pjp,pj in y direction
!
      do ii=1,mx
        q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
      end do
      end do
      return
      end if
    end subroutine lint2

    subroutine cubt2(nx,ny,p,mx,my,q,intpol,jy,dym,dy,dyp,  &
      dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
      implicit none
      integer nx,ny,mx,my,intpol(2),jy(my),ix(mx)
      integer jsave,j,jj,ii
      real p(nx,ny),q(mx,my)
      real pjm(mx),pj(mx),pjp(mx),pjpp(mx)
      real dym(my),dy(my),dyp(my),dypp(my)
      real dxm(mx),dx(mx),dxp(mx),dxpp(mx)
      if (intpol(1).eq.1) then
!
!     linear in x
!
      jsave = -3
      do jj=1,my
!
!       load closest four j lines containing interpolate on xx mesh
!       for j-1,j,j+1,j+2 in pjm,pj,pjp,pjpp
!
      j = jy(jj)
      if (j.eq.jsave) then
!
!       j pointer has not moved since last pass (no updates or interpolation)
!
      else if (j.eq.jsave+1) then
!
!       update j-1,j,j+1 and interpolate j+2
!
        do ii=1,mx
          pjm(ii) = pj(ii)
          pj(ii) = pjp(ii)
          pjp(ii) = pjpp(ii)
        end do
        call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
      else if (j.eq.jsave+2) then
!
!     update j-1,j and interpolate j+1,j+2
!
        do ii=1,mx
          pjm(ii) = pjp(ii)
          pj(ii) = pjpp(ii)
        end do
        call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
        call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
      else if (j.eq.jsave+3) then
!
!       update j-1 and interpolate j,j+1,j+2
!
        do ii=1,mx
          pjm(ii) = pjpp(ii)
        end do
        call lint1(nx,p(1,j),mx,pj,ix,dx)
        call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
        call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
      else
!
!       interpolate all four j-1,j,j+1,j+2
!
        call lint1(nx,p(1,j-1),mx,pjm,ix,dx)
        call lint1(nx,p(1,j),mx,pj,ix,dx)
        call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
        call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
      end if
!
!     save j pointer for next pass
!
      jsave = j
!
!     cubically interpolate q(ii,jj) from pjm,pj,pjp,pjpp in y direction
!
      do ii=1,mx
        q(ii,jj) = dym(jj)*pjm(ii)+dy(jj)*pj(ii)+dyp(jj)*pjp(ii)+   &
                     dypp(jj)*pjpp(ii)
      end do
      end do
      return

      else
!
!     cubic in x
!
      jsave = -3
      do jj=1,my
!
!       load closest four j lines containing interpolate on xx mesh
!       for j-1,j,j+1,j+2 in pjm,pj,pjp,pjpp
!
        j = jy(jj)
        if (j.eq.jsave) then
!
!         j pointer has not moved since last pass (no updates or interpolation)
!
        else if (j.eq.jsave+1) then
!
!         update j-1,j,j+1 and interpolate j+2
!
          do ii=1,mx
            pjm(ii) = pj(ii)
            pj(ii) = pjp(ii)
            pjp(ii) = pjpp(ii)
          end do
          call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
        else if (j.eq.jsave+2) then
!
!         update j-1,j and interpolate j+1,j+2
!
          do ii=1,mx
            pjm(ii) = pjp(ii)
            pj(ii) = pjpp(ii)
          end do
          call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
        else if (j.eq.jsave+3) then
!
!         update j-1 and interpolate j,j+1,j+2
!
          do ii=1,mx
            pjm(ii) = pjpp(ii)
          end do
          call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
        else
!
!         interpolate all four j-1,j,j+1,j+2
!
          call cubt1(nx,p(1,j-1),mx,pjm,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
          call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
        end if
!
!       save j pointer for next pass
!
      jsave = j
!
!       cubically interpolate q(ii,jj) from pjm,pj,pjp,pjpp in y direction
!
      do ii=1,mx
        q(ii,jj) = dym(jj)*pjm(ii)+dy(jj)*pj(ii)+dyp(jj)*pjp(ii)+  &
                     dypp(jj)*pjpp(ii)
      end do
      end do
      return
      end if
    end subroutine cubt2

!=======================================================================

END module regrid
