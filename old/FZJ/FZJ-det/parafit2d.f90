MODULE parafit2d

  implicit none

  integer ::  nx_pf2, ny_pf2

  real, dimension(5,5) ::  mat_pf2

  real*8, dimension(:,:), allocatable ::  x1_pf2, y1_pf2, &
                                          x2_pf2, y2_pf2

  contains

SUBROUTINE para_fit_2d(nx,ny,f,a,yfit)

  implicit none

  integer,                intent(in) ::  nx, ny
  real, dimension(nx,ny), intent(in) ::  f

  real, dimension(5),     intent(out) ::  a
  real, dimension(nx,ny), intent(out) ::  yfit

  integer ::  i, j
  logical ::  l_calc

  real                     ::  mat4(5,5)
  real*8                   ::  vec(5)
  real*8, dimension(nx,ny) ::  x1, y1, x2, y2

  l_calc = .TRUE.
  if ( nx == nx_pf2 .and. ny == ny_pf2 )  l_calc = .FALSE.

  if ( l_calc ) then

    do i=1, nx
      x1(i,:) = 2.d0/dble(nx-1)*dble(i-1) - 1.d0
    enddo
    do j=1, ny
      y1(:,j) = 2.d0/dble(ny-1)*dble(j-1) - 1.d0
    enddo
    x2(:,:) = x1(:,:)*x1(:,:)
    y2(:,:) = y1(:,:)*y1(:,:)

    mat_pf2(:,5) = (/sum(x2),sum(x1),sum(y2),sum(y1),dble(nx*ny)/)
    mat_pf2(:,4) = (/sum(x2*y1),sum(x1*y1),sum(y2*y1),sum(y2),sum(y1)/)
    mat_pf2(:,3) = (/sum(x2*y2),sum(x1*y2),sum(y2*y2),sum(y1*y2),sum(y2)/)
    mat_pf2(:,2) = (/sum(x2*x1),sum(x2),sum(y2*x1),sum(y1*x1),sum(x1)/)
    mat_pf2(:,1) = (/sum(x2*x2),sum(x1*x2),sum(y2*x2),sum(y1*x2),sum(x2)/)
    mat_pf2(:,:) = mat_pf2(:,:)/dble(nx*ny)

    mat4 = real(mat_pf2)
    call inv_matrix(5,mat4)
    mat_pf2 = dble(mat4)

    ! for the next time
    nx_pf2 = nx  ;  ny_pf2 = ny
    if ( allocated(x1_pf2) )  &
       deallocate( x1_pf2, y1_pf2, x2_pf2, y2_pf2 )
    allocate( x1_pf2(nx,ny), y1_pf2(nx,ny), &
              x2_pf2(nx,ny), y2_pf2(nx,ny) )
    x1_pf2(:,:) = x1(:,:)  ;  y1_pf2(:,:) = y1(:,:)
    x2_pf2(:,:) = x2(:,:)  ;  y2_pf2(:,:) = y2(:,:)

  else

    x1(:,:) = x1_pf2(:,:)  ;  y1(:,:) = y1_pf2(:,:)
    x2(:,:) = x2_pf2(:,:)  ;  y2(:,:) = y2_pf2(:,:)

  end if

  vec = (/sum(f*x2),sum(f*x1),sum(f*y2),sum(f*y1),sum(f)/)
  vec(:) = vec(:)/dble(nx*ny)

  do j=1, 5
    a(j) = real(sum(mat_pf2(:,j)*vec(:)))
  enddo

  yfit(:,:) = real( a(1)*x2(:,:) + a(2)*x1(:,:) + &
                    a(3)*y2(:,:) + a(4)*y1(:,:) + a(5) )

  RETURN

END subroutine para_fit_2d


SUBROUTINE inv_matrix(n,mat)

  implicit none

  integer, intent(in) ::  n

  real, dimension(n,n), intent(inout) ::  mat

  integer ::  i, indx(n)
  real    ::  d, y(n,n)

  y(:,:) = 0.
  do i=1, n
    y(i,i) = 1.
  enddo
  call ludcmp(mat,n,n,indx,d)
  do i=1, n
    call lubksb(mat,n,n,indx,y(1,i))
  enddo

  mat(:,:) = y(:,:)

  RETURN

END subroutine inv_matrix

END module parafit2d


      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END subroutine

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END subroutine

