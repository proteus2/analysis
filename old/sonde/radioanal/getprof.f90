PROGRAM get_prof

  use netcdfio
  use regress
  
  implicit none

! data properties --------------------------------------------------------------
  integer, parameter                    ::  maxl     = 100000
  real,    parameter                    ::  MIS_in   = -999.
! filtering properties ---------------------------------------------------------
  real,    parameter                    ::  limit_in = 50.                 ! [m]
  integer, parameter                    ::  itv      = 10
! interpolation ----------------------------------------------------------------
  real,    parameter                    ::  dz       = 50.                 ! [m]
! analysis range ---------------------------------------------------------------
  real,    parameter                    ::  zbst     = 17.e3               ! [m]
  real,    parameter                    ::  ztst     = 30.e3
  real,    parameter                    ::  zbuf     = 50.
  integer, parameter                    ::  order_p  = 3 
!-------------------------------------------------------------------------------
! maxl     : maximum number of lines to read in the input data
! MIS_in   : missing value assigned in the input data
! limit_in : length used to filter the original data in process of data-reading
!            when the height interval is too short (i.e., shorter than limit_in)
!            The filtered data is made so that it has intervals (at least)
!            longer than limit_in.
! itv      : itv*limin_in is used to cutout the data in process of data-reading.
!            It is needed when the data is bad (i.e., One of the intervals is 
!            longer than itv*limin_in) so that the interpolation may produce 
!            large error).
!            Too small itv may cause short-range profile.
!            Too large itv may cause large error in the interpolation process.
! dz       : interval for interpolating
! zbst     : bottom height for analysis
! ztst     : top height for analysis
! zbuf     : buffer height for polynomial fitting to obtain basic states
! order_p  : order of polynomials used to obtain basic states
!-------------------------------------------------------------------------------

  integer                         ::  nlev, nzt, nzw, nzp, nn, nst

  real, allocatable, dimension(:) ::  hmt, hmw, hmp, tm, um, vm, pm

  real, allocatable, dimension(:) ::  ft, fu, fv, fz, fp

  real, allocatable, dimension(:) ::  zst, tst, ust, vst
  real, allocatable, dimension(:) ::  tstbas, ustbas, vstbas, rstbas, bvst
  real, allocatable, dimension(:) ::  tstprt, ustprt, vstprt

  character (len=256)             ::  ifname, tmpch
  character (len=128)             ::  tit1, tit2, tit3, tit4
  character (len=128)             ::  wfn1, wfn2, wfn3, wfn4
  character (len=4)               ::  arg1, arg2, arg3, arg4
  character (len=5)               ::  stid

! local var.s
  integer                         ::  i,j,k
  integer                         ::  istat, tmpi
  integer                         ::  kzb, kzt
  real, allocatable, dimension(:) ::  h, t, w, deg, p
  real, allocatable, dimension(:) ::  pst, pbas
  real                            ::  p0s, pbasavg, pavg

! buffer-related var.s
  integer                         ::  nst2, kzb2, kzt2, k1
  real                            ::  zb2, zt2
  real, allocatable, dimension(:) ::  z2, t2, u2, v2
  real, allocatable, dimension(:) ::  tbas2, ubas2, vbas2, bv2
  real, allocatable, dimension(:) ::  tprt2, uprt2, vprt2

! variables for generalized linear regression
  integer, parameter              ::  ma = order_p+1
  integer                         ::  ia(ma)
  real                            ::  chisq, a(ma),covar(ma,ma),afunc(50)
  real, allocatable, dimension(:) ::  sig2

! command-line arguments
  integer  ::  iargc
  external     getarg, iargc

! constants
  real, parameter ::  rd = 287.0, cp = 1004.0, g = 9.8
  real, parameter ::  MIS = 1.e32
  real            ::  pi

!-------------------------------------------------------------------------------
!  GET THE ARGUMENTS     
!-------------------------------------------------------------------------------

  call getarg(1,stid)
  call getarg(2,arg1)
  call getarg(3,arg2)
  call getarg(4,arg3)

!-------------------------------------------------------------------------------
!  INPUT DATA NAME       
!-------------------------------------------------------------------------------

  write(ifname,'(a,a2,a2,a2,a)')                              &
       '/export30/atmos_data/SOD/'//trim(stid)//              &
       '/RS_AQC_'//trim(stid)//'_2007',arg1,arg2,arg3,'00.txt'
  print*
  print*, trim(ifname)

!-------------------------------------------------------------------------------
!  ONPUT DATA NAME       
!-------------------------------------------------------------------------------

  write(tit1,'(a)') 'Original Temperature '

  write(wfn1,'(a,a5,a,a2,a2,a2,a)') &
       '/usr/users/miok/KEOP/',stid,'/ort',arg1,arg2,arg3,'.nc'

  write(tit2,'(a)') 'Original WIND'

  write(wfn2,'(a,a5,a,a2,a2,a2,a)') &
       '/usr/users/miok/KEOP/',stid,'/orw',arg1,arg2,arg3,'.nc'

  write(tit3,'(a)') 'Interpolated Var.s'

  write(wfn3,'(a,a5,a,a2,a2,a2,a)') &
       '/usr/users/miok/KEOP/',stid,'/int',arg1,arg2,arg3,'.nc'

  write(tit4,'(a)') 'Basic-state and perturbation Var.s'

  write(wfn4,'(a,a5,a,a2,a2,a2,a)') &
       '/usr/users/miok/KEOP/',stid,'/st',arg1,arg2,arg3,'.nc'

!-------------------------------------------------------------------------------

  pi = 4.*atan(1.0)

  if (iargc() /= 4) then
    print*, 'PROGRAM  STOP'
    print*, 'analyod <stid> <month> <day> <hour>'
    STOP
  end if                                        
                        
!-------------------------------------------------------------------------------
!  READ INPUT DATA
!-------------------------------------------------------------------------------

  open(11,file=ifname,action='read',form='formatted')
  read(11,*) tmpch                                             ! table-head line
  nlev = 0
  do i=1, maxl
    read(11,*,iostat=istat) tmpch
    if (istat == 0) then
      nlev = nlev + 1
    else
      EXIT
    end if
  enddo
  close(11)

  allocate( h(nlev),p(nlev),t(nlev),w(nlev),deg(nlev) )
  h(:) = MIS_in ; p(:)=MIS_in ; t(:)=MIS_in ; w(:)=MIS_in ; deg(:)=MIS_in

  open(11,file=ifname,action='read',form='formatted')
  read(11,*) tmpch                                             ! table-head line
  do i=1, nlev
    read(11, '(13x,f8.2,2x,f8.1,4x,f6.1,15x,f5.1,4x,f6.1)', iostat=istat) &
        p(i), h(i), t(i), w(i), deg(i)
    if ( p(i)==h(i) .and. p(i)==t(i) ) then  ! i.e., absence of data in the line
      h  (i) = MIS_in
      p  (i) = MIS_in
      t  (i) = MIS_in
      w  (i) = MIS_in
      deg(i) = MIS_in
    end if
  enddo
  close(11)

  ! confirm the height culumn not to have missing values
  tmpi = 0
  do i=1, nlev
    if (h(i) /= MIS_in) then
      tmpi = tmpi + 1
      h  (tmpi) = h  (i)
      p  (tmpi) = p  (i)
      t  (tmpi) = t  (i)
      w  (tmpi) = w  (i)
      deg(tmpi) = deg(i)
    end if
  enddo
  nlev = tmpi


  tmpi = 2
  do i=2, nlev

    ! use the data sparsely, if it's interval is too short.
    ! Reversed data (i.e., h(i) > h(i+1)) are also filtered.
    if ( (h(i)-h(tmpi-1)) >= limit_in ) then
      tmpi = tmpi + 1
      if (tmpi .gt. i)  tmpi = i
    end if
    h  (tmpi) = h  (i)
    p  (tmpi) = p  (i)
    t  (tmpi) = t  (i)
    w  (tmpi) = w  (i)
    deg(tmpi) = deg(i)

    ! cutout the data, if it's interval is too long.
    if ( (h(tmpi)-h(tmpi-1)) > (itv*limit_in) ) then
      tmpi = tmpi - 1
      print*, 'DATA INTERVAL IS TOO LONG near Z =', h(tmpi)
      print*, 'INTERVAL :', h(tmpi+1)-h(tmpi)
      EXIT
    end if

  enddo
  if ( (h(tmpi)-h(tmpi-1)) < limit_in ) tmpi = tmpi - 1

  nlev = tmpi

  if (nlev < 20)  then
    print*, 'The number of valid data is smaller then 20'
    STOP
  end if


  ! missing values
  h  (nlev+1:) = MIS
  p  (nlev+1:) = MIS
  t  (nlev+1:) = MIS
  w  (nlev+1:) = MIS
  deg(nlev+1:) = MIS

  do i=1, nlev
    if (p  (i) == MIS_in)  p  (i) = MIS
    if (t  (i) == MIS_in)  t  (i) = MIS
    if (w  (i) == MIS_in)  w  (i) = MIS
    if (deg(i) == MIS_in)  deg(i) = MIS
  end do
 
!-------------------------------------------------------------------------------
!  REMOVE MISSING VALUES
!-------------------------------------------------------------------------------

  ! w.r.t. T

  tmpi = 0
  do i=1, nlev
    if (t(i) /= MIS) then
      tmpi = tmpi + 1
    end if
  end do

  nzt = tmpi

  if (nzt == 0) then
    print*, 'TEMPERATURE DATA DOES NOT EXIST' 
    STOP       
  end if

  allocate(hmt(nzt))         ; hmt = 0.    
  allocate(tm (nzt))         ; tm  = 0.    

  tmpi = 0
  do i=1, nlev
    if (t(i) /= MIS) then
      tmpi = tmpi + 1
      tm (tmpi) = t(i) + 273.
      hmt(tmpi) = h(i)
    end if
  end do

  deallocate(t)

!-------------------------------------------------------------------------------

  ! w.r.t. wind

  tmpi = 0                          
  do i=1, nlev
    if (w(i) /= MIS) then
      tmpi = tmpi + 1
    end if
  end do

  nzw = tmpi

  if (nzw == 0) then
    print*, 'WIND DATA DOES NOT EXIST' 
    STOP       
  end if

  allocate(hmw(nzw))         ; hmw = 0.
  allocate(um (nzw))         ; um  = 0.
  allocate(vm (nzw))         ; vm  = 0.

  tmpi = 0
  do i=1, nlev
    if (w(i) /= MIS) then
      tmpi = tmpi + 1
      um (tmpi) = -w(i)*sin(deg(i)*pi/180.)
      vm (tmpi) = -w(i)*cos(deg(i)*pi/180.)
      hmw(tmpi) = h(i)
    end if
  end do

  deallocate(w,deg)

!-------------------------------------------------------------------------------

  ! w.r.t. p

  tmpi = 0
  do i=1, nlev
    if (p(i) /= MIS) then
      tmpi = tmpi + 1
    end if
  end do

  nzp = tmpi

  if (nzp == 0) then
    print*, 'PRESSURE DATA DOES NOT EXIST'
    STOP
  end if

  allocate(hmp(nzp))         ; hmp = 0.
  allocate(pm (nzp))         ; pm  = 0.

  tmpi = 0
  do i=1, nlev
    if (p(i) /= MIS) then
      tmpi = tmpi + 1
      pm (tmpi) = p(i)
      hmp(tmpi) = h(i)
    end if
  end do

  deallocate(p)

!-------------------------------------------------------------------------------

  deallocate(h)

!-------------------------------------------------------------------------------
!  OUTPUT THE FILTERED RAWDATA
!-------------------------------------------------------------------------------

  call out1d(trim(wfn1),1,(/'T'/),tm,'Z',nzt,hmt,trim(tit1))
  
  call out1d(trim(wfn2),2,(/'U','V'/),(/um,vm/),'Z',nzw,hmw,trim(tit2))

!-------------------------------------------------------------------------------
!  INTERPORATION 
!-------------------------------------------------------------------------------

  ! The extrapolated data below the bottom of raw data is not valid, 
  ! but contained for convenience (i.e., All profiles start from z = 0).

  if (hmt(nzt) < hmw(nzw)) then
    nn = int(hmt(nzt)/dz)+1
  else
    nn = int(hmw(nzw)/dz)+1
  end if

  allocate(fz(nn)) 
  do k=1, nn 
    fz(k) = dz * (k-1)
  end do

  ! interporation of temperature data
  allocate(ft(nn))
  call cs_int(nzt,hmt,tm,nn,fz,ft)
  deallocate(hmt,tm)

  ! interporation of wind data
  allocate(fu(nn))
  call cs_int(nzw,hmw,um,nn,fz,fu)
  deallocate(um)

  allocate(fv(nn))
  call cs_int(nzw,hmw,vm,nn,fz,fv)
  deallocate(vm)

  deallocate(hmw)

  ! interporation of pressure data
  allocate(fp(nn))
  call cs_int(nzp,hmp,pm,nn,fz,fp)
  deallocate(hmp,pm)

!-------------------------------------------------------------------------------
!  OUTPUT INTERPOLATED DATA
!-------------------------------------------------------------------------------
  
  call out1d(trim(wfn3),3,(/'T','U','V'/),(/ft,fu,fv/),'Z',nn,fz,trim(tit3))


! ===== YOU NEED TO CHECK THE INTERPORATED DATA ================================


! END OF THE DATA TREATMENT

!===============================================================================


  if (fz(nn) < ztst) then
    print*, 'Top of the data is lower than critical value.'
    print*, 'Analysis does not performed.'
    STOP
  end if


!-------------------------------------------------------------------------------
!  SELECT DATA RANGE FOR ANALYSIS
!-------------------------------------------------------------------------------

  zb2 = zbst - zbuf
  zt2 = min(fz(nn),ztst+zbuf)
         
  kzb = int(zbst/dz + 0.001) + 1
  kzt = int(ztst/dz + 0.001) + 1

  kzb2 = int(zb2/dz + 0.001) + 1
  kzt2 = int(zt2/dz + 0.001) + 1

  nst  = kzt  - kzb  + 1
  nst2 = kzt2 - kzb2 + 1

  allocate(zst(nst))          ; zst = 0.
  allocate(tst(nst))          ; tst = 0.
  allocate(ust(nst))          ; ust = 0.
  allocate(vst(nst))          ; vst = 0.

  allocate(tstbas(nst))       ; tstbas = 0.
  allocate(ustbas(nst))       ; ustbas = 0.
  allocate(vstbas(nst))       ; vstbas = 0.
  allocate(rstbas(nst))       ; rstbas = 0.
  allocate(bvst  (nst))       ; bvst   = 0.

  allocate(tstprt(nst))       ; tstprt = 0.
  allocate(ustprt(nst))       ; ustprt = 0.
  allocate(vstprt(nst))       ; vstprt = 0.

  allocate(z2(nst2))          ; z2 = 0.
  allocate(t2(nst2))          ; t2 = 0.
  allocate(u2(nst2))          ; u2 = 0.
  allocate(v2(nst2))          ; v2 = 0.

  allocate(tbas2(nst2))       ; tbas2 = 0.
  allocate(ubas2(nst2))       ; ubas2 = 0.
  allocate(vbas2(nst2))       ; vbas2 = 0.
  allocate(bv2  (nst2))       ; bv2   = 0.

  allocate(tprt2(nst2))       ; tprt2 = 0.
  allocate(uprt2(nst2))       ; uprt2 = 0.
  allocate(vprt2(nst2))       ; vprt2 = 0.

  allocate(sig2(nst2))        ; sig2 = 0.

  allocate(pst (nst))         ; pst  = 0.
  allocate(pbas(nst))         ; pbas = 0.

  do k=1, nst
    zst(k) = (k-1)*dz + zbst
    pst(k) = fp(k-1+kzb)
  end do

  ! include buffer 
  do k=1, nst2
    z2(k) = (k-1)*dz + zb2
    t2(k) = ft(k-1+kzb2)
    u2(k) = fu(k-1+kzb2)
    v2(k) = fv(k-1+kzb2)
  end do

!-------------------------------------------------------------------------------
!  OBTAIN THE BASIC STATE AND PERTURBATION USING POLYNOMIAL FITTING
!-------------------------------------------------------------------------------

  ! fit T
  ia(1:ma) = 1
  a(1:ma) = 1.
  sig2(1:nst2) = 1.
  call lfit(z2,t2,sig2,nst2,a,ia,ma,covar,ma,chisq)

  do k=1, nst2 
    tbas2(k) = 0.
    call funcs(z2(k),afunc,ma)
    do i=1, ma
      tbas2(k) = tbas2(k) + a(i)*afunc(i)
    end do
    tprt2(k) = t2(k) - tbas2(k)
  end do 

  ! fit u
  ia(1:ma) = 1
  a(1:ma) = 1.
  sig2(1:nst2) = 1.
  call lfit(z2,u2,sig2,nst2,a,ia,ma,covar,ma,chisq)

  do k=1, nst2
    ubas2(k) = 0.
    call funcs(z2(k),afunc,ma)
    do i=1, ma
      ubas2(k) = ubas2(k) + a(i)*afunc(i)
    end do
    uprt2(k) = u2(k) - ubas2(k)
  end do

  ! fit v
  ia(1:ma) = 1
  a(1:ma) = 1.
  sig2(1:nst2) = 1.
  call lfit(z2,v2,sig2,nst2,a,ia,ma,covar,ma,chisq)
    
  do k=1, nst2
    vbas2(k) = 0.
    call funcs(z2(k),afunc,ma)
    do i=1, ma
      vbas2(k) = vbas2(k) + a(i)*afunc(i)
    end do
    vprt2(k) = v2(k) - vbas2(k)
  end do

  ! Brunt-Vaisala freq.
  do k=2, nst2-1
    bv2(k) = ((tbas2(k+1)-tbas2(k-1))/(2.0*dz)+g/cp)*g/tbas2(k)
    if(bv2(k) .lt. 0.) then
      bv2(k) = 0.
    else
      bv2(k) = sqrt(bv2(k))
    end if
  end do
  bv2(1)    = bv2(2)
  bv2(nst2) = bv2(nst2-1)

  ! exclude buffer
  do k=1, nst2
    if (z2(k) == zst(1))  k1 = k
  end do

  do k=1, nst
    tst   (k) = t2   (k+k1-1)
    ust   (k) = u2   (k+k1-1)
    vst   (k) = v2   (k+k1-1)
    tstbas(k) = tbas2(k+k1-1)
    ustbas(k) = ubas2(k+k1-1)
    vstbas(k) = vbas2(k+k1-1)
    bvst  (k) = bv2  (k+k1-1)
    tstprt(k) = tprt2(k+k1-1)
    ustprt(k) = uprt2(k+k1-1)
    vstprt(k) = vprt2(k+k1-1)
  end do

  ! calculate pbar using hydrostatic relation
  pbas(1) = 1.
  pbasavg = 0.
  do k=2, nst
    pbas(k) = pbas(k-1) * exp(-g/rd*2./(tstbas(k)+tstbas(k-1))*dz)
    pbasavg = pbasavg + log(pbas(k))
  enddo
  pbasavg = pbasavg / nst

  ! determine p0
  pavg = 0.
  do k=1, nst
    pavg = pavg + log(pst(k))
  enddo
  pavg = pavg / nst

  p0s = exp(pavg-pbasavg)

  pbas(:) = pbas(:) * p0s

  if (abs(pst(1)-pbas(1))/pbas(1).gt.0.2 .or. &
      abs(pst(nst)-pbas(nst))/pbas(nst).gt.0.2) then
    print*, 'PRESSURE DATA and HYDROSTATIC PRESSURE are largely diff.'
    print*, 'sfc:', pst(1), pbas(1)
    print*, 'top:', pst(nst), pbas(nst)
  end if

  ! calculate rhobar
  rstbas(:) = 100.*pbas(:)/rd/tstbas(:)


  deallocate(z2,t2,u2,v2,pst)
  deallocate(tbas2,ubas2,vbas2,pbas,bv2)
  deallocate(tprt2,uprt2,vprt2)
  deallocate(sig2)

  print *, 'nz_anal : ', nst


!-------------------------------------------------------------------------------
!  OUPUT BASIC-STATE AND PERTURBATION VAR.S
!-------------------------------------------------------------------------------

  call out1d(trim(wfn4),11,(/'T   ','U','V','Tbas','Ubas','Vbas','RHO0','N',   &
                             'Tper','Uper','Vper'/),                           &
             (/tst,ust,vst,tstbas,ustbas,vstbas,rstbas,bvst,                   &
               tstprt,ustprt,vstprt/),                                         &
             'Z',nst,zst,trim(tit4))


  STOP

END program


SUBROUTINE cs_int(nxi,xi,yi,nxo,xo,yo)

  use cubicspline

  implicit none

  integer, intent(in)  :: nxi, nxo
  real,    intent(in)  :: xi(nxi), yi(nxi), xo(nxo)
  real,    intent(out) :: yo(nxo)

  integer :: i,j,k
  real    :: tmp(nxi)


  call spline(xi,yi,nxi,1.e30,1.e30,tmp)
  do k=1, nxo
    call splint(xi,yi,tmp,nxi,xo(k),yo(k))
  enddo


  RETURN

END subroutine

