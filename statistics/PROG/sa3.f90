program spectral_analysis_3

  use netcdfio
  use pwrspd
  use cubicspline

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxa = 320, nya = 160, nxo = 320, nyo = 320
  integer, parameter :: nt1all = 1200, nt2all = 960
  integer, parameter :: ts1 = 1, ts2 = 1, nt1 = nt1all, nt2 = nt2all
  integer, parameter :: nfu1 = 5        ! CTL (pre-ind.)
  integer, parameter :: nfu2 = 4        ! CO2
  integer, parameter :: nxp = 134, nypa0 = 144, nypa = 80
!--------------------------------------------------------------
  real,    parameter :: dx = 360./nxo, dt = 1., dy = 1.
!--------------------------------------------------------------
  integer, parameter :: nkkp = nxp/2+1, nkka = nxa/2+1, nll = nypa/2+1
  integer, parameter :: noo1 = nt1/2*2+1, noo2 = nt2/2*2+1

  integer :: i,j,n, nf, ncid, ii
  real    :: pi

  real, dimension(nxo,nyo,nt1) :: to1
  real, dimension(nxo,nyo,nt2) :: to2
  real, dimension(nxa,nt1)     :: ua1
  real, dimension(nxa,nt2)     :: ua2

  real, dimension(nxp,nt1)     :: topac1
  real, dimension(nxp,nt2)     :: topac2
  real, dimension(nypa0,nt1)   :: sst0p1, sst0a1
  real, dimension(nypa0,nt2)   :: sst0p2, sst0a2
  real, dimension(nypa,nt1)    :: sstp1, ssta1
  real, dimension(nypa,nt2)    :: sstp2, ssta2

  real, dimension(nkkp,noo1)   :: psdp1, psdp1r
  real, dimension(nkkp,noo2)   :: psdp2, psdp2r
  real, dimension(nkka,noo1)   :: psda1, psda1r
  real, dimension(nkka,noo2)   :: psda2, psda2r
  real, dimension(nll,noo1)    :: psdpac1, psdatl1, psdocn1r
  real, dimension(nll,noo2)    :: psdpac2, psdatl2, psdocn2r

  real    :: lona(nxa), lata(nya), lono(nxo), lato(nyo), t1(nt1), t2(nt2)
  real    :: lat(nypa)
  real    :: dsinphia(nya), dsinphio(nyo)
  real    :: kkp(nkkp), kka(nkka), oo1(noo1), oo2(noo2), ll(nll)
  integer, dimension(nyo) :: isttp, iendp, istta, ienda

  integer :: tmpi, ntsf
  real*8  :: tmp1(nxa), tmp2(nya), tmp3(nxo), tmp4(nyo)
  real*8  :: tmp5(2,nya), tmp6(2,nyo)
  real    :: tmpr1, tmpr2
  real    :: tmp1d1(nt1), tmp1d2(nt1), tmp1d3(nt2), tmp1d4(nt2)
  real    :: tmp1d5(nxp), tmp1d6(nxp), tmp1d7(nxa), tmp1d8(nxa)
  real    :: tmp1d9(nypa0)
  real, dimension(:,:,:), allocatable :: tmp3d

  character*128 :: fdir(2)
  character*128 :: fn_to1, fn_ua1(nfu1)
  character*128 :: fn_to2, fn_ua2(nfu2)


  pi = acos(-1.)

! read data ------------------------------------------------------------
  fdir(1) = '/export30/kyh/MIROC/pictl/'
  fdir(2) = '/export30/kyh/MIROC/CO2/'
  write(fn_to1,'(a)') trim(fdir(1))//'tos_O1.nc'
  write(fn_ua1(1),'(a)') trim(fdir(1))//'ua_A1_1_20.nc'
  write(fn_ua1(2),'(a)') trim(fdir(1))//'ua_A1_21_40.nc'
  write(fn_ua1(3),'(a)') trim(fdir(1))//'ua_A1_41_60.nc'
  write(fn_ua1(4),'(a)') trim(fdir(1))//'ua_A1_61_80.nc'
  write(fn_ua1(5),'(a)') trim(fdir(1))//'ua_A1_81_100.nc'
  write(fn_to2,'(a)') trim(fdir(2))//'tos_O1.nc'
  write(fn_ua2(1),'(a)') trim(fdir(2))//'ua_A1_1_20.nc'
  write(fn_ua2(2),'(a)') trim(fdir(2))//'ua_A1_21_40.nc'
  write(fn_ua2(3),'(a)') trim(fdir(2))//'ua_A1_41_60.nc' 
  write(fn_ua2(4),'(a)') trim(fdir(2))//'ua_A1_61_80.nc' 

! axis
  call opennc(fn_ua1(1),ncid)
  call dget1d(ncid,'lon',nxa,tmp1)
  call dget1d(ncid,'lat',nya,tmp2)
  call dget2d(ncid,'lat_bnds',2,nya,tmp5)
  call closenc(ncid)
  lona = real(tmp1)
  lata = real(tmp2)
  dsinphia(:) = real(dsin(pi/180.*tmp5(2,:))-dsin(pi/180.*tmp5(1,:)))

  call opennc(fn_to1,ncid)
  call dget1d(ncid,'lon',nxo,tmp3)
  call dget1d(ncid,'lat',nyo,tmp4)
  call dget2d(ncid,'lat_bnds',2,nyo,tmp6)
  call closenc(ncid)
  lono = real(tmp3)
  lato = real(tmp4)
  dsinphio(:) = real(dsin(pi/180.*tmp6(2,:))-dsin(pi/180.*tmp6(1,:)))

  do n=1, nt1
    t1(n) = real(n-1)/12.
  enddo
  do n=1, nt2
    t2(n) = real(n-1)/12.
  enddo

! variables
  call opennc(fn_to1,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts1,nt1,to1)
  call closenc(ncid)
  call opennc(fn_to2,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts2,nt2,to2)
  call closenc(ncid)

  tmpi = 0
  do nf=1, nfu1
    call opennc(fn_ua1(nf),ncid)
    call dilen(ncid,'time',ntsf)
    allocate(tmp3d(nxa,2,ntsf))             ;  tmp3d = 0.
    call geta4d(ncid,'ua',1,nxa,80,2,10,1,1,ntsf,tmp3d(:,1:2,:))
    call closenc(ncid)
    do n=1, ntsf
      tmpi = tmpi + 1
      ua1(:,tmpi) = (tmp3d(:,1,n)+tmp3d(:,2,n)) / 2.
    enddo
    deallocate(tmp3d)
  enddo

  tmpi = 0
  do nf=1, nfu2
    call opennc(fn_ua2(nf),ncid)
    call dilen(ncid,'time',ntsf)
    allocate(tmp3d(nxa,2,ntsf))             ;  tmp3d = 0.
    call geta4d(ncid,'ua',1,nxa,80,2,10,1,1,ntsf,tmp3d(:,1:2,:))
    call closenc(ncid)
    do n=1, ntsf
      tmpi = tmpi + 1
      ua2(:,tmpi) = (tmp3d(:,1,n)+tmp3d(:,2,n)) / 2.
    enddo
    deallocate(tmp3d)
  enddo


! (1) SST over the Pacific ocean ---------------------------------------
  topac1 = 0.
  do n=1, nt1
  do i=1, nxp
    tmpr1 = 0.
    tmpr2 = 0.
    do j=152, 169
      if (to1(i+107,j,n) .ne. 1.e20) then
        tmpr1 = tmpr1 + to1(i+107,j,n) * dsinphio(j)
        tmpr2 = tmpr2 + dsinphio(j)
      end if
    enddo
    topac1(i,n) = tmpr1 / tmpr2
  enddo
  enddo

  topac2 = 0.
  do n=1, nt2
  do i=1, nxp
    tmpr1 = 0.
    tmpr2 = 0.
    do j=152, 169
      if (to2(i+107,j,n) .ne. 1.e20) then
        tmpr1 = tmpr1 + to2(i+107,j,n) * dsinphio(j)
        tmpr2 = tmpr2 + dsinphio(j)
      end if
    enddo
    topac2(i,n) = tmpr1 / tmpr2
  enddo
  enddo

  call out2d2('../sa3pac.nc',1,(/'pic'/),topac1,'lon1',nxp,lono(108:107+nxp), &
              't1',nt1,t1, 1,(/'CO2'/),topac2,'lon2',nxp,lono(108:107+nxp), &
              't2',nt2,t2, 'Pacific SST')

  ! 2-D linear trend
  tmp1d1 = 0.
  do i=1, nxp
    tmp1d1(:) = tmp1d1(:) + topac1(i,:)/nxp
  enddo
  call detrend1(nt1, t1, tmp1d1, tmp1d2, 0)
  do n=1, nt1
    topac1(:,n) = topac1(:,n) - tmp1d2(n)
  enddo

  tmp1d5 = 0.
  do n=1, nt1
    tmp1d5(:) = tmp1d5(:) + topac1(:,n)/nt1
  enddo
  call detrend1(nxp, lono(108:107+nxp), tmp1d5, tmp1d6, 0)
  do i=1, nxp
    topac1(i,:) = topac1(i,:) - tmp1d6(i)
  enddo

  tmp1d3 = 0.
  do i=1, nxp
    tmp1d3(:) = tmp1d3(:) + topac2(i,:)/nxp
  enddo
  call detrend1(nt2, t2, tmp1d3, tmp1d4, 0)
  do n=1, nt2
    topac2(:,n) = topac2(:,n) - tmp1d4(n)
  enddo

  tmp1d5 = 0.
  do n=1, nt2
    tmp1d5(:) = tmp1d5(:) + topac2(:,n)/nt2
  enddo
  call detrend1(nxp, lono(108:107+nxp), tmp1d5, tmp1d6, 0)
  do i=1, nxp
    topac2(i,:) = topac2(i,:) - tmp1d6(i)
  enddo


  call psd2d(nxp,nt1,dx,dt,topac1,kkp,oo1,psdp1r)
  call psd2d(nxp,nt2,dx,dt,topac2,kkp,oo2,psdp2r)
  do n=1, noo1
    psdp1(:,n) = psdp1r(:,noo1+1-n)
  enddo
  do n=1, noo2
    psdp2(:,n) = psdp2r(:,noo2+1-n)
  enddo

  call out2d2('../sa3pacps.nc',1,(/'PS_pic'/),psdp1,'k1',nkkp,kkp, &
              'o1',noo1,oo1, 1,(/'PS_CO2'/),psdp2,'k2',nkkp,kkp, &
              'o2',noo2,oo2, 'Pacific SST')

! (2) 200mb zonal wind -------------------------------------------------
  call out2d2('../sa3ua.nc',1,(/'pic'/),ua1,'lon1',nxa,lona, &
              't1',nt1,t1, 1,(/'CO2'/),ua2,'lon2',nxa,lona, &
              't2',nt2,t2, 'ua')

  ! 1-D linear trend
  tmp1d1 = 0.
  do i=1, nxa
    tmp1d1(:) = tmp1d1(:) + ua1(i,:)/nxa
  enddo
  call detrend1(nt1, t1, tmp1d1, tmp1d2, 0)
  do n=1, nt1
    ua1(:,n) = ua1(:,n) - tmp1d2(n)
  enddo

  tmp1d3 = 0.
  do i=1, nxa
    tmp1d3(:) = tmp1d3(:) + ua2(i,:)/nxa
  enddo
  call detrend1(nt2, t2, tmp1d3, tmp1d4, 0)
  do n=1, nt2
    ua2(:,n) = ua2(:,n) - tmp1d4(n)
  enddo


  call psd2d(nxa,nt1,dx,dt,ua1,kka,oo1,psda1r)
  call psd2d(nxa,nt2,dx,dt,ua2,kka,oo2,psda2r)
  do n=1, noo1
    psda1(:,n) = psda1r(:,noo1+1-n)
  enddo
  do n=1, noo2 
    psda2(:,n) = psda2r(:,noo2+1-n)
  enddo

  call out2d2('../sa3uaps.nc',1,(/'PS_pic'/),psda1,'k1',nkka,kka, &
              'o1',noo1,oo1, 1,(/'PS_CO2'/),psda2,'k2',nkka,kka, &
              'o2',noo2,oo2, 'Zonal wind')

! (3) N-S spectral analysis --------------------------------------------
  do j=1, nypa
    lat(j) = (j-40.5) * dy
  enddo

  isttp(89:113)  = 132 
  isttp(114:143) = 120
  isttp(144:152) = 114
  isttp(153:161) = 108
  isttp(162:170) = 102
  isttp(171:180) = 96
  isttp(181:217) = 94
  isttp(218:232) = 126
  iendp(89:176)  = 264
  iendp(177:187) = 246
  iendp(188:232) = 233

  istta(89:96)   = 259
  istta(97:178)  = 265
  istta(179:187) = 246
  istta(188:192) = 239
  istta(193:232) = 233
  ienda(89:200)  = 20
  ienda(201:232) = 314
  do j=89, 232
    if ( ienda(j) .lt. istta(j) )  ienda(j) = ienda(j) + nxo
  enddo


  do n=1, nt1
    do j=1, nypa0

      tmpi  = 0
      tmpr1 = 0.
      do i=isttp(j+88), iendp(j+88)
        if (to1(i,j+88,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to1(i,j+88,n)
        end if
      enddo
      sst0p1(j,n) = tmpr1 / tmpi

      tmpi  = 0
      tmpr1 = 0.
      do i=istta(j+88), ienda(j+88)
        ii = i
        if (ii .gt. nxo)  ii = ii - nxo
        if (to1(ii,j+88,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to1(ii,j+88,n) 
        end if
      enddo
      sst0a1(j,n) = tmpr1 / tmpi

    enddo

    call spline(lato(89:88+nypa0),sst0p1(:,n),nypa0,1.e+30,1.e+30,tmp1d9)
    do j=1, nypa
      call splint(lato(89:88+nypa0),sst0p1(:,n),tmp1d9,nypa0,lat(j),sstp1(j,n))
    enddo
    call spline(lato(89:88+nypa0),sst0a1(:,n),nypa0,1.e+30,1.e+30,tmp1d9)
    do j=1, nypa
      call splint(lato(89:88+nypa0),sst0a1(:,n),tmp1d9,nypa0,lat(j),ssta1(j,n))
    enddo

  enddo


  do n=1, nt2
    do j=1, nypa0

      tmpi  = 0
      tmpr1 = 0.
      do i=isttp(j+88), iendp(j+88)
        if (to2(i,j+88,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to2(i,j+88,n) 
        end if 
      enddo 
      sst0p2(j,n) = tmpr1 / tmpi 
 
      tmpi  = 0 
      tmpr1 = 0. 
      do i=istta(j+88), ienda(j+88)
        ii = i 
        if (ii .gt. nxo)  ii = ii - nxo
        if (to2(ii,j+88,n) .ne. 1.e20) then 
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to2(ii,j+88,n)
        end if
      enddo
      sst0a2(j,n) = tmpr1 / tmpi

    enddo

    call spline(lato(89:88+nypa0),sst0p2(:,n),nypa0,1.e+30,1.e+30,tmp1d9)
    do j=1, nypa
      call splint(lato(89:88+nypa0),sst0p2(:,n),tmp1d9,nypa0,lat(j),sstp2(j,n))
    enddo
    call spline(lato(89:88+nypa0),sst0a2(:,n),nypa0,1.e+30,1.e+30,tmp1d9)
    do j=1, nypa
      call splint(lato(89:88+nypa0),sst0a2(:,n),tmp1d9,nypa0,lat(j),ssta2(j,n))
    enddo

  enddo

  call out2d2('../sa3ocn2.nc',2,(/'Pac_pic','Atl_pic'/),(/sstp1,ssta1/), &
              'lat1',nypa,lat,'t1',nt1,t1, &
              2,(/'Pac_CO2','Atl_CO2'/),(/sstp2,ssta2/),'lat2',nypa,lat, &
              't2',nt2,t2, 'Pacific and Atlantic Ocean')

  ! 1-D linear trend
  tmp1d1 = 0.
  do j=1, nypa
    tmp1d1(:) = tmp1d1(:) + sstp1(j,:)/nypa
  enddo
  call detrend1(nt1, t1, tmp1d1, tmp1d2, 0)
  do n=1, nt1
    sstp1(:,n) = sstp1(:,n) - tmp1d2(n)
  enddo

  tmp1d1 = 0.
  do j=1, nypa
    tmp1d1(:) = tmp1d1(:) + ssta1(j,:)/nypa
  enddo
  call detrend1(nt1, t1, tmp1d1, tmp1d2, 0)
  do n=1, nt1
    ssta1(:,n) = ssta1(:,n) - tmp1d2(n)
  enddo

  tmp1d3 = 0.
  do j=1, nypa
    tmp1d3(:) = tmp1d3(:) + sstp2(j,:)/nypa
  enddo
  call detrend1(nt2, t2, tmp1d3, tmp1d4, 0)
  do n=1, nt2
    sstp2(:,n) = sstp2(:,n) - tmp1d4(n)
  enddo

  tmp1d3 = 0.
  do j=1, nypa
    tmp1d3(:) = tmp1d3(:) + ssta2(j,:)/nypa
  enddo
  call detrend1(nt2, t2, tmp1d3, tmp1d4, 0)
  do n=1, nt2
    ssta2(:,n) = ssta2(:,n) - tmp1d4(n)
  enddo


  call psd2d(nypa,nt1,dy,dt,sstp1,ll,oo1,psdocn1r)
  call psd2d(nypa,nt2,dy,dt,sstp2,ll,oo2,psdocn2r) 
  do n=1, noo1
    psdpac1(:,n) = psdocn1r(:,noo1+1-n)
  enddo 
  do n=1, noo2  
    psdpac2(:,n) = psdocn2r(:,noo2+1-n)
  enddo

  call psd2d(nypa,nt1,dy,dt,ssta1,ll,oo1,psdocn1r)
  call psd2d(nypa,nt2,dy,dt,ssta2,ll,oo2,psdocn2r)
  do n=1, noo1
    psdatl1(:,n) = psdocn1r(:,noo1+1-n)
  enddo
  do n=1, noo2   
    psdatl2(:,n) = psdocn2r(:,noo2+1-n)
  enddo

  call out2d2('../sa3ocn2ps.nc',2,(/'Pac_pic','Atl_pic'/), &
              (/psdpac1,psdatl1/),'l1',nll,ll,'o1',noo1,oo1, &
              2,(/'Pac_CO2','Atl_CO2'/),(/psdpac2,psdatl2/),'l2',nll,ll, &
              'o2',noo2,oo2, 'PSD for two oceans')


end program


subroutine detrend1(nt, t, x, trend, opt)

  use regress

  integer,             intent(in)    :: nt, opt
  real, dimension(nt), intent(in)    :: t
  real, dimension(nt), intent(inout) :: x
  real, dimension(nt), intent(out)   :: trend

  integer :: i,j
  real    :: a(2), afunc(50), sig(nt), covar(12,12), chisq

  a(:)   = 0.0         ! coeff. for ia=0
  sig(:) = 1.0 
  call lfit(t,x,sig,nt,a,(/1,1/),2,covar,12,chisq)
  trend = 0.
  do i=1, nt 
    call funcs(t(i),afunc,2)
    do j=1, 2 
      trend(i) = trend(i) + a(j)*afunc(j)
    end do
  end do
  if (opt .eq. 1) then 
    x(:) = x(:) - trend(:)
  end if


end subroutine detrend1 


subroutine detrend_season(nt,x)

  integer, intent(in)    :: nt
  real,    intent(inout) :: x(nt)

  integer :: i,j,n, nn(12)
  real    :: ave(12)

  n   = 0
  nn  = 0
  ave = 0.
  do i=1, nt
    n = n + 1
    ave(n) = ave(n) + x(i)
    nn(n) = nn(n) + 1
    if (n .eq. 12)  n = 0
  enddo
  ave(:) = ave(:) / nn(:)

  n = 0
  do i=1, nt
    n = n + 1
    x(i) = x(i) - ave(n)
    if (n .eq. 12)  n = 0
  enddo

end subroutine detrend_season

