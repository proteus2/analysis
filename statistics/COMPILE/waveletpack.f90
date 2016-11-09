MODULE waveletpack

!--------------------------------------------------------
!
! pad0  :  padding with zeros. (for non-periodic data)
! ma    :  ( order of polynomial + 1 ) for de-trending
! param :  >= 5 for Morlet wavelet (default : 6)
!
!--------------------------------------------------------

  use netcdfio
  use regress    , only: lfit  ,funcs

  public :: wlncpack, wlpack

  contains


subroutine wlncpack(ofname,ndata,dt,x,ydata,ma,nj,s0,subscale,w_ftn,param, &
               pad0,savgopt,avgscale,sigopt,redopt,siglvl,reconopt)

  implicit none

  integer, intent(in) :: ndata, ma, nj, subscale, w_ftn, pad0
  integer, intent(in) :: savgopt, sigopt, redopt, reconopt
  real,    intent(in) :: dt, s0, param, avgscale(2), siglvl
  real, dimension(ndata), intent(in) :: x, ydata
  character(len=*),       intent(in) :: ofname

  real, dimension(nj,ndata+1) :: psd
  real, dimension(ndata)      :: pwr_savg
  real, dimension(nj)         :: fscale, ps_theor, signif, sig_glob
  real                        :: mean_savg, sig_savg

  integer :: i,j, ndpad, javg1, javg2, ncid, ia(ma)
  double complex, dimension(ndata,nj) :: wcoef
  double precision, dimension(ndata) :: ydata8, yrecon, coi
  double precision, dimension(nj)    :: wscale, fscale8
  double precision, dimension(nj)    :: dof, ps_theor8, signif8, mean_savg8
  double precision, dimension(nj)    :: sig_glob8, dof_savg, sig_savg8
  double precision :: mean_recon, var_recon, rms_recon
  double precision :: pi, lag1, ymean, variance, cdelta, psi0
  double precision :: dt8, dj, param8, siglvl8
  real, dimension(ndata) :: sig_det, func, trend
  real :: a(ma), afunc(50), chisq_det, cov(ma+10,ma+10)
  real :: lag1_4, x_out(ndata+1)


  pi = 4.d0 * atan(1.d0)
  dj = 1.d0/subscale
  ydata8 = dble(ydata)
  dt8 = dble(dt)
  param8 = dble(param)
  
  ndpad = int(2.d0**(nint(log(dble(ndata))/log(2.d0))+1))
  if (pad0 .eq. 0) then
    ndpad = ndata
  endif
  
  
! Remove trend in data -------------------------------------------------
  trend = 0.
  if (ma .gt. 1) then
    ia(:) = 1           ! calculate the coeff.(=a) if ia.ne.0
    a(:)  = 0.0         ! coeff. for ia=0
    sig_det(:) = 1.0
    call lfit(x,ydata,sig_det,ndata,a,ia,ma,cov,ma+10,chisq_det)
    do i=1, ndata
      call funcs(x(i),afunc,ma)
      do j=1, ma
        trend(i) = trend(i) + a(j)*afunc(j)
      end do
    end do
  end if

  ydata8 = ydata8 - trend

  print*, ' Edge of the data : ', real(ydata8(1)), real(ydata8(ndata))

! Wavelet transform & calculate the power spectrum ---------------------
  call wavelet(ndata,ydata8,dt8,w_ftn,param8,dble(s0),dj,nj,ndpad, &
               wcoef,wscale,fscale8,coi)
  fscale = real(fscale8)

  psd = 0.
  do i=1, ndata
    psd(:,i) = real(abs(wcoef(i,:))**2)
  enddo

  ! Global averaged power spectrum
  do i=1, ndata
    psd(:,ndata+1) = psd(:,ndata+1) + psd(:,i)
  enddo
  psd(:,ndata+1) = psd(:,ndata+1) / ndata

  ! define constants (only four empirical values)
  cdelta = -1.d0  ;  psi0 = -1.d0
  if ( (w_ftn.eq.0) .and. (param8.eq.6.d0) ) then
    cdelta = 0.776d0  ;  psi0 = pi**(-0.25d0)
  else if ( (w_ftn.eq.1) .and. (int(param8).eq.4) ) then
    cdelta = 1.132d0  ;  psi0 = 1.079d0
  else if ( (w_ftn.eq.2) .and. (int(param8).eq.2) ) then
    cdelta = 3.541d0  ;  psi0 = 0.867d0
  else if ( (w_ftn.eq.2) .and. (int(param8).eq.6) ) then
    cdelta = 1.966d0  ;  psi0 = 0.884d0
  end if

! Scale-averaged power series ------------------------------------------
  pwr_savg = 0.
  if (savgopt .eq. 1) then

    dof_savg(1) = dble(avgscale(1))                ! smallest scale
    dof_savg(2) = dble(avgscale(2))                ! largest scale
    javg1 = 0  ;  javg2 = 0
    do j=1, nj
      if ( (fscale8(j).ge.dof_savg(1)) .and. (javg1.eq.0) )  javg1 = j
      if (fscale8(j) .le. dof_savg(2))  javg2 = j
    enddo
    do j=javg1, javg2
      pwr_savg(:) = pwr_savg(:) + psd(j,:)/real(wscale(j))
    enddo
    pwr_savg(:) = pwr_savg(:)*real(dj*dt8/cdelta)

  end if

! Significant test -----------------------------------------------------
  lag1 = 0.d0
  if (redopt .eq. 1) then
    call autocor1(ndata,real(ydata8),lag1_4)
  end if
  lag1 = dble(lag1_4)

  if (sigopt .eq. 1) then
    siglvl8 = dble(siglvl)
    ! regular chi-square test
    call wave_signif(0,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof,ps_theor8,signif8, &
            ymean,variance,cdelta,psi0)

    ! significance test of the global averaged spectrum
    dof(:) = ndata - wscale(:)
    call wave_signif(1,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof,ps_theor8,sig_glob8, &
            ymean,variance,cdelta,psi0)

    ps_theor = real(ps_theor8)
    signif   = real(signif8)
    sig_glob = real(sig_glob8)

    ! significance test of the scale-averaged time series
    if (savgopt .eq. 1) then
      call wave_signif(2,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof_savg,mean_savg8,sig_savg8, &
            ymean,variance,cdelta,psi0)

      mean_savg = real(mean_savg8(1))
      sig_savg = real(sig_savg8(1))
      print*, ' Scale-avg. Degrees of freedom = ', real(dof_savg(1))
    end if
  end if

! Variance spectrum --> PSD [ (var)/(cycle/unit_time) ] ----------------
  psd = psd * 2.*dt
  ps_theor = ps_theor * 2.*dt
  signif = signif * 2.*dt
  sig_glob = sig_glob * 2.*dt

! Output ---------------------------------------------------------------
  x_out(:ndata)  = x(:)
  x_out(ndata+1) = x(ndata) + dt*10.

  if (savgopt .eq. 0) then
    call out2d(trim(ofname)//'.nc',1,(/'PS'/),psd,'FScale',nj,fscale, &
               'Xout',ndata+1,x_out,'Wavelet local and global-avg. PSD')
  else
    call out2d1d(trim(ofname)//'.nc',1,(/'PS'/),psd,'FScale',nj,fscale, &
           'Xout',ndata+1,x_out,1,(/'PS_savg'/),pwr_savg,'X',ndata,x, &
           'Wavelet local, global-avg., and scale-avg. PSD')
  end if

  if (sigopt .eq. 1) then
    if (savgopt .eq. 0) then
      call out1d(trim(ofname)//'sig.nc', &
             3,(/'PS_theor','Sig','Sig_glob'/), &
             (/ps_theor,signif,sig_glob/),'FScale',nj,fscale, &
             'Significance test')
    else
      call out1d2(trim(ofname)//'sig.nc', &
             3,(/'PS_theor','Sig','Sig_glob'/), &
             (/ps_theor,signif,sig_glob/),'FScale',nj,fscale, &
             1,(/'Scale_avg'/),(/mean_savg,sig_savg/), &
             'confi_lvl',2,(/0.5,1.-siglvl/),'Significance test')
    end if
  end if

! Reconstruction -------------------------------------------------------
  if (reconopt .eq. 1) then

  if (sigopt .eq. 0) then
    ymean = 0.d0
    do i=1, ndata
      ymean = ymean + ydata8(i)
    enddo
    ymean = ymean / ndata
      variance = 0.d0
    do i=1, ndata
      variance = variance + (ydata8(i) - ymean)**2
    enddo
    variance = variance/(dble(ndata))
  end if

  var_recon = 0.d0
  do i=1, ndata
  do j=1, nj
    var_recon = var_recon + (abs(wcoef(i,j))**2)/wscale(j)
  enddo
  enddo
  var_recon = var_recon*dj*dt8/(cdelta*ndata)
  print'(a,f19.10)', ' Original mean          =', ymean
  print'(a,f19.10)', ' Original variance      =', variance
  print'(a,f19.10)', ' Reconstructed variance =', var_recon
  print'(a,f14.5,a)', ' Var. Ratio             =', var_recon/variance, &
                        '      (low due to zero-padding)'

  mean_recon = 0.d0
  rms_recon = 0.d0
  yrecon = 0.d0
  do i=1, ndata
    do j=1, nj
      yrecon(i) = yrecon(i) + (dble(wcoef(i,j)))/sqrt(wscale(j))
    enddo
    yrecon(i) = yrecon(i)*dj*sqrt(dt8)/(cdelta*psi0)
    rms_recon = rms_recon + (ydata8(i)-ymean-yrecon(i))**2
    mean_recon = mean_recon + yrecon(i)
  enddo
  mean_recon = mean_recon/ndata
  rms_recon = sqrt(rms_recon/ndata)

  print*, 'RMS difference of reconstructed time series'
  print'(a,f14.5)', '                        =', real(rms_recon)
  print*

  end if

end subroutine wlncpack

subroutine wlpack(ndata,dt,x,ydata,ma,nj,s0,subscale,w_ftn,param,pad0, &
               savgopt,avgscale,sigopt,redopt,siglvl,reconopt,outopt,ofname, &
               fscale,psd,ps_theor,signif,sig_glob, &
               pwr_savg,mean_savg,sig_savg)

  implicit none

  integer, intent(in) :: ndata, ma, nj, subscale, w_ftn, pad0
  integer, intent(in) :: savgopt, sigopt, redopt, reconopt, outopt
  real,    intent(in) :: dt, s0, param, avgscale(2), siglvl
  real, dimension(ndata), intent(in) :: x, ydata
  character(len=*),       intent(in) :: ofname

  real, dimension(nj,ndata+1), intent(out) :: psd
  real, dimension(ndata),      intent(out) :: pwr_savg
  real, dimension(nj), intent(out) :: fscale, ps_theor, signif, sig_glob
  real,                intent(out) :: mean_savg, sig_savg

  integer :: i,j, ndpad, javg1, javg2, ncid, ia(ma)
  double complex, dimension(ndata,nj) :: wcoef
  double precision, dimension(ndata) :: ydata8, yrecon, coi
  double precision, dimension(nj)    :: wscale, fscale8
  double precision, dimension(nj)    :: dof, ps_theor8, signif8, mean_savg8
  double precision, dimension(nj)    :: sig_glob8, dof_savg, sig_savg8
  double precision :: mean_recon, var_recon, rms_recon
  double precision :: pi, lag1, ymean, variance, cdelta, psi0
  double precision :: dt8, dj, param8, siglvl8
  real, dimension(ndata) :: sig_det, func, trend
  real :: a(ma), afunc(50), chisq_det, cov(ma+10,ma+10)
  real :: lag1_4, x_out(ndata+1)


  pi = 4.d0 * atan(1.d0)
  dj = 1.d0/subscale
  ydata8 = dble(ydata)
  dt8 = dble(dt)
  param8 = dble(param)

  ndpad = int(2.d0**(nint(log(dble(ndata))/log(2.d0))+1))
  if (pad0 .eq. 0) then
    ndpad = ndata
  endif


! Remove trend in data -------------------------------------------------
  trend = 0.
  if (ma .gt. 1) then
    ia(:) = 1           ! calculate the coeff.(=a) if ia.ne.0
    a(:)  = 0.0         ! coeff. for ia=0
    sig_det(:) = 1.0
    call lfit(x,ydata,sig_det,ndata,a,ia,ma,cov,ma+10,chisq_det)
    do i=1, ndata
      call funcs(x(i),afunc,ma)
      do j=1, ma
        trend(i) = trend(i) + a(j)*afunc(j)
      end do
    end do
  end if

  ydata8 = ydata8 - trend

  print*, ' Edge of the data : ', real(ydata8(1)), real(ydata8(ndata))

! Wavelet transform & calculate the power spectrum ---------------------
  call wavelet(ndata,ydata8,dt8,w_ftn,param8,dble(s0),dj,nj,ndpad, &
               wcoef,wscale,fscale8,coi)
  fscale = real(fscale8)

  psd = 0.
  do i=1, ndata
    psd(:,i) = real(abs(wcoef(i,:))**2)
  enddo

  ! Global averaged power spectrum
  do i=1, ndata
    psd(:,ndata+1) = psd(:,ndata+1) + psd(:,i)
  enddo
  psd(:,ndata+1) = psd(:,ndata+1) / ndata

  ! define constants (only four empirical values)
  cdelta = -1.d0  ;  psi0 = -1.d0
  if ( (w_ftn.eq.0) .and. (param8.eq.6.d0) ) then
    cdelta = 0.776d0  ;  psi0 = pi**(-0.25d0)
  else if ( (w_ftn.eq.1) .and. (int(param8).eq.4) ) then
    cdelta = 1.132d0  ;  psi0 = 1.079d0
  else if ( (w_ftn.eq.2) .and. (int(param8).eq.2) ) then
    cdelta = 3.541d0  ;  psi0 = 0.867d0
  else if ( (w_ftn.eq.2) .and. (int(param8).eq.6) ) then
    cdelta = 1.966d0  ;  psi0 = 0.884d0
  end if

! Scale-averaged power series ------------------------------------------
  pwr_savg = 0.
  if (savgopt .eq. 1) then

    dof_savg(1) = dble(avgscale(1))                ! smallest scale
    dof_savg(2) = dble(avgscale(2))                ! largest scale
    javg1 = 0  ;  javg2 = 0
    do j=1, nj
      if ( (fscale8(j).ge.dof_savg(1)) .and. (javg1.eq.0) )  javg1 = j
      if (fscale8(j) .le. dof_savg(2))  javg2 = j
    enddo
    do j=javg1, javg2
      pwr_savg(:) = pwr_savg(:) + psd(j,:)/real(wscale(j))
    enddo
    pwr_savg(:) = pwr_savg(:)*real(dj*dt8/cdelta)

  end if

! Significant test -----------------------------------------------------
  lag1 = 0.d0
  if (redopt .eq. 1) then
    call autocor1(ndata,real(ydata8),lag1_4)
  end if
  lag1 = dble(lag1_4)

  if (sigopt .eq. 1) then
    siglvl8 = dble(siglvl)
    ! regular chi-square test
    call wave_signif(0,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof,ps_theor8,signif8, &
            ymean,variance,cdelta,psi0)

    ! significance test of the global averaged spectrum
    dof(:) = ndata - wscale(:)
    call wave_signif(1,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof,ps_theor8,sig_glob8, &
            ymean,variance,cdelta,psi0)

    ps_theor = real(ps_theor8)
    signif   = real(signif8)
    sig_glob = real(sig_glob8)

    ! significance test of the scale-averaged time series
    if (savgopt .eq. 1) then
      call wave_signif(2,ndata,ydata8,dt8,w_ftn,param8,dj,nj, &
            wscale,fscale8,lag1,siglvl8,dof_savg,mean_savg8,sig_savg8, &
            ymean,variance,cdelta,psi0)

      mean_savg = real(mean_savg8(1))
      sig_savg = real(sig_savg8(1))
      print*, ' Scale-avg. Degrees of freedom = ', real(dof_savg(1))
    end if
  end if

! Variance spectrum --> PSD [ (var)/(cycle/unit_time) ] ----------------
  psd = psd * 2.*dt
  ps_theor = ps_theor * 2.*dt
  signif = signif * 2.*dt
  sig_glob = sig_glob * 2.*dt

! Output ---------------------------------------------------------------
  if (outopt .eq. 1) then
    x_out(:ndata)  = x(:)
    x_out(ndata+1) = x(ndata) + dt*10.

    if (savgopt .eq. 0) then
      call out2d(trim(ofname)//'.nc',1,(/'PS'/),psd,'FScale',nj,fscale, &
                 'Xout',ndata+1,x_out,'Wavelet local and global-avg. PSD')
    else
      call out2d1d(trim(ofname)//'.nc',1,(/'PS'/),psd,'FScale',nj,fscale, &
             'Xout',ndata+1,x_out,1,(/'PS_savg'/),pwr_savg,'X',ndata,x, &
             'Wavelet local, global-avg., and scale-avg. PSD')
    end if

    if (sigopt .eq. 1) then
      if (savgopt .eq. 0) then
        call out1d(trim(ofname)//'sig.nc', &
               3,(/'PS_theor','Sig','Sig_glob'/), &
               (/ps_theor,signif,sig_glob/),'FScale',nj,fscale, &
               'Significance test')
      else
        call out1d2(trim(ofname)//'sig.nc', &
               3,(/'PS_theor','Sig','Sig_glob'/), &
               (/ps_theor,signif,sig_glob/),'FScale',nj,fscale, &
               1,(/'Scale_avg'/),(/mean_savg,sig_savg/), &
               'confi_lvl',2,(/0.5,1.-siglvl/),'Significance test')
      end if
    end if

  end if

! Reconstruction -------------------------------------------------------
  if (reconopt .eq. 1) then

  if (sigopt .eq. 0) then
    ymean = 0.d0
    do i=1, ndata
      ymean = ymean + ydata8(i)
    enddo
    ymean = ymean / ndata
      variance = 0.d0
    do i=1, ndata
      variance = variance + (ydata8(i) - ymean)**2
    enddo
    variance = variance/(dble(ndata))
  end if

  var_recon = 0.d0
  do i=1, ndata
  do j=1, nj
    var_recon = var_recon + (abs(wcoef(i,j))**2)/wscale(j)
  enddo
  enddo
  var_recon = var_recon*dj*dt8/(cdelta*ndata)
  print'(a,f19.10)', ' Original mean          =', ymean
  print'(a,f19.10)', ' Original variance      =', variance
  print'(a,f19.10)', ' Reconstructed variance =', var_recon
  print'(a,f14.5,a)', ' Var. Ratio             =', var_recon/variance, &
                        '      (low due to zero-padding)'

  mean_recon = 0.d0
  rms_recon = 0.d0
  yrecon = 0.d0
  do i=1, ndata
    do j=1, nj
      yrecon(i) = yrecon(i) + (dble(wcoef(i,j)))/sqrt(wscale(j))
    enddo
    yrecon(i) = yrecon(i)*dj*sqrt(dt8)/(cdelta*psi0)
    rms_recon = rms_recon + (ydata8(i)-ymean-yrecon(i))**2
    mean_recon = mean_recon + yrecon(i)
  enddo
  mean_recon = mean_recon/ndata
  rms_recon = sqrt(rms_recon/ndata)

  print*, 'RMS difference of reconstructed time series'
  print'(a,f14.5)', '                        =', real(rms_recon)
  print*

  end if


  return

end subroutine wlpack

END module waveletpack


subroutine autocor1(n,x,cor1)

  implicit none

  integer, intent(in) :: n
  real, intent(in)    :: x(n)
  real, intent(out)   :: cor1
  integer :: i
  real    :: x0(n-1), x1(n-1), x0mean, x1mean
  real    :: tmp, tmp2

  do i=1, n-1
    x0(i) = x(i)
    x1(i) = x(i+1)
  enddo
  tmp = 0.
  do i=1, n
    tmp = tmp + x(i)
  enddo
  x0mean = (tmp-x(n)) / (n-1)
  x1mean = (tmp-x(1)) / (n-1)
  x0 = x0 - x0mean
  x1 = x1 - x1mean
  cor1 = 0.
  tmp = 0.  ;  tmp2 = 0.
  do i=1, n-1
    cor1 = cor1 + x0(i)*x1(i)
    tmp  = tmp  + x0(i)**2
    tmp2 = tmp2 + x1(i)**2
  enddo
  cor1 = cor1 / sqrt(tmp*tmp2)

  return

end subroutine autocor1

