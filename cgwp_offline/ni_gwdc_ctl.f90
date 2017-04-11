
! GWDC parameters ::::::::::::::::::::::::::::::::::::::::::::::::::::::
      integer,               parameter ::  nphi    = 2
      integer,               parameter ::  nc      = 30
      real,                  parameter ::  c_max   = 60.
!yh+2
!      real,                  parameter ::  cfactor = 200.
!yh-2
      real, dimension(nphi), parameter ::  phi_dir = (/45.,135./)
!
!     nphi    - Number of wave-propagation directions considered
!
!     nc      - Number of positive phase speeds in the discrete spectrum
!               Total number is nc*2+1 (including zero and negatives).
!
!     c_max   - Maximum phase speed in the spectrum
!
!     cfactor - Conversion factor for magnitude of the cloud-top
!               momentum flux
!               (See the appendix in Song et al. (2007, JAS).)
!
!     phi_dir - Wave-propagation directions considered
!               They must be in [0,180).
!               45, 135 (deg) are chosen by Choi and Chun (2011, JAS).
!

      ! parameters in Song and Chun (2005, JAS)
      ! (see Song et al., 2007, JAS.)
      real, parameter ::  hscale = 5.e3
      real, parameter ::  tscale = 1200.
      real, parameter ::  lt     = tscale
      real, parameter ::  ah     = pi*1.e4*tscale*tscale

      ! Criterion of the maximum subgrid-scale heating rate [K/s]
      real, parameter ::  schm_c = 0.0001/86400.
!yh-2
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!yh+2
      real ::  cfactor, beta_wm

      cfactor = 125.0
!yh-2

