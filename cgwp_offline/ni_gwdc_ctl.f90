




! Subroutine NI_gwdc_ctl

      Subroutine NI_gwdc_ctl (                                          &

! parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, at_extremity, n_proc, n_procx, n_procy                          &
     &, neighbour, g_datastart, me                                      &

! model dimensions
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_levels                                        &

! model switches
     &, model_domain                                                    &
     &, Ltimer, l_gwdc                                                  &

! IN coordinate information
     &, eta_theta_levels, r_theta_levels                                &
     &, sec_theta_latitude                                              &

! IN time stepping information
     &, sub_timestep, Substep_Number, Num_Substeps                      &
     &, timestep_number                                                 &

! diagnostic info
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &  STASHwork                                                       &

! IN data fields
     &, u, v                                                            &
     &, p_layer_boundaries                                              &
     &, rho_wet, rho_wet_theta, theta                                   &
     &, exner_theta_levels                                              &
     &, scheat, kscbas, ksctop                                          &
     &, cumulus, l_shallow, l_congestus                                 &

! IN/OUT
     &, R_u, R_v                                                        &

! error information
     &, Error_code  )

!
! PURPOSE:  Interface to Atmospheric Physics GWDC Scheme.
!           1) Define several parameters used in the GWDC scheme.
!           2) Call subroutines of GWDC scheme and its diagnostics.
!
!
! CURRENT CODE OWNER:  Young-Ha Kim
!
!
! HISTORY:
!
!     Date    Comment
!   --------  -------
!   11/05/09  First written for the UM.                        Y.-H. Kim
!
!
! CODE DESCRIPTION:
!
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with intent IN. ie: input variables.

! Parallel setup variables

      Integer ::                                                        &
     &  halo_i                                                          &
                             ! size of halo in i direction
     &, halo_j                                                          &
                             ! size of halo in j direction
     &, off_x                                                           &
                             ! size of small halo in i direction
     &, off_y                                                           &
                             ! size of small halo in j direction
     &, global_row_length                                               &
                             ! number of points on a row
     &, global_rows                                                     &
                             ! number of global rows
     &, n_proc                                                          &
                             ! total number of processors
     &, n_procx                                                         &
                             ! number of processors in longitude
     &, n_procy                                                         &
                             ! number of processors in latitude
     &, neighbour(4)                                                    &
                             ! array with the IDs of the four neighbours
                             ! in the horizontal plane
     &, g_datastart (3,0:n_proc-1)                                      &
                             ! starting indices in 3-D global domain for
                             ! each processor
     &, me                   ! my processor number

      Logical ::                                                        &
     &  at_extremity(4)  ! indicates if this processor is at north, 
                         ! south, east or west of the processor grid

! Model dimensions

      Integer ::                                                        &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, wet_levels

! Model switches

      Integer ::                                                        &
     &  model_domain

!yh+2
!      Logical ::                                                        &
!     &  Ltimer                                                          &
!                      ! true then output some timing information
!     &, L_gwdc        ! switch for the convective GWD scheme

      Logical ::                                                        &
     &  Ltimer
                      ! true then output some timing information
      Real ::                                                           &
     &  L_gwdc(6)     ! switch for the convective GWD scheme
!yh-2

! Model parameters

      Real ::                                                           &
     &  timestep                                                        &
     &, sub_timestep

! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! TYPSTS starts
! CSUBMODL must be included before this file
!Applicable to all configurations (except MOS variables)
!STASH related variables for describing output requests and space
!management.
!LL
!LL   AUTHOR            Rick Rawlins
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2             Code creation for Dynamic allocation
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.5  Apr. 95   Sub-Models project.
!LL                  Dimensioning of various STASH arrays altered in
!LL                  accordance with internal model separation scheme.
!LL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
!LL                  longer required.
!LL                  S.J.Swarbrick
!LL
!
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS

      INTEGER :: MOS_MASK_LEN         ! Size of bit mask for MOS

      COMMON/DSIZE_AO/  MOS_MASK_LEN

! TYPSTSZ end
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface
! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      INTEGER :: MOS_OUTPUT_LENGTH
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP,MOS_OUTPUT_LENGTH

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
      INTEGER:: MOS_MASK(MOS_MASK_LEN)
! TYPSTS end

! Diagnostics info
       Real ::                                                          &
     & STASHwork(*) ! STASH workspace


! Data arrays

      Real ::                                                           &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, rho_wet(row_length, rows, model_levels)                         &
     &, rho_wet_theta(row_length, rows, model_levels-1)                 &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        model_levels)                                             &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &                     model_levels)                                &
     &, p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, scheat(row_length, rows, wet_levels)                            &
              ! subgrid-scale convective heating rate [K/s]
     &, schmax(row_length, rows)
              ! maximum subgrid-scale convective heating rate [K/s]

      Integer ::                                                        &
     &  kscbas(row_length, rows)                                        &
     &, ksctop(row_length, rows)

      Logical ::                                                        &
     &  cumulus    (row_length, rows)                                   &
     &, l_shallow  (row_length, rows)                                   &
     &, l_congestus(row_length, rows)

! Co-ordinate arrays
      Real ::                                                           &
     &  eta_theta_levels(0:model_levels)                                &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y)

! Time information for current timestep
      Integer ::                                                        &
     &  timestep_number                                                 &
     &, Substep_Number                                                  &
     &, Num_Substeps

! Arguments with intent in/out. ie: input variables changed on output.

      Real ::                                                           &
     &  R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

      Integer ::                                                        &
     &  Error_code

! Global variables

!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------

! Local variables

      Logical ::                                                        &
        ! switches for STASH
     &  gwdc_drag_u_on                                                  &
     &, gwdc_drag_v_on                                                  &
     &, gwdc_mom_east_on                                                &
     &, gwdc_mom_west_on                                                &
     &, gwdc_mom_north_on                                               &
     &, gwdc_mom_south_on                                               &
     &, gwdc_mom_e_ctop_on                                              &
     &, gwdc_mom_w_ctop_on                                              &
     &, gwdc_mom_n_ctop_on                                              &
     &, gwdc_mom_s_ctop_on                                              &
     &, gwdc_heatmax_on                                                 &
     &, gwdc_znwcq_on                                                   &
     &, gwdc_spec_on

      Integer ::                                                        &
     &  i,j,k, ij     ! loop counters

      Integer, parameter ::                                             &
     &  t_shallow       = 1                                             &
     &, t_congestus     = 2                                             &
     &, t_deep          = 3

      Integer, parameter ::                                             &
     &  ctype_c = t_deep       ! Criterion of convection type
                               ! for parameterizing GWDC
                               !   1, all subgrid-scale convection
                               !   2, congestus or deep convection
                               !   3, deep convection only
                               ! 3 is recommended (for computation time)

      Integer, parameter ::                                             &
     &  PNorth = 1                                                      &
     &, PSouth = 3

! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

! Local data arrays

      Logical ::                                                        &
     &  gwdc_col(row_length, rows) ! columns at which GWDC is calculated

      Integer ::                                                        &
     &  nc_gwdc                                                         &
                                     ! number of GWDC columns
     &, type_conv(row_length, rows)  ! convection type

      Integer, Dimension(:), Allocatable ::                             &
     &  igwdc                                                           &
                      ! i-direction index for 1-D array of GWDC columns
     &, jgwdc         ! j-direction index for 1-D array of GWDC columns

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostics not requested
      Real, Dimension(:,:,:), Allocatable ::                            &
     &  gwdc_drag_u                                                     &
     &, gwdc_drag_v                                                     &
     &, gwdc_mom_east                                                   &
     &, gwdc_mom_west                                                   &
     &, gwdc_mom_north                                                  &
     &, gwdc_mom_south                                                  &
     &, gwdc_znwcq

      Real, Dimension(:,:), Allocatable ::                              &
     &  gwdc_mom_e_ctop                                                 &
     &, gwdc_mom_w_ctop                                                 &
     &, gwdc_mom_n_ctop                                                 &
     &, gwdc_mom_s_ctop                                                 &
     &, gwdc_heatmax

      Real, Dimension(:,:,:,:), Allocatable ::                          &
     &  gwdc_spec

! External routines

      External timer
      External gw_conv, diagnostics_gwdc


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

      ! parameters for Warner and McIntyre scheme (see UMDP34.)
      ! The code of gw_conv is written assuming (p,s,t) as (5/3,1,3).
      real, parameter ::  mstar   = 2.0*pi/4.3e3
      real, parameter ::  p_wm    = 5.0/3.0
      real, parameter ::  s_wm    = 1.0
      real, parameter ::  t_wm    = 3.0
!yh+2
!      real, parameter ::  beta_wm = 2.0
!
!      ! Criterion of the maximum subgrid-scale heating rate [K/s]
!      ! normalized by cfactor
!      real, parameter ::  schm_c = 0.0001/86400. * (125./cfactor)
!yh-2
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!yh+2
      real ::  cfactor, beta_wm, schm_c

      cfactor = l_gwdc(2)
      beta_wm = l_gwdc(3)
      schm_c = 0.0001/86400. * (125./cfactor)
      write(6,*) 'cfactor, beta_wm : ', cfactor, beta_wm
!yh-2

! ----------------------------------------------------------------------
! Section GWDC.1  Set stash diagnostic switches and
!                     columns at which GWDC will be calculated.
! ----------------------------------------------------------------------

! STASH DIAGNOSTIC SWITCHES

      gwdc_drag_u_on     = sf(501,6)
      gwdc_drag_v_on     = sf(502,6)
      gwdc_mom_east_on   = sf(503,6)
      gwdc_mom_west_on   = sf(504,6)
      gwdc_mom_north_on  = sf(505,6)
      gwdc_mom_south_on  = sf(506,6)
      gwdc_mom_e_ctop_on = sf(507,6)
      gwdc_mom_w_ctop_on = sf(508,6)
      gwdc_mom_n_ctop_on = sf(509,6)
      gwdc_mom_s_ctop_on = sf(510,6)
      gwdc_heatmax_on    = sf(511,6)
      gwdc_znwcq_on      = sf(512,6)
      gwdc_spec_on       = sf(301,6)  ! 301-483


      if (error_code  ==  0) then

! DEPENDS ON: timer
      if ( Ltimer )  Call timer ('Convective GW Drag ',3)

      ! Initialize diagnostics
      if ( gwdc_drag_u_on ) then
        allocate( gwdc_drag_u(row_length,rows,model_levels) )
        gwdc_drag_u(:,:,:) = 0.0
      end if

      if ( gwdc_drag_v_on ) then
        allocate( gwdc_drag_v(row_length,n_rows,model_levels) )
        gwdc_drag_v(:,:,:) = 0.0
      end if

      if ( gwdc_mom_east_on ) then
        allocate( gwdc_mom_east(row_length,rows,model_levels) )
        gwdc_mom_east(:,:,:) = 0.0
      end if

      if ( gwdc_mom_west_on ) then
        allocate( gwdc_mom_west(row_length,rows,model_levels) )
        gwdc_mom_west(:,:,:) = 0.0
      end if

      if ( gwdc_mom_north_on ) then
        allocate( gwdc_mom_north(row_length,rows,model_levels) )
        gwdc_mom_north(:,:,:) = 0.0      
      end if

      if ( gwdc_mom_south_on ) then
        allocate( gwdc_mom_south(row_length,rows,model_levels) )
        gwdc_mom_south(:,:,:) = 0.0      
      end if

      if ( gwdc_mom_e_ctop_on ) then
        allocate( gwdc_mom_e_ctop(row_length,rows) )
        gwdc_mom_e_ctop(:,:) = 0.0      
      end if

      if ( gwdc_mom_w_ctop_on ) then
        allocate( gwdc_mom_w_ctop(row_length,rows) )
        gwdc_mom_w_ctop(:,:) = 0.0      
      end if

      if ( gwdc_mom_n_ctop_on ) then
        allocate( gwdc_mom_n_ctop(row_length,rows) )
        gwdc_mom_n_ctop(:,:) = 0.0      
      end if

      if ( gwdc_mom_s_ctop_on ) then
        allocate( gwdc_mom_s_ctop(row_length,rows) )
        gwdc_mom_s_ctop(:,:) = 0.0      
      end if

      if ( gwdc_heatmax_on ) then
        allocate( gwdc_heatmax(row_length,rows) )
        gwdc_heatmax(:,:) = 0.0      
      end if

      if ( gwdc_znwcq_on ) then
        allocate( gwdc_znwcq(row_length,rows,model_levels) )
        gwdc_znwcq(:,:,:) = -999.
      end if

      if ( gwdc_spec_on ) then
        allocate( gwdc_spec(row_length,rows,model_levels,(nc*2+1)*3) )
        gwdc_spec(:,:,:,:) = 0.0
      end if


! COLUMNS AT WHICH GWDC WILL BE PARAMETERIZED

      ! maximum subgrid-scale convective heating rate
      do j=1, rows
      do i=1, row_length
        schmax(i,j) = maxval( scheat(i,j,kscbas(i,j):ksctop(i,j)) )
      enddo
      enddo

      ! convection type
      type_conv(:,:) = 0
      do j=1, rows
      do i=1, row_length
!        if ( schmax(i,j) > 0.0 .and. ksctop(i,j) > 0 ) then
        if ( cumulus(i,j) ) then
          if ( l_congestus(i,j) ) then
            type_conv(i,j) = t_congestus    ! congestus conv.
          else if ( l_shallow(i,j) ) then
            type_conv(i,j) = t_shallow      ! shallow conv.
          else
            type_conv(i,j) = t_deep         ! deep conv.
          end if
        end if
      enddo
      enddo

      if (model_domain == mt_global) then
        if ( at_extremity(PSouth) ) then
          type_conv(:,1) = 0
          schmax   (:,1) = 0.0
        end if
        if ( at_extremity(PNorth) ) then
          type_conv(:,rows) = 0
          schmax   (:,rows) = 0.0
        end if
      end if

      ! columns for calculating GWDC
      nc_gwdc = 0
      do j=1, rows
      do i=1, row_length
        if ( type_conv(i,j) >= ctype_c .and. schmax(i,j) > schm_c ) then
          gwdc_col(i,j) = .TRUE.
          nc_gwdc = nc_gwdc + 1
        else
          gwdc_col(i,j) = .FALSE.
        end if
      enddo
      enddo

      ! write if needed
      write(6,*) ' [GWDC] number of columns : ', nc_gwdc

      allocate( igwdc(nc_gwdc), jgwdc(nc_gwdc) )

      ! 2-D indices -> 1-D index
      ij = 0
      do j=1, rows
      do i=1, row_length
        if ( gwdc_col(i,j) ) then
          ij = ij + 1
          igwdc(ij) = i
          jgwdc(ij) = j
        end if
      enddo
      enddo

! ----------------------------------------------------------------------
! Section GWDC.2  Call convective gravity wave drag scheme
! ----------------------------------------------------------------------

! DEPENDS ON: gw_conv
      Call gw_conv(row_length,rows,n_rows,off_x,off_y                   &
     &, halo_i,halo_j,n_proc,g_datastart,me                             &
     &, model_domain,at_extremity                                       &
     &, model_levels,wet_levels                                         &
     &, nc_gwdc,nphi,nc,c_max,phi_dir,hscale,tscale,cfactor,lt,ah       &
     &, mstar,p_wm,s_wm,t_wm,beta_wm                                    &
     &, u,v,theta,exner_theta_levels,rho_wet,rho_wet_theta              &
     &, eta_theta_levels,r_theta_levels,sec_theta_latitude              &
     &, sub_timestep, Substep_Number, Num_Substeps, timestep_number     &
     &, scheat,schmax,kscbas,ksctop,igwdc,jgwdc                         &
     &, R_u, R_v,l_gwdc                                                 &
! diagnostics
     &, gwdc_drag_u,gwdc_drag_v,gwdc_drag_u_on,gwdc_drag_v_on           &
     &, gwdc_mom_east,gwdc_mom_west                                     &
     &, gwdc_mom_north,gwdc_mom_south                                   &
     &, gwdc_mom_east_on,gwdc_mom_west_on                               &
     &, gwdc_mom_north_on,gwdc_mom_south_on                             &
     &, gwdc_mom_e_ctop,gwdc_mom_w_ctop                                 &
     &, gwdc_mom_n_ctop,gwdc_mom_s_ctop                                 &
     &, gwdc_mom_e_ctop_on,gwdc_mom_w_ctop_on                           &
     &, gwdc_mom_n_ctop_on,gwdc_mom_s_ctop_on                           &
     &, gwdc_heatmax,gwdc_heatmax_on                                    &
     &, gwdc_znwcq,gwdc_znwcq_on                                        &
     &, gwdc_spec,gwdc_spec_on                                          &
     &  )

      deallocate( igwdc, jgwdc )

! DEPENDS ON: timer
      if ( Ltimer )  Call timer ('Convective GW Drag ',4)

! ----------------------------------------------------------------------
! Section GWDC.3  Call GWDC diagnostics
! ----------------------------------------------------------------------

      if ( sf(0,6) ) then  ! diagnostics requested this timestep

! DEPENDS ON: timer
        if ( Ltimer )  Call timer ('Diags   ',3)

! DEPENDS ON: diagnostics_gwdc
        Call diagnostics_gwdc(                                          &
     &                        row_length, rows, model_levels, n_rows    &
     &,                       off_x, off_y, at_extremity                &
     &,                       gwdc_drag_u, gwdc_drag_v                  &
     &,                       gwdc_mom_east, gwdc_mom_west              &
     &,                       gwdc_mom_north, gwdc_mom_south            &
     &,                       gwdc_mom_e_ctop, gwdc_mom_w_ctop          &
     &,                       gwdc_mom_n_ctop, gwdc_mom_s_ctop          &
     &,                       gwdc_heatmax, gwdc_znwcq,                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &                        STASHwork                                 &
     &                       )

        ! Deallocate diagnostic variables
        if ( gwdc_drag_u_on )  deallocate( gwdc_drag_u )
        if ( gwdc_drag_v_on )  deallocate( gwdc_drag_v )

        if ( gwdc_mom_east_on  )  deallocate( gwdc_mom_east  )
        if ( gwdc_mom_west_on  )  deallocate( gwdc_mom_west  )
        if ( gwdc_mom_north_on )  deallocate( gwdc_mom_north )
        if ( gwdc_mom_south_on )  deallocate( gwdc_mom_south )

        if ( gwdc_mom_e_ctop_on )  deallocate( gwdc_mom_e_ctop )
        if ( gwdc_mom_w_ctop_on )  deallocate( gwdc_mom_w_ctop )
        if ( gwdc_mom_n_ctop_on )  deallocate( gwdc_mom_n_ctop )
        if ( gwdc_mom_s_ctop_on )  deallocate( gwdc_mom_s_ctop )

        if ( gwdc_heatmax_on )  deallocate( gwdc_heatmax )
        if ( gwdc_znwcq_on   )  deallocate( gwdc_znwcq   )

        if ( gwdc_spec_on ) then
! DEPENDS ON: diagnostics_gwdc2
          Call diagnostics_gwdc2(                                       &
     &                           row_length, rows, model_levels, n_rows &
     &,                          off_x, off_y, at_extremity             &
     &,                          nc, gwdc_spec,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &                           STASHwork                              &
     &                          )
          deallocate( gwdc_spec )
        end if

! DEPENDS ON: timer
        if ( Ltimer )  Call timer ('Diags   ',4)

      end if  ! on sf(0,6)

      end if  ! on error_code == 0


      Return
      END SUBROUTINE NI_gwdc_ctl

