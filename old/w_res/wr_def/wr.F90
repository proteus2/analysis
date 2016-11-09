    program w_definition

    use temeq
    use avg

    implicit none

    include 'netcdf.inc'
    include 'c_math.inc'

    integer, parameter ::  nx = 144, ny = 73, nz = 23
    integer, parameter ::  mon1 = 6, mon2 = 8
    integer, parameter ::  year1 = 1994, year2 = 1995
!------------------------------------------------------------------------------------
!   Grids and data needed for calculation
!------------------------------------------------------------------------------------

    real, dimension(ny)        ::  lat
    real, dimension(nz)        ::  lev

    real, dimension(:,:,:,:), allocatable ::  u, v, t, ome
    real, dimension(:,:,:),   allocatable ::  epfy, epfz, epdrag

    real, dimension(:,:,:),   allocatable ::  vmrd, wmr, vadv_u, wadv_u, cor
    real, dimension(:,:,:),   allocatable ::  vadv_t,wadv_t, wmrd, wmrm
    real, dimension(:,:,:),   allocatable ::  wmd, wmm, wm
    real, dimension(:,:),     allocatable ::  wmr1, wm1

    integer ::  nday, year, month
    integer ::  n, imon, iyr, nn, nmon, nyear
    real    ::  mon(12)


    uniform_y = .TRUE.
    uniform_z = .FALSE.
    l_ybdy = (/.FALSE.,.FALSE./)
    l_zbdy = (/.TRUE. ,.TRUE. /)

    mon = (/1,2,3,4,5,6,7,8,9,10,11,12/)*1.

    nmon  = mon2  - mon1  + 1
    nyear = year2 - year1 + 1

    allocate( wmr1(ny,nz),       wm1(ny,nz)       )
    allocate( wmr (ny,nz,nmon ), wm (ny,nz,nmon ) )
    allocate( wmrm(ny,nz,nyear), wmm(ny,nz,nyear) )

    nn = 0

    N_MON:  DO imon=mon1, mon2
    N_YR:   DO iyr=year1, year2


    nn = nn + 1

    year = iyr
    month = imon
    print*, year, month

print*, 'read'
    call readdata

    allocate( wmrd(ny,nz,nday), vmrd(ny,nz,nday), wmd(ny,nz,nday) )

print*, 'calculation'
    call utend_tem_hydp(nx,ny,nz,nday,lat,lev*100.,u,v,ome,t,1, &
                        vmrd,wmrd,vadv_u,wadv_u,cor,epdrag,epfy,epfz, &
                        vadv_t,wadv_t,wmd)

    deallocate( u, v, t, ome )

    call avg_d3((/ny,nz,nday,1/),wmrd,1., wmrm(:,:,iyr-year1+1))
    call avg_d3((/ny,nz,nday,1/),wmd ,1., wmm (:,:,iyr-year1+1))

    deallocate( wmrd, vmrd, wmd )

    ENDDO  N_YR

    call avg_d3((/ny,nz,nyear,1/),wmrm,1., wmr(:,:,imon-mon1+1))
    call avg_d3((/ny,nz,nyear,1/),wmm ,1., wm (:,:,imon-mon1+1))

    ENDDO  N_MON

    call avg_d3((/ny,nz,nmon,1/),wmr,1., wmr1)
    call avg_d3((/ny,nz,nmon,1/),wm ,1., wm1 )

    call dump


    CONTAINS

!------------------------------------------------------------------------------

    subroutine readdata 
  
    character(len=256) :: rfn1

    integer :: i,n
    integer :: istat, ncid1
    integer :: latid, levid, dimid
    integer :: varid(4)
    integer :: lev4(nz)
    real, dimension(:,:,:), allocatable :: tmp3
    integer*2, dimension(:,:,:,:), allocatable ::  tmp
    real*8, dimension(4) :: sf, ao


    write(rfn1,'(a,i4.4,a,i4.4,a,i2.2,a)')   &
       '/data6/kyh/ERA40/6hr/',year,'/'// &
       'ERA40-',year,'.',month,'.nc'

    istat = nf_open(trim(rfn1),nf_nowrite,ncid1)

    istat = nf_inq_varid(ncid1,'latitude',latid)
    istat = nf_inq_varid(ncid1,'levelist',levid)
    istat = nf_inq_varid(ncid1,'u'    ,varid(1))
    istat = nf_inq_varid(ncid1,'v'    ,varid(2))
    istat = nf_inq_varid(ncid1,'t'    ,varid(3))
    istat = nf_inq_varid(ncid1,'w'    ,varid(4))

    istat = nf_inq_dimid(ncid1,'time',dimid)
    istat = nf_inq_dimlen(ncid1,dimid, nday)
    nday = nday / 4

    do i=1, 4 
      istat = nf_get_att_double(ncid1,varid(i),"scale_factor",sf(i))
      istat = nf_get_att_double(ncid1,varid(i),"add_offset"  ,ao(i))
    enddo

    istat = nf_get_var_real(ncid1,latid, lat)
    istat = nf_get_var_int (ncid1,levid, lev4)
    lev = float(lev4)

    allocate( u(nx,ny,nz,nday), v  (nx,ny,nz,nday) )
    allocate( t(nx,ny,nz,nday), ome(nx,ny,nz,nday) )

    allocate( tmp(nx,ny,nz,nday*4) )
    allocate( tmp3(nx,ny,nz) )
    do i=1, 4
      istat = nf_get_var_int2(ncid1,varid(i), tmp)
      do n=1, nday
        tmp3(:,:,:) = (float(tmp(:,:,:,n*4-3))+float(tmp(:,:,:,n*4-2))+ &
                       float(tmp(:,:,:,n*4-1))+float(tmp(:,:,:,n*4  ))) &
                      /4.*real(sf(i))+real(ao(i))
        select case ( i )
          case ( 1 )
            u(:,:,:,n) = tmp3(:,:,:)
          case ( 2 )
            v(:,:,:,n) = tmp3(:,:,:)
          case ( 3 )
            t(:,:,:,n) = tmp3(:,:,:)
          case ( 4 )
            ome(:,:,:,n) = tmp3(:,:,:)
        end select
      enddo
    enddo
    deallocate( tmp3 )
    deallocate( tmp )

    istat = nf_close(ncid1)

    return
    end subroutine readdata

!-----------------------------------------------------------------------

    subroutine dump

    character(len=100) :: wfn1
  
    integer :: istat ,ncid
    integer :: latdid,zdid  ,daydid
    integer :: latid ,zid   ,dayid
    integer :: epdid1, epdid2

    write(wfn1,'(A,i4.4,a,i2.2,a)') 'res/era_wdef.nc'

    istat = nf_create(trim(wfn1),nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',ny  ,latdid)
    istat = nf_def_dim(ncid,'p'  ,nz  ,zdid  )
    istat = nf_def_var(ncid,'lat'    ,nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'p'      ,nf_real,1,zdid  ,zid  )
    istat = nf_def_var(ncid,'w_res'  ,nf_real,2,(/latdid,zdid/),epdid1)
    istat = nf_def_var(ncid,'w'      ,nf_real,2,(/latdid,zdid/),epdid2)
    istat = nf_put_att_real(ncid,epdid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,epdid2,'_FillValue',nf_real,1,1.e+20)
    istat = nf_enddef(ncid)
    istat = nf_put_var_real(ncid,latid ,lat)
    istat = nf_put_var_real(ncid,zid   ,lev)
    istat = nf_put_var_real(ncid,epdid1,wmr1)
    istat = nf_put_var_real(ncid,epdid2,wm1 )
    istat = nf_close(ncid)

    return
    end subroutine dump

    end program w_definition
