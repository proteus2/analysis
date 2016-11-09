    program w_randel

    use temeq
    use avg

    implicit none

    include 'netcdf.inc'
    include 'c_math.inc'

    integer, parameter ::  nx = 144, ny = 73
    integer, parameter ::  nz = 23
    integer, parameter ::  nyo = 62
    integer, parameter ::  mon1 = 6, mon2 = 8
    integer, parameter ::  year1 = 1994, year2 = 1995
    real,    parameter ::  dt = 86400.
    integer            ::  ilat_avg(2)
    data ilat_avg/31,43/

!    integer ::  nyo = ny - (ilat_avg(2)-ilat_avg(1)-1)
    real, dimension(nyo) ::  lato
    real, dimension(ny)  ::  lat
    real, dimension(nz)  ::  lev

    real, dimension(:,:,:), allocatable ::  um
    real, dimension(:,:,:), allocatable ::  epdrag, waa_all

    real, dimension(:,:,:,:), allocatable ::  w_all
    real, dimension(:,:,:,:), allocatable ::  wa_all

    integer ::  nday, year, month, nyear, nmon
    integer ::  n, imon, iyr, nn
    real    ::  mon(mon2-mon1+1)


    uniform_y = .FALSE.
    l_ybdy = (/.FALSE.,.FALSE./)

    nyear = year2-year1+1
    nmon  = mon2-mon1+1

    do n=1, nmon
      mon(n) = mon1+n-1.
    enddo

    allocate( waa_all(nyo,nz,3) )
    allocate( wa_all(nyo,nz,nmon,3) )
    allocate( w_all(nyo,nz,3,nyear) )

    N_MON:  DO imon=mon1, mon2

    nn = 0

    N_YR:   DO iyr=year1, year2

    nn = nn + 1

    year = iyr
    month = imon
    print*, year, month

print*, 'reading'
    call readdata

print*, 'calculation'
    call wmr_dynbal_hydp(ny,nz,nday,dt,lat,lev*100.,um,1,(/epdrag/), &
                         ilat_avg, lato,nyo,w_all(:,:,:,nn))

    deallocate( um )
    deallocate( epdrag )

    ENDDO  N_YR

print*, 'averaging'
    call avg_d4((/nyo,nz,3,nyear/),w_all,1., wa_all(:,:,imon-mon1+1,:))

    ENDDO  N_MON

    deallocate( w_all )

    call avg_d3((/nyo,nz,nmon,3/),wa_all,1., waa_all)

    call dump


    CONTAINS

!-----------------------------------------------------------------------

    subroutine readdata 
 
    use avg
 
    character(len=256) :: rfn1, rfn2

    integer :: j     ,k     ,n
    integer :: istat ,ncid1 ,ncid2
    integer :: latid ,levid ,dimid
    integer :: uid   ,epdid
    integer :: lev4(nz)
    integer*2, dimension(:,:,:,:), allocatable ::  tmp
    real*8  :: sf, ao


    write(rfn1,'(a,i4.4,a,i4.4,a,i2.2,a)')   &
       '/data6/kyh/ERA40/6hr/',year,'/'// &
       'ERA40-',year,'.',month,'.nc'

    write(rfn2,'(a,i4.4,a,i2.2,a)')      &
       '/data3/kyh/analy/epflux/res/'//  &
       'epflx-ERA40-.',year,'-',month,'.nc'

    istat = nf_open(trim(rfn1),nf_nowrite,ncid1)
    istat = nf_open(trim(rfn2),nf_nowrite,ncid2)

    istat = nf_inq_varid(ncid1,'latitude',latid)
    istat = nf_inq_varid(ncid1,'levelist',levid)
    istat = nf_inq_varid(ncid1,'u'    ,uid  )
    istat = nf_inq_varid(ncid2,'epd'     ,epdid)

    istat = nf_inq_dimid(ncid1,'time',dimid)
    istat = nf_inq_dimlen(ncid1,dimid, nday)
    nday = nday / 4

    istat = nf_get_var_real(ncid1,latid, lat)
    istat = nf_get_var_int (ncid1,levid, lev4)
    lev = float(lev4)

    istat = nf_get_att_double(ncid1,uid,"scale_factor",sf)
    istat = nf_get_att_double(ncid1,uid,"add_offset"  ,ao)

    allocate( um(ny,nz,nday), epdrag(ny,nz,nday) )

    allocate( tmp(nx,ny,nz,nday*4) )
    istat = nf_get_var_int2(ncid1,uid  , tmp)
    do n=1, nday
    do k=1, nz
    do j=1, ny
      um(j,k,n) = sum( (float(tmp(:,j,k,n*4-3))+float(tmp(:,j,k,n*4-2))+ &
                        float(tmp(:,j,k,n*4-1))+float(tmp(:,j,k,n*4  ))) &
                     ) / (4.*nx)
    enddo
    enddo
    enddo
    um(:,:,:) = um(:,:,:)*real(sf)+real(ao)

    deallocate( tmp )

    istat = nf_get_var_real(ncid2,epdid,epdrag)

    istat = nf_close(ncid1)
    istat = nf_close(ncid2)

    epdrag = epdrag / 86400.

    return
    end subroutine readdata

!-----------------------------------------------------------------------

    subroutine dump

    character(len=100) :: wfn1
  
    integer :: istat ,ncid
    integer :: latdid,zdid  ,lrdid
    integer :: latid ,zid   ,lrid
    integer :: wid1  ,wid2  ,wid3  ,lid1

    write(wfn1,'(A,i4.4,a,i2.2,a)') 'res/wres_dyn.nc'

    istat = nf_create(trim(wfn1),nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',nyo ,latdid)
    istat = nf_def_dim(ncid,'p'  ,nz  ,zdid  )
    istat = nf_def_dim(ncid,'lat_rng',nyo+1,lrdid)
    istat = nf_def_var(ncid,'lat'    ,nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'p'      ,nf_real,1,zdid  ,zid  )
    istat = nf_def_var(ncid,'lat_rng',nf_real,1,lrdid ,lrid)
    istat = nf_def_var(ncid,'w_tot'  ,nf_real,2,(/latdid,zdid/),wid1)
    istat = nf_def_var(ncid,'w_ut'   ,nf_real,2,(/latdid,zdid/),wid2)
    istat = nf_def_var(ncid,'w_ep'   ,nf_real,2,(/latdid,zdid/),wid3)
    istat = nf_enddef(ncid)

    istat = nf_put_att_real(ncid,lrid,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wid2,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wid3,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_var_real(ncid,latid,lato)
    istat = nf_put_var_real(ncid,zid  ,lev )
    istat = nf_put_var_real(ncid,lrid ,lat(ilat_avg))
    istat = nf_put_var_real(ncid,wid1 ,waa_all(:,:,1))
    istat = nf_put_var_real(ncid,wid2 ,waa_all(:,:,2))
    istat = nf_put_var_real(ncid,wid3 ,waa_all(:,:,3))
    istat = nf_close(ncid)

    return
    end subroutine dump

    end program w_randel
