
    integer, parameter :: start_day=9, start_hour=13
    real, parameter :: lat = 34.
    character(len=2), parameter ::  speriod='p3'

    integer, parameter :: sx=43
    integer, parameter :: ex=129
    integer, parameter :: sy=72
    integer, parameter :: ey=158 

!== READ WRF OUTPUT =====================================================
 
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0912-1018/new_data/new/temperature/'//ncfilename, ncid)
