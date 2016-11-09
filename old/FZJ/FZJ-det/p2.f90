
    integer, parameter :: start_day=8, start_hour=7
    real, parameter :: lat = 27.
    character(len=2), parameter ::  speriod='p2'

    integer, parameter :: sx=32
    integer, parameter :: ex=116
    integer, parameter :: sy=70
    integer, parameter :: ey=154 

!== READ WRF OUTPUT =====================================================
 
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0806-0912/new_data/new/temperature/'//ncfilename, ncid)
