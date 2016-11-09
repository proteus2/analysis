
    integer, parameter :: start_day=7, start_hour=1
    real, parameter :: lat = 20.
    character(len=2), parameter ::  speriod='p1'

    integer, parameter :: sx=32
    integer, parameter :: ex=116
    integer, parameter :: sy=32
    integer, parameter :: ey=116 

!== READ WRF OUTPUT =====================================================
 
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
