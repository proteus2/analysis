MODULE m_etc

  implicit none

! common variables

  character(len=128) ::  f_namelist
  character(len=6)   ::  stid
  character(len=4)   ::  year
  character(len=2)   ::  month, date, hour


CONTAINS


SUBROUTINE input_arg

  integer  ::  iargc
  external ::  getarg, iargc

  if (iargc() /= 6)  call stop_message(  &
        'PROGRAM STOP - wrong arguments\n'//  &
        ' getprof <namelist file> <stid> <year> <month> <date> <hour>')

  call getarg(1,f_namelist)
  call getarg(2,stid )       ! station ID
  call getarg(3,year )
  call getarg(4,month)
  call getarg(5,date )
  call getarg(6,hour )

END subroutine input_arg

SUBROUTINE stop_message(mess)

  character(len=*), intent(in) ::  mess

  write(6,'(a)') mess

  STOP

END subroutine stop_message

END module m_etc

