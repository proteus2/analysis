
 module params

 implicit none

 double precision      :: zb , zt , ntr, nst, ubt, vbt, cqx, cqy
 double precision      :: zs , ubb, vbb, ub0, vb0
 double precision      :: ut, ub, u0, cq

 integer            , parameter      :: nc     = 241
 integer            , parameter      :: nd     = 1!36
 double precision   , parameter      :: dc     = 0.5
 double precision   , parameter      :: dd     = 5.
 double precision, parameter         :: epsabs = 0.01
 double precision, parameter         :: epsrel = 0.01
 integer                             :: key = 6
 double precision   , dimension(nc)  :: cpgrd

 double precision, parameter ::  delh = 5.d3, delt = 1200.

 end module params
