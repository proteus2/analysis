undef("disp_eq")
function disp_eq(kwn:numeric, freq:numeric, l_radian:logical, n_bv:numeric)
local disp, kwn, freq, l_radian, n_bv, two_pi, r_e, beta, nk, no, k, o,  \
      oo, okbi, nbeta, okbi1okbi, c_n, c2_n
begin

 two_pi = 2.*(4.*atan(1.))
 r_e = 6371.e3
 beta = 2.*(7.292116e-5)/r_e

 nk = dimsizes(kwn)
 no = dimsizes(freq)

 if ( l_radian ) then
   ; [rad] -> [rad/m]  ;  [rad/day] -> [rad/s]
   k = conform_dims((/no,nk/), kwn/(two_pi*r_e), 1)
   o = conform_dims((/no,nk/), freq/86400.     , 0)
 else
   ; [cyc] -> [rad/m]  ;  [cyc/day] -> [rad/s]
   k = conform_dims((/no,nk/), kwn/r_e           , 1)
   o = conform_dims((/no,nk/), freq*two_pi/86400., 0)
 end if
 
 oo = o^2
 okbi = o*k/beta
 nbeta = n_bv*beta
 okbi1okbi = okbi*(1.+okbi)

 c_n = ispan(1,4,1) + 0.5   ; for n = 1 to 4
 c2_n = c_n^2

 disp = new((/2,3,no,nk/),float,1.e30)
 disp = 1.e3


 k@_FillValue = 0.  ; not to divide by 0

 ; Kelvin / symm
 disp(0,0,:,:) = two_pi/1.e3*o/(n_bv*k)

 ; MRG / anti-symm
 disp(1,0,:,:) = two_pi/1.e3*oo/(nbeta*(1.+okbi))

 ; IG n=1 / symm
 disp(0,1,:,:) = two_pi/1.e3*oo/  \
      (nbeta*(c_n(0)+sqrt(c2_n(0)+okbi1okbi)))

 ; IG n=2 / anti-symm
 disp(1,1,:,:) = two_pi/1.e3*oo/  \
      (nbeta*(c_n(1)+sqrt(c2_n(1)+okbi1okbi)))

 okbi1okbi@_FillValue = 0.  ; not to divide by 0

 ; Rossby n=1 / symm
 disp(0,2,:,:) = two_pi/1.e3*oo/  \
      (nbeta*(c_n(0)-sqrt(c2_n(0)+okbi1okbi)))

 ; Rossby n=2 / anti-symm
 disp(1,2,:,:) = two_pi/1.e3*oo/  \
      (nbeta*(c_n(1)-sqrt(c2_n(1)+okbi1okbi)))

 disp = where(disp .lt. 0., disp@_FillValue, disp)

 disp(:,2,0,:) = 0.
 disp(0,2,:,:) = where(okbi .gt. -1., disp(0,2,:,:), disp@_FillValue)
 disp(1,2,:,:) = where(okbi .gt. -1., disp(1,2,:,:), disp@_FillValue)

 disp@Kelvin  = (/0,0/)
 disp@MRG     = (/1,0/)
 disp@IG1     = (/0,1/)
 disp@IG2     = (/1,1/)
 disp@Rossby1 = (/0,2/)
 disp@Rossby2 = (/1,2/)

return(disp)
end
