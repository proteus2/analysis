
undef("qbo_phase_ind")

function qbo_phase_ind(season:string, nph:integer, filename:string,  \
                       varname:string, p1:numeric, p2:numeric, lat0:numeric,  \
                       minmax_p1:numeric, minmax_p2:numeric,  \
                       y0:integer, y1:integer, y2:integer)

local vout, ncri, umin_p1, umax_p1, umin_p2, umax_p2, nma, iseason, imon0,  \
      nt, it0, it9, u, f1, uc, ua, x, y, deg, deg00, deg0, phs, tmpi,  \
      time_tag, nt_phase, i_phases, n

begin

;--------------------------------------- parameters
 ncri = 2

 umin_p1 = min(minmax_p1)
 umax_p1 = max(minmax_p1)
 umin_p2 = min(minmax_p2)
 umax_p2 = max(minmax_p2)

 if ( any( season .eq. (/"DJF","MAM","JJA","SON"/) ) ) then
   nma = 3
   iseason = ind( season .eq. (/"DJF","MAM","JJA","SON"/) )
   imon0 = (/ (/0,1,11/), (/2,3,4/), (/5,6,7/),(/8,9,10/) /)
 end if

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 u = new((/ncri,nt/),float,1.e20)
 f1 = addfile(filename,"r")
 u(0,:) = tofloat( dim_avg( f1->$varname$(it0:it9,{p1},{-lat0:lat0}) ) )
 u(1,:) = tofloat( dim_avg( f1->$varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) )
 if (dimsizes(p2) .eq. 2) then
   u(1,:) = dim_avg( f1->$varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) -  \
            dim_avg( f1->$varname$(it0:it9,{p2(1)},{-lat0:lat0}) )
 end if
 delete(f1)

;--------------------------------------- process 0
 uc = new((/ncri,12/),float,1.e20)
 do n=0, 11
   uc(:,n) = dim_avg(u(:,n::12))
 end do

 ua = u
 do n=0, 11
   ua(:,n::12) = u(:,n::12) - conform_dims((/ncri,nt/12/),uc(:,n),0)
 end do


 x = ua(0,:) - 0.5*tofloat(umin_p1+umax_p1)
 y = (ua(1,:) - 0.5*(umin_p2+umax_p2))*  \
     (tofloat(umax_p1-umin_p1)/tofloat(umax_p2-umin_p2))

 deg = where(x .ne. 0., atan(y/x)*180./3.141592, 90.)
 deg = where(abs(y) .gt. abs(x) .and. x*y .gt. 0.,  \
             90. - atan(x/y)*180./3.141592, deg)
 deg = where(abs(y) .gt. abs(x) .and. x*y .lt. 0.,  \
             -(90. - (atan(abs(x/y))*180./3.141592)), deg)
 deg = where(x .eq. 0. .and. y .lt. 0., -90., deg)
 deg = where(x .lt. 0., deg+180., deg)      ; -90 <= d < 270
 deg = 90. - deg                            ; -180 < d <= 180 (clockwise)
 deg = where(deg .lt. 0., deg+360., deg)    ;   0 <= d < 360

 deg00 = 45.
 deg0 = fspan(deg00,deg00+360.-(360./tofloat(nph)),nph)

 phs = where(deg .lt. deg0(0) .or. deg .ge. deg0(nph-1), 1, 999)
 do i=0, nph-2
   phs = where(deg .ge. deg0(i) .and. deg .lt. deg0(i+1), i+2, phs)
 end do

 tmpi = ispan(0,nt-1,1)
 time_tag = tmpi
 time_tag = -999
 do n=0, nma-1
   time_tag = where(mod(tmpi,12) .eq. imon0(iseason,n), iseason, time_tag)
 end do
 delete(tmpi)

 nt_phase = new(nph,integer)
 i_phases = new((/nph,nt/),integer,-999)
 do n=0, nph-1
   tmpi = ind( phs .eq. n+1 .and. time_tag .eq. iseason )
   nt_phase(n) = dimsizes(tmpi)
   i_phases(n,:nt_phase(n)-1) = tmpi
   delete(tmpi)
 end do

 vout = i_phases(:,:max(nt_phase)-1)

 vout@nt_phase = nt_phase

 vout@i_allphases = ind( time_tag .eq. iseason )

 vout@phases = phs(vout@i_allphases)
 vout@degree = deg(vout@i_allphases)

 vout@i_phase1 = i_phases(0,:nt_phase(0)-1)
 vout@i_phase2 = i_phases(1,:nt_phase(1)-1)
 if (nph .gt. 2) then
   vout@i_phase3 = i_phases(2,:nt_phase(2)-1)
   if (nph .gt. 3) then
     vout@i_phase4 = i_phases(3,:nt_phase(3)-1)
     if (nph .gt. 4) then
       vout@i_phase5 = i_phases(4,:nt_phase(4)-1)
       if (nph .gt. 5) then
         vout@i_phase6 = i_phases(5,:nt_phase(5)-1)
         if (nph .gt. 6) then
           vout@i_phase7 = i_phases(6,:nt_phase(6)-1)
           vout@i_phase8 = i_phases(7,:nt_phase(7)-1)
         end if
       end if
     end if
   end if
 end if

 return(vout)

end

