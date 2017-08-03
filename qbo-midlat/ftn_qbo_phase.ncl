
undef("qbo_phase_ind")

function qbo_phase_ind(season:string, filename:string, varname:string,  \
                       p1:numeric, p2:numeric, lat0:numeric,  \
                       y0:integer, y1:integer, y2:integer)

local vout, season, filename, varname, p1, p2, lat0, y0, y1, y2,  \
      ncri, nma, iseason, imon0, nt, it0, it9, its, u, f1, uc, ua,  \
      nph, phs, deg,  \
      nt_phase, i_phases, n, tmpi

begin

;--------------------------------------- parameters
 ncri = 2

 if ( any( season .eq. (/"DJF","MAM","JJA","SON"/) ) ) then
   nma = 3
   iseason = ind( season .eq. (/"DJF","MAM","JJA","SON"/) )
   imon0 = (/ (/0,1,11/), (/2,3,4/), (/5,6,7/),(/8,9,10/) /)
 end if

;--------------------------------------- read data
 nt = (y2-y1+1)*nma
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 its = new(nt,integer)
 do n=0, nma-1
   its(n::nma) = ispan(it0+imon0(iseason,n),it9,12)
 end do

 u = new((/ncri,nt/),float,1.e20)
 f1 = addfile(filename,"r")
 tmpi = dimsizes( f1->$varname$ )
 if (dimsizes(tmpi) .eq. 2) then
   if (tmpi(0) .lt. tmpi(1)) then
     u(0,:) = tofloat( f1->$varname$({p1},its) )
     u(1,:) = tofloat( f1->$varname$({p2(0)},its) )
     if (dimsizes(p2) .eq. 2) then
       u(1,:) = tofloat( f1->$varname$({p2(0)},its) -  \
                         f1->$varname$({p2(1)},its) )
     end if
   end if
 else
   u(0,:) = tofloat( dim_avg( f1->$varname$(its,{p1},{-lat0:lat0}) ) )
   u(1,:) = tofloat( dim_avg( f1->$varname$(its,{p2(0)},{-lat0:lat0}) ) )
   if (dimsizes(p2) .eq. 2) then
     u(1,:) = tofloat( dim_avg( f1->$varname$(its,{p2(0)},{-lat0:lat0}) ) -  \
                       dim_avg( f1->$varname$(its,{p2(1)},{-lat0:lat0}) ) )
   end if
 end if
 delete([/f1,tmpi/])

;--------------------------------------- process 0
 uc = new((/ncri,nma/),float,1.e20)
 do n=0, nma-1
   uc(:,n) = dim_avg(u(:,n::nma))
 end do

 ua = u
 do n=0, nma-1
   ua(:,n::nma) = u(:,n::nma) - conform_dims((/ncri,nt/nma/),uc(:,n),0)
 end do

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ( season .eq. "DJF" .and. ( p1 .eq. 50 .and. p2 .eq. 20 ) ) then

 nph = 3

 phs = new(nt,integer)
 phs = 0
 phs = where( ua(0,:) .lt. -10. .and. ua(1,:) .gt.   5. , 1, phs)
 phs = where( ua(0,:) .gt.   0. .and. ua(1,:) .gt.   0. , 2, phs)
 phs = where( ua(0,:) .gt.   0. .and. ua(1,:) .lt. -10. , 3, phs)

 deg = phs  ; dummy

else

end if

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 nt_phase = new(nph,integer)
 i_phases = new((/nph,nt/),integer,-999)
 do n=0, nph-1
   tmpi = ind( phs .eq. n+1 )
   nt_phase(n) = dimsizes(tmpi)
   i_phases(n,:nt_phase(n)-1) = tmpi
   delete(tmpi)
 end do

 vout = i_phases(:,:max(nt_phase)-1)

 vout@n_phase = nph

 vout@nt_phase = nt_phase

 vout@phases = phs
 vout@degree = deg

 vout@ua1 = ua(0,:)
 vout@ua2 = ua(1,:)

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

