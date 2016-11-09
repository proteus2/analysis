undef("eofqbo")
function eofqbo(u:numeric, opteof:integer, optets:integer, l_normalize:logical)
local tsn, u, opteof, optets,  \
      ud, optEOF, optETS, v_eof, im, ts, n_nonmiss, sd, da0
begin

 ; remove mean
 ud = u - conform_dims(dimsizes(u),dim_avg(u),0)

 optEOF = True
 optEOF@jopt = opteof  ; covariance / correlation

 optETS = True
 optETS@jopt = optets

 v_eof = eofunc(ud,2,optEOF)
 do im=0, 1
   if (sum(v_eof(im,:)) .lt. 0.) then
     v_eof(im,:) = v_eof(im,:)*(-1.0)
   end if
 end do

 print("variance by the two EOFs: "+sum(v_eof@pcvar)+"%")

 ts = eofunc_ts(ud,v_eof,optETS)

 tsn = ts
 if ( l_normalize ) then
   n_nonmiss = dimsizes(ind(.not. ismissing(ts(0,:))))
   sd = dim_stddev(ts)*tofloat(n_nonmiss-1)/tofloat(n_nonmiss)
   tsn = ts/conform_dims(dimsizes(ts),sd,0)
 end if

 da0 = sqrt( tsn(0,:)^2 + tsn(1,:)^2 )
 print("variance of distances (not shifted) : "+avg((da0 - avg(da0))^2))
 
 tsn@eof = v_eof

return(tsn)
end

;-------------------------------------------------------------------------------
 
undef("shift_eof")
function shift_eof(tsn:numeric, g100res:numeric, n_iter:integer)
local tsn_s, tsn, g100res, n_iter,  \
      v_shift, ntsn, npt, npth, ndima, dv, da,  \
      i_shift, j_shift, gres, iter_s, ji
begin

 tsn_s = tsn
 v_shift = (/0.,0./)

 ntsn = dimsizes(tsn(0,:))
 npt = 201
 npth = npt/2
 ndima = (/npt,npt,ntsn/)
 dv = new((/npt,npt/),float)
 da = new(ndima,float)

 i_shift = conform_dims(ndima,tofloat(ispan(-npth,npth,1)),1)
 j_shift = conform_dims(ndima,tofloat(ispan(-npth,npth,1)),0)
 
 gres = g100res  ; grid resolution

 do iter_s=0, n_iter-1
   print("iteration "+(iter_s+1))

   ; distances
   da = sqrt( (conform_dims(ndima,tsn_s(0,:),2) - i_shift*gres)^2 +  \
              (conform_dims(ndima,tsn_s(1,:),2) - j_shift*gres)^2 )
   ; dv: variance of distances
   da = da - conform_dims(ndima,dim_avg(da),(/0,1/))
   dv = dim_avg(da^2)

   ji = ind_resolve(minind(ndtooned(dv)),dimsizes(dv))

   v_shift = v_shift + tofloat( (/ji(0,1),ji(0,0)/) - npth )*gres
   tsn_s(0,:) = tsn_s(0,:) - tofloat(ji(0,1) - npth)*gres
   tsn_s(1,:) = tsn_s(1,:) - tofloat(ji(0,0) - npth)*gres
   print("shift pt. : ("+v_shift(0)+","+v_shift(1)+")")
   print("variance of distances : "+dv(ji(0,0),ji(0,1)))

   gres = gres*0.1
 end do

 tsn_s@xy_shift = v_shift
 tsn_s@var_of_d = dv(ji(0,0),ji(0,1))

return(tsn_s)
end

;-------------------------------------------------------------------------------
 
undef("phsamp_eof")
function phsamp_eof(tsn)
local phs_amp, tsn, ampm, pi, phase, aphase, i_nonmiss, i0_nm, i9_nm, n, p_rate
begin
 
 amp = sqrt( tsn(0,:)^2 + tsn(1,:)^2 )
 ampm = avg(amp)
 ; The amplitudes are normalized such that their average becomes sqrt(2).
 amp = amp*(sqrt(2.)/ampm)

 pi = 4.*atan(1.)

 phase = atan2(tsn(1,:),tsn(0,:))/(2.*pi)
 phase = where(phase .lt. 0., phase+1., phase)

 aphase = phase

 i_nonmiss = ind(.not. ismissing(aphase))
 i0_nm = min(i_nonmiss)
 i9_nm = max(i_nonmiss)
 do n=i0_nm+1, i9_nm
   if ( aphase(n) - aphase(n-1) .lt. -0.5 ) then
     aphase(n:) = aphase(n:) + 1.
   end if
   if ( aphase(n) - aphase(n-1) .gt. 0.5 ) then
     aphase(n:) = aphase(n:) - 1.
   end if
 end do

 p_rate = new(dimsizes(aphase),float,aphase@_FillValue)
 p_rate(i0_nm+1:i9_nm-1) = (aphase(i0_nm+2:i9_nm) - aphase(i0_nm:i9_nm-2))/2.

 phs_amp = aphase
 phs_amp@phase = phase
 phs_amp@aphase = aphase
 phs_amp@p_rate = p_rate
 phs_amp@ampl = amp

return(phs_amp)
end

