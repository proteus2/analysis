 param0_cexp    = (/"CGW"/)
 param0_expname = (/"tc" /)
; param0_cexp    = (/"CGW","CTL"/)
; param0_expname = (/"tc" ,"tk" /)
 param0_fdir_u0 = "$DATD/AOL60"+param0_cexp
 param0_fname_u0 = param0_fdir_u0+"/xm.uan"+param0_expname+"a.pj.nc"
 param0_u0varname = "u_p_uv"
 param0_lat0 = 10.
 param0_p1 = 70
 param0_p2 = 30

 param0_lowlev = False
 param0_varname_t = "T_p_uv"
 param0_varname_z = "ht_p_uv"
 param0_varname_uv = (/"u","v"/)+"_p_uv"
 param0_p_v = 500
; param0_lowlev = True
 param0_varname_ll = "theta"
 param0_z_v_ll = 1.2e3

 param0_y0 = 1951
 param0_y1 = 1951
 param0_y2 = 2000 ; 2097

 param0_nma = 3
 param0_nph = 4

 if (param0_p1 .eq. 70) then
   param0_umin_p1 = -12
   param0_umax_p1 = 8
 end if
 if (param0_p2 .eq. 30) then
   param0_umin_p2 = -22
   param0_umax_p2 = 18
 end if

