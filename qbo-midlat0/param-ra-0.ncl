; param0_cexp    = "qbo_obs"
; param0_expname = "qbo_sonde"
; param0_u0varname = "u"
; param0_p1 = 50  ;50 ;70
; param0_p2 = 20  ;20 ;30
 
 param0_cexp    = "ERA-I"
 param0_expname = "era-int"
 param0_u0varname = "U_GDS0_ISBL"
 param0_p1 = 70  ;50 ;70
 param0_p2 = 30  ;20 ;30
 
; param0_cexp    = "ERA-I-nr"
; param0_expname = "era-int_f"
; param0_u0varname = "u"
; param0_p1 = 22 ;23  ; 22 (54)  ;  23 (66)
; param0_p2 = 17 ;20  ; 19 (29)  ;  20 (36)
 param0_pp = (/19.,23.,29.,36.,44.,54.,66.,80.,96./)  ; 17 ~ 25
 param0_ppi0 = 17

 param0_cexp_v = (/"ERA-I"/)
 param0_expname_v = (/"era-int"/)
 param0_fname_t = "t"
 param0_fname_z = "gp"
 param0_fname_uv = (/"u","v"/)
 param0_fname_w = "ome"
 param0_fname_pr = "precip"
 param0_varname_t = "T_GDS0_ISBL"
 param0_varname_z = "Z_GDS0_ISBL"
 param0_varname_uv = (/"U","V"/)+"_GDS0_ISBL"
 param0_varname_w = "W_GDS0_ISBL"
 param0_varname_pr = "TP_GDS0_SFC"
 param0_varname_t2m = "t2"
 param0_p_v = 70

 param0_y0 = 1979
 param0_y9 = 2012
 param0_fdir_u0 = "$DATD/clim/"+param0_cexp
 param0_fname_u0 = param0_fdir_u0+"/"+param0_expname+".u_ypt."+param0_y0+"-"+param0_y9+".nc"
 if (param0_expname .eq. "qbo_sonde") then
   param0_y0 = 1953
   param0_y9 = 2016
   param0_fdir_u0 = "$HOME/portal/"+param0_cexp
   param0_fname_u0 = param0_fdir_u0+"/"+param0_expname+".nc"
 end if

 param0_lat0 = 10.
 param0_y1 = 1979 ;1979 ;1958 ;1979
 param0_y2 = 2012 ;2012 ;2015 ;2012

 param0_nma = 3
 param0_nph = 4

 if (param0_p1 .eq. 70) then
   param0_umin_p1 = -16
   param0_umax_p1 = 10
 end if
 if (param0_p2 .eq. 30) then
   param0_umin_p2 = -30
   param0_umax_p2 = 30
 end if
 if (param0_p1 .eq. 50) then
   param0_umin_p1 = -25
   param0_umax_p1 = 15
 end if
 if (param0_p2 .eq. 20) then
   param0_umin_p2 = -28
   param0_umax_p2 = 32
 end if
 if (param0_expname .eq. "era-int_f") then
   ;                 (/19.,23.,29.,36.,44.,54.,66.,80.,96./)
   tmp_param0_umin = (/-28,-99,-30,-30,-99,-25,-18,-99,-99/)
   tmp_param0_umax = (/32 ,99 ,30 ,26 ,99 ,15 ,12 ,99 ,99 /)

   param0_umin_p1 = tmp_param0_umin(param0_p1-param0_ppi0)
   param0_umax_p1 = tmp_param0_umax(param0_p1-param0_ppi0)
   param0_umin_p2 = tmp_param0_umin(param0_p2-param0_ppi0)
   param0_umax_p2 = tmp_param0_umax(param0_p2-param0_ppi0)
 end if
 
