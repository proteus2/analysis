 par1_season = "DJF"
 par1_y1 = 1979 ;1958 ;1979
 par1_y2 = 2012 ;2015 ;2012
 par1_expgroup_v = "ERA-I"

 par1_p_v = 70

 par1_nma = 3
 par1_nph = 4

 if (par1_expgroup_v .eq. "ERA-I") then
   par1_expname_v   = "era-int"
   par1_fname_t     = "t"
   par1_fname_z     = "gp"
   par1_fname_uv    = (/"u","v"/)
   par1_fname_w     = "ome"
   par1_fname_pr    = "precip"
   par1_varname_t   = "T_GDS0_ISBL"
   par1_varname_z   = "Z_GDS0_ISBL"
   par1_varname_uv  = (/"U","V"/)+"_GDS0_ISBL"
   par1_varname_w   = "W_GDS0_ISBL"
   par1_varname_pr  = "TP_GDS0_SFC"
   par1_varname_t2m = "t2"
   par1_y0 = 1979
   par1_y9 = 2012
 end if

 par1_fname_u0 = "$DATD/clim/"+par1_expgroup_v+"/"+  \
                 par1_expname_v+".u_ypt."+par1_y0+"-"+par1_y9+".nc"

 par1_case_v = "ra"

 if ( any( par1_season .eq. (/"DJF","MAM","JJA","SON"/) ) ) then
   par1_nma = 3
   par1_iseason = ind( par1_season .eq. (/"DJF","MAM","JJA","SON"/) )
   par1_imon0 = (/ (/0,1,11/), (/2,3,4/), (/5,6,7/),(/8,9,10/) /)
 end if

