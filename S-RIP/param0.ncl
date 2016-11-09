
param0_ra_p = (/"ERAI","MERRA","MERRA2","CFSR","JRA55"/)
param0_ra_m = (/"ERAIm","MERRA2m","JRA55_m"/)
;param0_ra = param0_ra_p
;param0_ra = param0_ra_m
param0_ra = array_append_record(param0_ra_p,param0_ra_m,0)
param0_var = "T"  ;  T / u / v / w (omega)

param0_yyyymm1 = 198101
param0_yyyymm2 = 201012

param0_lat_b = 15.  ; 5.
param0_p     = 100
param0_p_rng = (/125,1.5/)

param0_i_wav = 0

param0_wavs_k_a1   = (/1 ,1 /)
param0_wavs_k_a2   = (/20,60/)
param0_wavs_ome_a1 = (/0.03333,0.75/)
param0_wavs_ome_a2 = (/0.4    ,1.5 /)

param0_nbv     = 2.4e-2
param0_h_scale = 6.6e3

param0_ndays_spec = 90
param0_wgt_spec = 3./2.  ; (3 mon) / (2 mon)

