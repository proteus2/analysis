
param0_fname = "dat/qbo_sonde0.nc"
;param0_fname = "dat/qbo_sonde_hr0.nc"
param0_vname = "u"

param0_t1 = 1979    ; 1979.01 -
param0_t2 = 2014.0  ; - 2013.12

param0_l_deseasonalize = True
param0_i_des_harmonics = 1   ; 0: total  /  1: 1th harmonic  /  2: harmonics 1+2
param0_i_des_stations  = 0   ; 0: entire period  /  1: each station

param0_i_runavg = 0

;fix+
param0_t1_1 = min((/param0_t1,1956/))
param0_t2_1 = 1967.667   ; - 1967.08
param0_t1_2 = 1967.667   ; 1967.09 -
param0_t2_2 = 1976.0     ; - 1975.12
param0_t1_3 = 1976.0     ; 1976.01 -
;fix-
param0_t2_3 = param0_t2

