   begin

     infname="T799-0.25/od_oper_an_ml_20060707_20060806"
     outfname="t_od_oper_an_ml_20060707_20060806"

; get variable names from grib file
     grib_in  = addfile(infname+".grb","r")   
     names    = getfilevarnames(grib_in); extract all variable names 

; create output netcdf file
     ncdf_out = addfile(outfname+".nc" ,"c")       ; create output netCDF file

; specify desired GRIB variables
     grib_names = (/ "g0_lat_4", "g0_lon_5", "T_GDS0_HYBL" /)
     ncl_names  = (/ "g0_lat_4", "g0_lon_5", "T_GDS0_HYBL" /)

; loop through variables and output each to netcdf
     do i = 0, dimsizes(ncl_names)-1  
       ncdf_out->$ncl_names(i)$ = grib_in->$grib_names(i)$
     end do   

   end 
