;23456
      begin

;      grib_in = addfile("fnl_040601_00_00.grb","r")
 
      do tt=0,0

      grib_in = addfile("/export32/kjs/fnl_data/2004/06/fnl_040601_" + sprinti("%2.2i",tt) + "_00.grb","r")
      names   = getfilevarnames(grib_in)
      print(names)

      lev    = grib_in->lv_ISBL3
      lat    = grib_in->lat_3
      lon    = grib_in->lon_3
      zonvel = grib_in->U_GRD_3_ISBL_10
      mervel = grib_in->V_GRD_3_ISBL_10
      temp   = grib_in->TMP_3_ISBL_10
      geop   = grib_in->HGT_3_ISBL_10

      print (lev)
      print (lat)
      print (lon)

      nlev=dimsizes(lev)
      nlon=dimsizes(lon)
      nlat=dimsizes(lat)
 
;      system("rm /export15/kjs/raytrace/INPUT/GRIB/ddata01.nc")
      ncdf = addfile("fnl_040601_" + sprinti("%2.2i",tt) + "_00.nc","c")

      fAtt            =  True    
      fAtt@title      =  "GRIB-to-Netcdf of FNL-data" 
      fAtt@sourcefile =  "fnl_040601_00_00.grb"
      fileattdef (ncdf,fAtt)
 
      dimNames =  (/"lev", "lat", "lon"/)
      dimSizes =  (/nlev, nlat, nlon/)
      dimUnlim =  (/False, False, False/)
      filedimdef (ncdf, dimNames, dimSizes, dimUnlim)

      filevardef    (ncdf, "lev", typeof(lev) , "lev")
      filevarattdef (ncdf, "lev", lev)
      
      filevardef    (ncdf, "lat", typeof(lat), "lat")
      filevarattdef (ncdf, "lat", lat)

      filevardef    (ncdf, "lon", typeof(lon), "lon")
      filevarattdef (ncdf, "lon", lon)
      
      filevardef (ncdf, "V", typeof(mervel), (/"lev", "lat", "lon"/) )
      filevarattdef(ncdf, "V", mervel)

      filevardef (ncdf, "U", typeof(zonvel), (/"lev", "lat", "lon"/) )
      filevarattdef(ncdf, "U", zonvel)

      filevardef (ncdf, "T", typeof(temp), (/"lev", "lat", "lon"/) )
      filevarattdef(ncdf, "T", temp)

      filevardef (ncdf, "GEOP", typeof(geop), (/"lev", "lat", "lon"/) )
      filevarattdef(ncdf, "GEOP", geop)

      ncdf->lev   =  (/lev/) 
      ncdf->lat   =  (/lat/) 
      ncdf->lon   =  (/lon/) 

      ncdf->V     =  (/mervel/)
      ncdf->U     =  (/zonvel/)
      ncdf->T     =  (/temp/)
      ncdf->GEOP  =  (/geop/)

      delete(grib_in)
      delete(lev)
      delete(lon)
      delete(lat)
      delete(zonvel)
      delete(mervel)
      delete(temp)
      delete(geop)
      delete(ncdf)

      end do

      end
