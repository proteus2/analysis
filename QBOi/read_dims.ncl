
load "dat/coord_6hr"

;-- head of getlon_q ---------------------------------------------------
undef("getlon_q")
function getlon_q(model:string, filev:file, varname:string)
local lon, model, filev, varname, lon1, dlon, tag, ndim, idim, dimname
;-- body of getlon_q ---------------------------------------------------
begin
lon__MRI_1 = 0.
dlon__MRI = 1.125

  if ( isdefined("lon__"+str_sub_str(model,"-","_")+"_1") ) then
    if (model .eq. "MRI") then
      lon1 = lon__MRI_1
      dlon = dlon__MRI
    end if
    lon = (lon1 + ispan(0,round(360./dlon,3)-1,1))*dlon
    lon = where(lon .ge. 360., lon-360., lon)
    lon@units = "degrees_east"
    lon@long_name = "lonitude"
  else
    tag = 0
    ndim = dimsizes( dimsizes( filev->$varname$ ) )
    do idim=0, ndim-1
      dimname = filev->$varname$!idim
      if ( str_lower(dimname) .eq. "lon" .or.  \
           str_lower(dimname) .eq. "lonitude" .or.  \
           str_lower(dimname) .eq. "g0_lon_5" .or.  \
           str_lower(dimname) .eq. "g0_lon_3" ) then
        tag = 1
        lon = filev->$dimname$
        break
      end if
    end do
    if (tag .eq. 0) then
      do idim=0, ndim-1
        dimname = filev->$varname$!idim
        if ( isatt(filev->$dimname$, "axis") ) then
        if ( filev->$dimname$@axis .eq. "X" ) then
          tag = 1
          lon = filev->$dimname$
          break
        end if
        end if
      end do
    end if
    if (tag .eq. 0) then
      print("LONGITUDE NOT DEFINED")
      exit
    end if
  end if

  return(lon)

end

;-- head of getlat_q ---------------------------------------------------
undef("getlat_q")
function getlat_q(model:string, filev:file, varname:string)
local lat, model, filev, varname, tag, ndim, idim, dimname, ny
;-- body of getlat_q ---------------------------------------------------
begin
  if ( isdefined("lat__"+str_sub_str(model,"-","_")) ) then
    if (model .eq. "MRI") then
      lat = lat__MRI
    end if
    lat@units = "degrees_north"
    lat@long_name = "latitude"
  else
    tag = 0
    ndim = dimsizes( dimsizes( filev->$varname$ ) )
    do idim=0, ndim-1
      dimname = filev->$varname$!idim
      if ( str_lower(dimname) .eq. "lat" .or.  \
           str_lower(dimname) .eq. "latitude" .or.  \
           str_lower(dimname) .eq. "g0_lat_4" .or.  \
           str_lower(dimname) .eq. "g0_lat_2" ) then
        tag = 1
        lat = filev->$dimname$
        break
      end if
    end do
    if (tag .eq. 0) then
      do idim=0, ndim-1
        dimname = filev->$varname$!idim
        if ( isatt(filev->$dimname$, "axis") ) then
        if ( filev->$dimname$@axis .eq. "Y" ) then
          tag = 1
          lat = filev->$dimname$
          break
        end if
        end if
      end do
    end if
    if (tag .eq. 0) then
      print("LATITUDE NOT DEFINED")
      exit
    end if
    if ( isfilevar(filev, dimname+"_bnds") ) then
      lat@bnd1 = filev->$(dimname+"_bnds")$(:,0)
      lat@bnd2 = filev->$(dimname+"_bnds")$(:,1)
    end if
  end if

  if ( .not. isatt(lat, "bnd1") ) then
    lat@bnd1 = new(dimsizes(lat),typeof(lat))
    lat@bnd2 = new(dimsizes(lat),typeof(lat))
    ny = dimsizes(lat)
    lat@bnd1(1:ny-1) = 0.5*(lat(0:ny-2) + lat(1:ny-1))
    lat@bnd2(0:ny-2) = lat@bnd1(1:ny-1)
    lat@bnd1(0) = lat(0)*2. - lat@bnd1(1)
    lat@bnd2(ny-1) = lat(ny-1)*2. - lat@bnd2(ny-2)
  end if

  return(lat)

end

;-- head of getlev_q ---------------------------------------------------
undef("getlev_q")
function getlev_q(model:string, filev:file, varname:string)
local lev, model, filev, varname, tag, ndim, idim, dimname
;-- body of getlev_q ---------------------------------------------------
begin
  if ( isdefined("lev__"+str_sub_str(model,"-","_")) ) then
    if (model .eq. "MRI") then
      lev = lev__MRI
      lev@units = "hPa"
      lev@long_name = "pressure"
    end if
    if (model .eq. "ERA-I-nr") then
      lev = lev__ERA_I_nr
      lev@units = "hPa"
      lev@long_name = "pressure"
    end if
  else
    tag = 0
    ndim = dimsizes( dimsizes( filev->$varname$ ) )
    do idim=0, ndim-1
      dimname = filev->$varname$!idim
      if ( str_lower(dimname) .eq. "plev" .or.  \
           str_lower(dimname) .eq. "lev" .or.  \
           str_lower(dimname) .eq. "level" .or.  \
           str_lower(dimname) .eq. "p" .or.  \
           str_lower(dimname) .eq. "z" .or.  \
           str_lower(dimname) .eq. "eta" .or.  \
           str_lower(dimname) .eq. "ht" .or.  \
           str_lower(dimname) .eq. "hybrid_ht" .or.  \
           str_lower(dimname) .eq. "lv_isbl1" ) then
        tag = 1
        lev = filev->$dimname$
        break
      end if
    end do
    if (tag .eq. 0) then
      do idim=0, ndim-1
        dimname = filev->$varname$!idim
        if ( isatt(filev->$dimname$, "axis") ) then
        if ( filev->$dimname$@axis .eq. "Z" .or.  \
             filev->$dimname$@axis .eq. "P" ) then
          tag = 1
          lev = filev->$dimname$
          break
        end if
        end if
      end do
    end if
    if (tag .eq. 0) then
      print("VERTICAL LEVEL NOT DEFINED")
      exit
    end if
  end if

  return(lev)

end

