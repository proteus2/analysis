
;-- head of lat_intg ---------------------------------------------------
undef("lat_intg")
function lat_intg(lat:numeric, lat_bnd1:numeric, lat_bnd2:numeric,  \
                  lata:numeric, var:numeric, idim:integer, wgt:numeric)
local varo, lat, lat_bnd1, lat_bnd2, lata, var, idim, wgt,  \
      lata_l, lata_u, dlat, j1, j2, sumdlat, dlat1b, dlat2b, ndims,  \
      vari, dlat1o, dlat2o, varb1i, varb2i, varb1o, varb2o, ub1, ub2
;-- body of lat_intg ---------------------------------------------------
begin
  lata_l = min(lata)
  lata_u = max(lata)
  dlat = abs(lat_bnd2 - lat_bnd1)
  if (lat(0) .lt. lat(1)) then
    j1 = min( ind( lat .ge. lata_l ) )
    j2 = max( ind( lat .le. lata_u ) )
  else
    j1 = min( ind( lat .le. lata_u ) )
    j2 = max( ind( lat .ge. lata_l ) )
  end if
  dlat(j1) = dlat(j1)*0.5
  dlat(j2) = dlat(j2)*0.5
  sumdlat = tofloat(sum(dlat(j1:j2)))

  dlat1b = abs(lata_l) - tofloat(abs(lat(j1)))
  dlat2b = abs(lata_u) - tofloat(abs(lat(j2)))

  ndims = dimsizes(var)
  ndims(idim) = j2-j1+1

  if (dimsizes(ndims) .eq. 1) then
    vari = var(j1:j2)
  else if (dimsizes(ndims) .eq. 2) then
    if (idim .eq. 0) then
      vari = var(j1:j2,:)
    else if (idim .eq. 1) then
      vari = var(:,j1:j2)
    end if
    end if
  else if (dimsizes(ndims) .eq. 3) then
    if (idim .eq. 1) then
      vari = var(:,j1:j2,:)
    else if (idim .eq. 2) then
      vari = var(:,:,j1:j2)
    end if
    end if
  else if (dimsizes(ndims) .eq. 4) then
    if (idim .eq. 2) then
      vari = var(:,:,j1:j2,:)
    else if (idim .eq. 3) then
      vari = var(:,:,:,j1:j2)
    end if
    end if
  end if
  end if
  end if
  end if
  if ( .not. isdefined("vari") ) then
    print("Check idim in lat_intg.")
    exit
  end if
 
  if (dimsizes(ndims) .eq. 1) then
    if (dimsizes(wgt) .eq. 1) then
      varo = tofloat( sum( vari*dlat(j1:j2) ) )
    else
      varo = tofloat( sum( vari*(dlat(j1:j2)*wgt(j1:j2)) ) )
    end if
  else 
    if (dimsizes(wgt) .eq. 1) then
      varo = tofloat( dim_sum_n( vari*  \
             conform_dims(ndims, dlat(j1:j2), idim), idim ) )
    else
      varo = tofloat( dim_sum_n( vari*  \
             conform_dims(ndims, dlat(j1:j2)*wgt(j1:j2), idim), idim ) )
    end if
  end if
 
  if ( dlat1b .ne. 0. .or. dlat2b .ne. 0. ) then
    dlat1o = tofloat(lat(j1) - lat(j1-1))
    dlat2o = tofloat(lat(j2+1) - lat(j2))

    sumdlat = sumdlat + (dlat1b + dlat2b)
;    sumdlat = lata*2.

    if (dimsizes(ndims) .eq. 1) then
      varb1i = var(:,j1,:)
      varb2i = var(:,j2,:)
      varb1o = var(:,j1-1,:)
      varb2o = var(:,j2+1,:)
    else if (dimsizes(ndims) .eq. 2) then
      if (idim .eq. 0) then
        varb1i = var(j1,:)
        varb2i = var(j2,:)
        varb1o = var(j1-1,:)
        varb2o = var(j2+1,:)
      else if (idim .eq. 1) then
        varb1i = var(:,j1)
        varb2i = var(:,j2)
        varb1o = var(:,j1-1)
        varb2o = var(:,j2+1)
      end if
      end if
    else if (dimsizes(ndims) .eq. 3) then
      if (idim .eq. 1) then
        varb1i = var(:,j1,:)
        varb2i = var(:,j2,:)
        varb1o = var(:,j1-1,:)
        varb2o = var(:,j2+1,:)
      else if (idim .eq. 2) then
        varb1i = var(:,:,j1)
        varb2i = var(:,:,j2)
        varb1o = var(:,:,j1-1)
        varb2o = var(:,:,j2+1)
      end if
      end if
    else if (dimsizes(ndims) .eq. 4) then
      if (idim .eq. 2) then
        varb1i = var(:,:,j1,:)
        varb2i = var(:,:,j2,:)
        varb1o = var(:,:,j1-1,:)
        varb2o = var(:,:,j2+1,:)
      else if (idim .eq. 3) then
        varb1i = var(:,:,:,j1)
        varb2i = var(:,:,:,j2)
        varb1o = var(:,:,:,j1-1)
        varb2o = var(:,:,:,j2+1)
      end if
      end if
    end if
    end if
    end if
    end if
    if ( .not. isdefined("varb1o") ) then
      print("Check j1 and j2 in lat_intg.")
      exit
    end if
 
;coslat !!
    ub1 = ( (0.5*dlat1b)*varb1o + (dlat1o - 0.5*dlat1b)*varb1i )/dlat1o
    ub2 = ( (0.5*dlat2b)*varb2o + (dlat2o - 0.5*dlat2b)*varb2i )/dlat2o

    varo = varo + tofloat(ub1*dlat1b + ub2*dlat2b)
  end if

  varo@sumdlat = sumdlat

;  varo = varo/sumdlat

  return(varo)

end

