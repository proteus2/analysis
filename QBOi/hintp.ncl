
;-- head of xrescale_q -------------------------------------------------
undef("xrescale_q")
function xrescale_q(v:numeric, lon0:numeric, wnmin:integer, wnmax:integer,  \
                    dlono:numeric)
local vo, v, lon0, wnmin, wnmax, dlono,  \
      twopi, narr, ndim, nxi, narra, nxo, lono, narr, xrad, fc, fc1,  \
      fc1c, fc1s, dimox, fcc, fcs
;-- body of xrescale_q -------------------------------------------------
begin
twopi = 6.2831853

  narr = dimsizes(v)
  ndim = dimsizes(narr)
  nxi = narr(ndim-1)
  narra = product(narr)/2

  nxo = toint(360./dlono)
  lono = ispan(0,nxo,1)*dlono

  narr(ndim-1) = nxo+1

  vo = new(narr,float)
  vo = 0.
  xrad = tofloat(ispan(0,nxo,1))/tofloat(nxo)*twopi
  xrad(nxo) = 0.
  xrad = xrad - tofloat(lon0)
  fc = ezfftf(v)
  fc1 = ndtooned(fc)
  fc1c = fc1(:narra-1)
  fc1s = fc1(narra:)

  dimox = ispan(0,ndim-2,1)
  do ik=wnmin-1, wnmax-1
    fcc = conform_dims(narr,onedtond(fc1c(ik::nxi/2),narr(dimox)),dimox)
    fcs = conform_dims(narr,onedtond(fc1s(ik::nxi/2),narr(dimox)),dimox)
    vo = vo + ( fcc*conform_dims(narr,cos(xrad*(ik+1)),ndim-1) +  \
                fcs*conform_dims(narr,sin(xrad*(ik+1)),ndim-1) )
  end do
  vo!(ndim-1) = "lon"
  vo&lon = lono

  vo@xbar = fc@xbar

  return(vo)

end

;-- head of yintp_q ----------------------------------------------------
undef("yintp_q")
function yintp_q(v:numeric, lat:numeric, dlato:numeric, bnd:numeric)
local vo, v, lat, dlato, bnd,  \
      ndim, nyh, lato, lata
;-- body of yintp_q ----------------------------------------------------
begin

  ndim = dimsizes(dimsizes(v))

  nyh = toint(bnd/dlato)
  lato = ispan(-nyh,nyh,1)*dlato

  lata = lat
  if (lat(1) .lt. lat(0)) then
    lata = -lat
  end if

  vo = linint1_n_Wrap(lata,v,False,lato,0,ndim-2)

  if (lat(1) .lt. lat(0)) then
    lato = -lato
  end if

  vo!(ndim-2) = "lat"
  vo&lat = lato

  return(vo)

end

