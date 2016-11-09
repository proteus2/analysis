load "$UTIL/ncl/header.ncl"

SONDE = "dat/"
iTim = ispan(36+2,755-2,1)
; differences:
; 1. de-seasonalization for each station/period
; 2. de-seasonalization by removing the 1st harmonic (1-yr sine function)
; 3. eofunc_ts_Wrap(,,opt_stand) with opt_stand = True and opt_stand@jopt = 1
; 4. dim_standardize_Wrap(,1)
; 5. shifting by (0.196, -0.067)

begin

files   = systemfunc ("ls "+SONDE+"qbo_sonde0.nc")
fin     = addfile (files, "r")
uwnd    = fin->u(:,:755)
uwnd    = runave_Wrap (uwnd, 5, 1)
do nl=0, dimsizes (uwnd&p)-1
  uwnd(nl,:)    = rmAnnCycle1D (uwnd(nl,:))
end do
;------------------------------------------------------------------------------
;yh  uwnd    := uwnd(:,iTim)
;yh+
tmp = uwnd(:,iTim)
delete(uwnd)
uwnd = tmp
delete(tmp)
;yh-
yrfrac  = uwnd&t
uwnd&p@units    = "hPa"
printVarSummary (uwnd)
;------------------------------------------------------------------------------
option          = True
option@jopt     = 1
eofs    = eofunc_Wrap (uwnd, 2, option)
eofs    = -1.*eofs
pcs     = eofunc_ts_Wrap (uwnd, eofs, False)
pcs     = dim_standardize_Wrap (pcs, 0)
pcs(0,:)= pcs(0,:) - 0.188
pcs(1,:)= pcs(1,:) + 0.107
printVarSummary (eofs)
printVarSummary (pcs)
;------------------------------------------------------------------------------
ampl    = pcs(0,:)
ampl    = sqrt (pcs(0,:)^2 + pcs(1,:)^2)
;------------------------------------------------------------------------------
pi      = 4.*atan(1.)
phse    = 0.5/pi*atan2(pcs(1,:),pcs(0,:))

end

