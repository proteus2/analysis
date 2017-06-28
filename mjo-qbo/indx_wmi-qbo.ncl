;==============================================================================
; load ncl libraries
;------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;==============================================================================
; set directory root
;------------------------------------------------------------------------------
CASE	= "f2000_qbo_uni_test"
MODEL   = "$DATD/mjo-qbo/"
OUTPUT	= "$DATD/mjo-qbo/"+CASE+"/"
FIGDIR	= "./fig/"
SCRIPT	= get_script_prefix_name ()
p0 = 70   ; 15  25  55  70  90  100
print ("run "+SCRIPT)
;------------------------------------------------------------------------------
Fout	= OUTPUT+"wmi_index."+CASE+".nc"
Fig	= FIGDIR+SCRIPT+"-p"+p0+"_"+CASE
;------------------------------------------------------------------------------
yr_sp = 6  ; spin-up years 
nyr = 24
if (CASE .eq. "f2000_cntl") then
  yr_sp = 0
  nyr = 10
end if
;------------------------------------------------------------------------------
; read ucomp
;------------------------------------------------------------------------------
files	= systemfunc ("ls "+MODEL+CASE+"/"+CASE+".var-mjo.remap.daily.nc")
fileq	= systemfunc ("ls "+MODEL+CASE+"/"+CASE+".qbo.monthly.nc")
fuo  	= addfile (files, "r")
fuq  	= addfile (fileq, "r")
;------------------------------------------------------------------------------
n1 = 365*yr_sp
n2 = (n1-1) + 365*nyr
u200	= fuo->U200(n1:n2,{-30:30},:)
printVarSummary (u200)
u850	= fuo->U850(n1:n2,{-30:30},:)
olr	= fuo->FLUT(n1:n2,{-30:30},:)
nm1 = 12*yr_sp
nm2 = (nm1-1) + 12*nyr
ntm = nm2 - nm1 + 1
uq0 = dim_avg_Wrap( fuq->U(nm1-1:nm2+1,:,{-10:10}) )
uq = transpose( uq0(1:ntm,:) )
uq = (/ transpose( (uq0(0:ntm-1,:) + uq0(1:ntm,:) + uq0(2:ntm+1,:))/3. ) /)
uq = uq - conform_dims(dimsizes(uq),dim_avg_Wrap(uq),0)
ndimq = dimsizes(uq)
nzq = ndimq(0)
pq = uq&$uq!0$
printVarSummary(uq)
uqsd05 = 0.5*dim_stddev(uq)
qboi = toint(uq*0.)
do k=0, nzq-1  ; DJF
  qboi(k,:) = where( uq(k,:) .gt. uqsd05(k) , 1 , qboi(k,:))
  qboi(k,:) = where( uq(k,:) .lt. -uqsd05(k), -1, qboi(k,:))
end do
;------------------------------------------------------------------------------
dims	= dimsizes (u200)
ntim	= dims(0)
nlat	= dims(1)
nlon	= dims(2)
;------------------------------------------------------------------------------
; compute daily climatologies and anomalies
;------------------------------------------------------------------------------
time	= u200&time
ymd	= ut_calendar (time, -2)
ymdhms	= ut_calendar (time, 0)
year	= floattointeger (ymdhms(:,0))
month	= floattointeger (ymdhms(:,1))
day	= floattointeger (ymdhms(:,2))
ddd	= day_of_year (year, month, day)
yyyyddd	= year*1000 + ddd
yrfrac	= yyyymmdd_to_yyyyfrac (ymd, 0.)
;------------------------------------------------------------------------------
timeq	= uq&time
ymdq	= ut_calendar (timeq, -2)
yrfracq	= yyyymmdd_to_yyyyfrac (ymdq, 0.)
;------------------------------------------------------------------------------
nhar	= 3
u200ClmDay	= clmDayTLL (u200, yyyyddd)
u200ClmDay_sm	= smthClmDayTLL (u200ClmDay, nhar)
u850ClmDay	= clmDayTLL (u850, yyyyddd)
u850ClmDay_sm	= smthClmDayTLL (u850ClmDay, nhar)
olrClmDay	= clmDayTLL (olr, yyyyddd)
olrClmDay_sm	= smthClmDayTLL (olrClmDay, nhar)
;------------------------------------------------------------------------------
u200Anom	= calcDayAnomTLL (u200, yyyyddd, u200ClmDay_sm)
u850Anom	= calcDayAnomTLL (u850, yyyyddd, u850ClmDay_sm)
olrAnom		= calcDayAnomTLL (olr, yyyyddd, olrClmDay_sm)
printVarSummary (olrAnom)
;------------------------------------------------------------------------------
; take average over latS to latN
;------------------------------------------------------------------------------
latS    = -15
latN    =  15
tmp	= u200Anom(time|:,lon|:,lat|:)
u2002D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
tmp	= u850Anom(time|:,lon|:,lat|:)
u8502D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
tmp	= olrAnom(time|:,lon|:,lat|:)
olr2D	= dim_avg_Wrap (tmp(:,:,{latS:latN}))
delete (tmp)
;------------------------------------------------------------------------------
; do band pass filtering
;------------------------------------------------------------------------------
ihp     = 2                             ; bpf=>band pass filter
nWgt    = 201
sigma   = 1.0                           ; Lanczos sigma
fca     = 1./100.
fcb     = 1./20.
wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
u200_bp = wgt_runave_Wrap (u2002D(lon|:,time|:), wgt, 0)
u850_bp = wgt_runave_Wrap (u8502D(lon|:,time|:), wgt, 0)
olr_bp	 = wgt_runave_Wrap (olr2D(lon|:,time|:), wgt, 0)
;olr3D_bp = wgt_runave_Wrap (olrAnom(lat|:,lon|:,time|:), wgt, 0)
;------------------------------------------------------------------------------
; remove means of band pass series
;------------------------------------------------------------------------------
u200_bp	= dim_rmvmean (u200_bp)
u850_bp	= dim_rmvmean (u850_bp)
olr_bp	= dim_rmvmean (olr_bp)
;olr3D_bp	= dim_rmvmean (olr3D_bp)
;------------------------------------------------------------------------------
; compute the temporal variance
;------------------------------------------------------------------------------
print ("compute variance and normalize")
var_u200	= dim_variance_Wrap (u200_bp)
var_u850	= dim_variance_Wrap (u850_bp)
var_olr		= dim_variance_Wrap (olr_bp)
;------------------------------------------------------------------------------
; compute the zonal mean of the temporal variance
;------------------------------------------------------------------------------
zavg_var_u200	= dim_avg_Wrap (var_u200)
zavg_var_u850	= dim_avg_Wrap (var_u850)
zavg_var_olr	= dim_avg_Wrap (var_olr)
;------------------------------------------------------------------------------
; normalize by sqrt (avg_var*)
;------------------------------------------------------------------------------
u200_bp	= u200_bp/sqrt(zavg_var_u200)
u850_bp	= u850_bp/sqrt(zavg_var_u850)
olr_bp	= olr_bp/sqrt(zavg_var_olr)
;------------------------------------------------------------------------------
; combine the normalized data into one variable
;------------------------------------------------------------------------------
cdata	= new ((/3*nlon,ntim/), double)
do i=0, nlon-1
  cdata(i,		:) = olr_bp(i,:)
  cdata(i+nlon,		:) = u200_bp(i,:)
  cdata(i+2*nlon,	:) = u850_bp(i,:)
end do
printVarSummary (cdata)
;------------------------------------------------------------------------------
; compute the combined EOF
;------------------------------------------------------------------------------
print ("do EOF analysis")
neof		= 2
eof_cdata	= eofunc (cdata, neof, False)
; --- sign of the EOFs ---
eof_olr1 = eof_cdata(0,0:nlon-1)
eof_olr1!0 = "lon"
eof_olr1&lon = olr_bp&lon
eof_u852 = eof_cdata(1,2*nlon:3*nlon-1)
eof_u852!0 = "lon"
eof_u852&lon = u850_bp&lon
if (sum(eof_olr1({60:180})) .gt. 0.) then
  eof_cdata(0,:) = eof_cdata(0,:)*(-1.)
end if
if (sum(eof_u852({60:180})) .lt. 0.) then
  eof_cdata(1,:) = eof_cdata(1,:)*(-1.)
end if
delete(eof_olr1)
delete(eof_u852)
; ------------------------
eof_ts_cdata	= eofunc_ts (cdata, eof_cdata, False)
printVarSummary (eof_cdata)
;------------------------------------------------------------------------------
; extract each variable
;------------------------------------------------------------------------------
nvar		= 3
ceof		= new ((/nvar,neof,nlon/), typeof (cdata))
do n=0, neof-1
  ceof(0,n,:)	= eof_cdata(n,0:nlon-1)
  ceof(1,n,:)	= eof_cdata(n,2*nlon:)
  ceof(2,n,:)	= eof_cdata(n,nlon:2*nlon-1)
end do
;------------------------------------------------------------------------------
; compute cross correlation of the EOFs; EOF 1 vs EOF 2
;------------------------------------------------------------------------------
mxlag		= 25
rlag_01		= esccr (eof_ts_cdata(0,:), eof_ts_cdata(1,:), mxlag)   ; (N,mxlag+1)
rlag_10		= esccr (eof_ts_cdata(1,:), eof_ts_cdata(0,:), mxlag)   ; (N,mxlag+1)
;------------------------------------------------------------------------------
ccr_12		= new ((/2*mxlag+1/), typeof(rlag_01))
ccr_12(mxlag:)	= rlag_10(0:mxlag)
ccr_12(0:mxlag)	= rlag_01(::-1)       ; reverse order
;------------------------------------------------------------------------------
; normalize the multivariate EOF 1&2 component time series
; compute (PC1^2+PC2^2): values > 1 indicate "strong" periods
;------------------------------------------------------------------------------
eof_ts_cdata(0,:)	= eof_ts_cdata(0,:)/stddev(eof_ts_cdata(0,:))
eof_ts_cdata(1,:)	= eof_ts_cdata(1,:)/stddev(eof_ts_cdata(1,:))
;------------------------------------------------------------------------------
rmm1		= eof_ts_cdata(0,:)
rmm2		= eof_ts_cdata(1,:)
indx            = sqrt (rmm1^2 + rmm2^2)
indx_smt        = runave(indx, 91, 0) ; 91-day running mean
indx_smt!0		= "year"
indx_smt&year	= yrfrac
;------------------------------------------------------------------------------
nGood   = num(.not.ismissing(indx))     ; # non-missing
nStrong = num(indx.ge.1.)
print("nGood="+nGood+"   nStrong="+nStrong+"   nOther="+(nGood-nStrong))
print ("mean amplitude = "+avg(indx))
;==============================================================================
; save mjo index
;------------------------------------------------------------------------------
print ("save the data")
indx!0		= "time"
indx&time	= olr&time
indx@long_name	= "MJO PC index"
indx@info	= "SQRT(PC1^2 + PC2^2)"
;------------------------------------------------------------------------------
rmm1@long_name		= "PC1"
rmm1@info		= "PC1/stddev(PC1)"
copy_VarCoords (indx, rmm1)
rmm2@long_name		= "PC2"
rmm2@info		= "PC2/stddev(PC2)"
copy_VarCoords (indx, rmm2)
;------------------------------------------------------------------------------
angBnd		= new ((/2,8/), "float")
angBnd(0,:)	= fspan (-180., 180.-45., 8)
angBnd(1,:)	= angBnd(0,:)+45.
r2d		= 180./(4.*atan(1.))
ang		= atan2(rmm2,rmm1)*r2d
phase		= new (ntim, "integer")
phase		= 0
do np=1, 8
  phase	= where (ang.ge.angBnd(0,np-1).and.ang.lt.angBnd(1,np-1), np, phase)
end do
;------------------------------------------------------------------------------
indx_mon = calculate_monthly_values(indx,"avg",0,False)
indx_mon(3    ) = indx_mon@_FillValue  ; can be different btw cases
indx_mon(ntm-4) = indx_mon@_FillValue
indx_mon = runave(indx_mon,3,0)
indx_wqbo = new((/nzq,(ntm/12)/),typeof(indx_mon))
indx_eqbo = new((/nzq,(ntm/12)/),typeof(indx_mon))
indx_nqbo = new((/nzq,(ntm/12)/),typeof(indx_mon))
do k=0, nzq-1  ; DJF
  indx_wqbo(k,:) = where(qboi(k,::12) .eq. 1 , indx_mon(::12), indx_mon@_FillValue)
  indx_eqbo(k,:) = where(qboi(k,::12) .eq. -1, indx_mon(::12), indx_mon@_FillValue)
  indx_nqbo(k,:) = where(qboi(k,::12) .eq. 0 , indx_mon(::12), indx_mon@_FillValue)
end do
print(round(pq,3)+": "+dim_avg(indx_wqbo)+" / "+dim_avg(indx_nqbo)+" / "+dim_avg(indx_eqbo))
;------------------------------------------------------------------------------
ndayavg = 5  ; coded below assuming as odd number
rmm1p		= runave(rmm1,ndayavg,0)
rmm2p		= runave(rmm2,ndayavg,0)
indxp       = sqrt (rmm1p^2 + rmm2p^2)
angp		= atan2(rmm2p,rmm1p)*r2d
phasecmp	= new (ntim, "integer")
phasecmp	= 0
do np=1, 8
  phasecmp = where (angp.ge.angBnd(0,np-1).and.angp.lt.angBnd(1,np-1), np, phasecmp)
end do
missv = -999
phasecmp@_FillValue = missv
phasecmp = where(indxp .le. 1., missv, phasecmp)
tmpi = nWgt/2 + ndayavg/2
tmp0 = phasecmp(tmpi)  ; keep the 1st non-missing value
tmp = phasecmp(1:) - phasecmp(:ntim-2)
phasecmp(1:) = where( tmp .eq. 1 .or. tmp .eq. 0 .or.  \
                      ( phasecmp(1:) .eq. 1 .and. phasecmp(:ntim-2) .eq. 8 ),  \
                      phasecmp(1:), missv )
phasecmp(tmpi) = tmp0
tmp1 = 1 + ind( (.not. ismissing(phasecmp(1:))) .and. ismissing(phasecmp(:ntim-2)) )
tmp2 = tmp1(0) + ind( (.not. ismissing(phasecmp(tmp1(0):ntim-2))) .and.  \
                      ismissing(phasecmp(tmp1(0)+1:)) )
do n=0, dimsizes(tmp1)-1
  if ( tmp2(n) - tmp1(n) .lt. (30 - ndayavg) ) then
    phasecmp(tmp1(n):tmp2(n)) = missv
  end if
end do
nday_active = num( .not. ismissing(phasecmp) )
print("Ndays of active MJO : "+nday_active+" ("+  \
      (tofloat(nday_active)/(ntim-nWgt/2*2-ndayavg/2*2)*100)+"%)")
;------------------------------------------------------------------------------
system ("rm "+Fout)
nCDF	= addfile (Fout, "c")
nCDF->indx	= indx
nCDF->rmm1	= rmm1
nCDF->rmm2	= rmm2
nCDF->phase	= phase
nCDF->phase_cmp	= phasecmp
;==============================================================================
; plot
;------------------------------------------------------------------------------
print ("plot the data")
wks     = gsn_open_wks ("eps", Fig)
plot    = new (4,graphic)
;------------------------------------------------------------------------------
; EOF plots
;------------------------------------------------------------------------------
res				= True
res@gsnDraw			= False
res@gsnFrame			= False
;------------------------------------------------------------------------------
resP				= True
resP@gsnPaperOrientation	= "portrait"
resP@gsnMaximize		= True
;------------------------------------------------------------------------------
kq0 = ind_nearest_coord(p0,pq,0)
res@vpHeightF = 0.08
res@xyLineThicknessF	= 1.
res@gsnYRefLine			= (/-uqsd05(kq0),uqsd05(kq0)/)
res@gsnYRefLineThicknessF = 0.3
res@gsnAboveYRefLineColor	= (/"transparent","red"/)
res@gsnBelowYRefLineColor	= (/"blue","transparent"/)
;res@xyExplicitLegendLabels	= (/"OLR", "U850", "U200"/)
;res@pmLegendDisplayMode		= "Always"
;res@pmLegendSide		= "Top"
;res@pmLegendParallelPosF        = 0.8
;res@pmLegendOrthogonalPosF      = -0.3
;res@pmLegendWidthF              = 0.2
;res@pmLegendHeightF		= 0.1
;lon = olr&lon
;res@trYMinF			= -0.16
;res@trYMaxF			= 0.16
;res@gsnLeftString  = "EOF "+(n+1)
uq@long_name = (uq@long_name)+" [m s:S:-1:N:]"
print(round(pq(kq0),3)+"-hPa wind")
plot(0) = gsn_csm_xy (wks, yrfracq - yr_sp, uq(kq0,:), res)
;------------------------------------------------------------------------------
delete (res@gsnYRefLine)
delete (res@gsnAboveYRefLineColor)
delete (res@gsnBelowYRefLineColor)
;delete (res@pmLegendDisplayMode)
;delete (res@xyExplicitLegendLabels)
;delete (res@xyLineThicknesses)
;------------------------------------------------------------------------------
; MJO "strong" index
;------------------------------------------------------------------------------
res@tiXAxisString		= "Year"
res@gsnYRefLine			= 1.0
res@gsnYRefLineColor		= "black"
res@xyMonoDashPattern		= True
res@xyLineColors		= (/"black", "blue"/)
res@xyLineThicknesses		= (/1, 2/)
;res@pmLegendDisplayMode		= "Always"
;res@pmLegendParallelPosF        = 0.75
;res@pmLegendWidthF              = 0.25
;res@pmLegendHeightF		= 0.1
;res@xyExplicitLegendLabels = (/"daily", "91-day runavg" /)
res@trYMinF		= 0.
res@trYMaxF		= 5.
res@tfPolyDrawOrder = "Predraw"
;res@gsnLeftString	= "MJO index"
indx@long_name		= "MJO index"
mjo_ind_plt			= new ((/2,ntim/), typeof(indx))
mjo_ind_plt(0,:)		= indx
mjo_ind_plt(1,:)		= (/indx_smt/)
plot(1) = gsn_csm_xy (wks, yrfrac - yr_sp, mjo_ind_plt, res)
dum = new(nyr,graphic)
resPG = True
resPG@gsFillColor = "yellow"
do i=1, nyr-1
  yr0pg = yyyymmdd_to_yyyyfrac( (i    )*10000 + 1201, 0. )
  yr1pg = yyyymmdd_to_yyyyfrac( (i + 1)*10000 + 0228, 0. )
  ypg = (/-1,-1,20,20,-1/)
  xpg = (/yr0pg,yr1pg,yr1pg,yr0pg,yr0pg/)
  dum(i) = gsn_add_polygon(wks,plot(1),xpg,ypg,resPG)
end do
;------------------------------------------------------------------------------
res@xyLineColors		= (/"black", "red"/)
res@xyLineThicknesses		= (/2, 1/)
res@trYMinF		= 0.
res@trYMaxF		= 2.
res@trXMaxF		= 25 - 2.7e-3
sctt = new((/2,nyr/),float)
sctt@long_name		= "MJO index"
do i=1, nyr-1
  djfc = yyyymmdd_to_yyyyfrac( (i + 1)*10000 + 0115, 0. ) + yr_sp
  sctt(0,i) = tofloat(indx_smt({djfc}))
  sctt(1,i) = ( uq(kq0,0+i*12) / uqsd05(kq0) + 3 ) /2.5
end do
print(escorc(sctt(0,:),uq(kq0,0::12)))
plot(2) = gsn_csm_xy (wks, ispan(1,nyr,1), sctt, res)
;------------------------------------------------------------------------------
ress = True
ress@vpHeightF = 0.3
ress@vpWidthF = 0.3
ress@gsnDraw			= False
ress@gsnFrame			= False
ress@xyMarkLineMode       = "Markers"  ; "MarkLines"
ress@xyMarkers            = 16
ress@tiXAxisString		= "MJO index"
;ress@gsnXRefLine			= 1.0
;ress@gsnXRefLineColor		= "black"
ress@xyMonoDashPattern		= True
;ress@xyLineColors		= (/"black", "blue"/)
;ress@xyLineThicknesses		= (/1, 2/)
;ress@pmLegendDisplayMode		= "Always"
;ress@pmLegendParallelPosF        = 0.75
;ress@pmLegendWidthF              = 0.25
;ress@pmLegendHeightF		= 0.1
;ress@xyExplicitLegendLabels = (/"daily", "91-day runavg" /)
;ress@trYMinF		= 0.
;res@trYMaxF		= 5.
;ress@tfPolyDrawOrder = "Predraw"
;ress@gsnLeftString	= "MJO index"
;indx@long_name		= "MJO index"
sctt = new((/2,nyr/),float)
do i=1, nyr-1
  djfc = yyyymmdd_to_yyyyfrac( (i + 1)*10000 + 0115, 0. ) + yr_sp
  sctt(0,i) = tofloat(indx_smt({djfc}))
  sctt(1,i) = uq(kq0,0+i*12)
end do
plot(3) = gsn_csm_xy (wks, sctt(0,:), sctt(1,:), ress)
;------------------------------------------------------------------------------
resP@gsnPanelRowSpec	= True
gsn_panel (wks,plot,(/1,1,1,2/),resP)
;==============================================================================
end

