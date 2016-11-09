load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

begin

 type = "x11"   ;new window
;type = "eps"    ;eps file
 wks = gsn_open_wks (type, "tx90p")
 gsn_define_colormap (wks, "gui_default")

;---------------------- FILES -----------------------------
fileMT  = systemfunc ("ls /home/yesmuldir/scratch/blocking/extreme_indices/no29_mxtemp.nc")
ft      = addfile (fileMT, "r")
temp    = ft->mx2t
printVarSummary(temp)

fileTH  = systemfunc ("ls /home/yesmuldir/scratch/blocking/extreme_indices/tx90p.nc")
fth     = addfile (fileTH, "r")
thsh    = fth->thsh
printVarSummary(thsh)

fileMon = systemfunc ("ls /home/cyoo/scratch/reanl/ecmwf/monthly/t2/t2.*.nc")
fmt     = addfiles (fileMon, "r")
t2      = fmt[:]->t2
t2      = lonFlip(t2)
printVarSummary(t2)
;==============================================================================
time    = thsh&time
ymd     = cd_calendar (time, 2)
ymdhms  = cd_calendar (time, 0)
year    = floattointeger (ymdhms(:,0))
month   = floattointeger (ymdhms(:,1))
day     = floattointeger (ymdhms(:,2))
ddd     = day_of_year (year, month, day)
;------------------------------------------------------------------------------
dys = where(temp.gt.thsh,1,0)             ;days which is greater than 90 th percentile 
copy_VarCoords(temp,dys)
printVarSummary(dys)
printMinMax(dys,0)
;==============================================================================
nlat = dimsizes(temp&latitude)
mlon = dimsizes(temp&longitude)

yrStrt = 1979
yrLast = 2014

   nyrs = yrLast-yrStrt+1
   nmos = 12

   tsum = new( (/nyrs,nmos,nlat,mlon/), "float")           ;percentage of number of daysfor every month
   do yr=yrStrt,yrLast
      nyr = yr-yrStrt
   do mo=1,nmos
        i = ind(year.eq.yr .and. mo.eq.month)
;        tsum(nyr,mo-1,:,:) = ((dim_sum_n(dys(i,:,:), 0))*100)/dimsizes(dys(i,0,0))  ; --> caution: int/int = int
        tsum(nyr,mo-1,:,:) = (tofloat(dim_sum_n(dys(i,:,:), 0))*100.)/dimsizes(i)
        delete(i)
   end do
   end do

printVarSummary(tsum)
printMinMax(tsum,0)
;========================================================
mon_per = reshape(tsum,(/nyrs*nmos,nlat,mlon/))
copy_VarCoords(t2,mon_per)
printVarSummary(mon_per)
printMinMax(mon_per,0)
;--------------------------------------------------------
;out      = "perc_monTX90p.nc"
;nCDF     = addfile (out, "c")
;nCDF->mon_percTX90  = mon_per
exit
;;------------------------ PLOT -------------------------
	
        res           = True
        res@cnFillOn  = True
        res@vpHeightF  = 0.3
        res@vpWidthF   = 0.75
        res@cnLinesOn = False

        plot = gsn_csm_contour_map_ce(wks,nmb_dys,res) 

end

