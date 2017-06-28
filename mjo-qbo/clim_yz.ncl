load "$UTIL/ncl/header.ncl"
xspace = 80  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
;axisfontratio = 0.8

load "read_cesm_var.ncl"

begin

;--------------------------------------- parameters
 expname = "f00qbo_uni"
; varname = "ICECLDF"
 varname = "CLDICESTR"
 yr1 = 2           ; the 1st year to analyze
 nyr = 5
 season = "DJF"   ;  DJF / MAM / ... / all (DJF to SON) / ann (J to D)
 lonbnd = (/30,210/)
 latbnd = 30
 zmax = 24

 mon0 = "0001-01"   ; initial month of the simulation

;--------------------------------------- read data
 dirname = "/grace/s2/kyh/cesm_scr/"+expname+"/run"

 yib = latbnd+3
 pib = 1.e3*exp(-(zmax+3)/7.)

 do yr=yr1, yr1+nyr-1
   files = get_file_indx(expname,dirname,"h0",season,yr)
   ti2 = files@t_indx_2
print(files@date_2(0)+"  /  "+files@date_2(1))
   f1 = addfiles(files,"r")
   ssmean = dim_avg_n_Wrap( dim_avg_Wrap(  \
            f1[:]->$varname$(ti2(0):ti2(1),{pib:},{-yib:yib},  \
                             {lonbnd(0):lonbnd(1)}) ), 0 )
   delete( [/ files, f1 /] )
   if (yr .eq. yr1) then
     cmean = ssmean
   else
     cmean = cmean + ssmean
   end if
 end do
 cmean = cmean/tofloat(nyr)
printVarSummary(cmean)

; files = systemfunc("ls "+dirname+"/"+expname+".cam.h0.*.nc")

;--------------------------------------- calculation

;--------------------------------------- cn level
; cnlev = tofloat(ispan(-10,10,1))
 cnlev = tofloat(ispan(1,25,1))
 if (varname .eq. "ICECLDF") then
   cnlev = cnlev*0.02
 else if (varname .eq. "CLDICESTR") then
   cnlev = cnlev*0.75
 else if (varname .eq. "CLOUD") then
   cnlev = cnlev*0.02
 else if ( varname .eq. "QRL" .or. varname .eq. "QRLC" ) then
   cnlev = (1.-cnlev(::-1))*0.15
 else if ( varname .eq. "QRS" .or. varname .eq. "QRSC" ) then
   cnlev = (cnlev-1.)*0.06
; else if (varname .eq. "uw") then
;;   cnlev = cnlev*0.1
;   cnlev = cnlev*0.2
; else if (varname .eq. "D") then
;   cnlev = cnlev*0.1
; end if
; end if
 else
   delete(cnlev)
   cnlev = -999
 end if
 end if
 end if
 end if
 end if

; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = (/0,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2/)
 if ( all(cnlev(:dimsizes(cnlev)-3) .lt. 0.) ) then
   cnfil = cnfil(::-1)
 end if

; lcol = (/2,3,4,5,6/)*0+1

;--------------------------------------- contour1
 wks = cr_wks("fig/clim_yz-"+varname+"-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 cmean&lev = 7.0*log(1.e3/cmean&lev)

 if (cmean@units .eq. "K/s") then
   cmean = cmean*86400.
   cmean@units = "K/day"
 else if (cmean@units .eq. "kg/kg") then
   cmean = cmean*1.e6
   cmean@units = "mg/kg"
 end if
 end if


 wi = 45
 he = 45
 xf = 10 + (/0,0,1,1,1,1/)*(wi+10)
 yf = 15 + (/0,1,0,1,2,3/)*(he+7)

 mstr = cmean@long_name
 xstr = "lat [deg]"
 ystr = "z [km]"
 lbopt = "["+cmean@units+"]"
 if ( cmean@units .eq. "1" .or. cmean@units .eq. "fraction" ) then
   lbopt = " "
 end if
 

 ipl = 0
 contour1 = cnshadeplot_ax(wks,cmean,cnlev,cnfil,lbopt)
; cnaddline_ax(contour1,wks,ul(iv,0,:,:),cnlev(iv,:),0,  \
;           "D",(/0,0,5/),"T",ltck(iv,:),"C",(/1,-1,1/))

 set_size(contour1,xf(ipl),yf(ipl),wi,he)
 axis_range(contour1,-latbnd,latbnd,0,zmax,False,False)
 axis_str(contour1,mstr,xstr,ystr)
; setvalues contour1
;  "tmXBMode" : "Explicit"
;  "tmXBValues" : ispan(0,360,90)
;  "tmXBLabels" : ispan(0,360,90)
;  "tmXBMinorValues" : ispan(0,360,15)
; end setvalues
; if (im .ne. 0) then
;   setvalues contour1
;    "tiMainOn"  : False
;    "tiXAxisOn" : False
;    "tiYAxisOn" : False
;   end setvalues
; end if
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.20
;  "pmLabelBarHeightF"        : 0.04
;  "lbLabelAutoStride" : False
;  "lbLabelFontHeightF" : 0.006
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.015
;  "lbTitleFontHeightF" : 0.008
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 txtstr = season
 dr_txt2(wks,txtstr(0),xf(ipl)+1,yf(ipl)+2,22,5.5,0.,1,1)


 frame(wks)

end

