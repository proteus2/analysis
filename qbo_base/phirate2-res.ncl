load "$UTIL/ncl/header.ncl"

xspace = 80  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 load "./param0.ncl"
 load "./param1_eof.ncl"

 eofopt = param1_eofopt
 etsopt = param1_etsopt
 shift_minvard = param1_shift_minvard

 fname = param0_fname
 vname = param0_vname

 t1 = param0_t1
 t2 = param0_t2

 l_deseasonalize = param0_l_deseasonalize

 i_runavg = param0_i_runavg

 if ( l_deseasonalize ) then
   i_des_harmonics = param0_i_des_harmonics
   i_des_stations = param0_i_des_stations
 end if

 if ( tofloat(toint(t1)) .ne. t1 ) then
   print("t1 should be integer.")
   exit
 end if

 ; for the title
 t2t = t2
 if ( tofloat(toint(t2)) .eq. t2 ) then
   t2t = t2t - 1
 end if

 wks = cr_wks("fig/phirate2-res_"+t1+"-"+t2t,600,0)
 gsn_define_colormap(wks,"BlGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = "phsamp_"+t1+"-"+t2t
 if ( l_deseasonalize ) then
   fname1 = fname1+"_d"+i_des_harmonics
   if (i_des_stations .ne. 0) then
     fname1 = fname1+"s"
   end if
 end if
 fname1 = fname1+"_m"+i_runavg
 fname1 = fname1+"_e"+eofopt+""+etsopt+""+shift_minvard

 if ( .not. isfilepresent("dat/"+fname1+".nc") ) then
   print("Excute 'qbo_eof-ampphs.ncl' first.")
   exit
 end if
 
 print("Reading "+"dat/"+fname1+".nc")
 f1 = addfile("dat/"+fname1+".nc","r")

 phi = f1->phase
 amp = f1->amplitude
 p_rate = f1->p_rate
 aphi = f1->a_phase
 t = f1->t

;--------------------------------------- process 0
 phi = where(phi .ge. 0.5, phi-1., phi)
 tm = (t - toint(t))*12

 tmy = (t - toint(t))

 pi = 4.*atan(1.)

 rftn0 = new((/5,dimsizes(t)/),float,p_rate@_FillValue)
 rftn0(0,:) = 1.
 rftn0(1,:) = cos(phi*(2.*pi))
 rftn0(2,:) = sin(phi*(2.*pi))
 rftn0(3,:) = cos(tmy*(2.*pi))
 rftn0(4,:) = sin(tmy*(2.*pi))

 rcoef0 = reg_multlin(p_rate,rftn0,False)

 c0 = rcoef0(0)
 cphi = sqrt(rcoef0(1)^2 + rcoef0(2)^2)
 ctmy = sqrt(rcoef0(3)^2 + rcoef0(4)^2)
 phi0r = atan2(rcoef0(2),rcoef0(1))
 tmy0r = atan2(rcoef0(4),rcoef0(3))
 phi0 = phi0r*180./pi
 tmy0 = tmy0r*180./pi
print(c0+"  "+cphi+"  "+ctmy+"  "+phi0+"  "+tmy0)

 p_rate_f = c0 + cphi*cos(phi*(2.*pi)-phi0r) +  \
                 ctmy*cos(tmy*(2.*pi)-tmy0r)

 ntb = 121
 tmyb = fspan(0,1,ntb)
 npb = 217  ; 201
 phib = fspan(-0.54,0.54,npb)

 p_rate_f_b = c0 +  \
     conform_dims((/ntb,npb/),cphi*cos(phib*(2.*pi)-phi0r),1) +  \
     conform_dims((/ntb,npb/),ctmy*cos(tmyb*(2.*pi)-tmy0r),0)

 pr_res = p_rate_f - p_rate

;--------------------------------------- dump

;--------------------------------------- cn level
 maxpr = 0.1
 ncolormap = 200
 
 cnlev = fspan(0,maxpr,ncolormap+1)
 cnfil = ispan(1,ncolormap+2,1)
 cnfil(0) = colind
 cnfil(ncolormap+1) = 1

 mcol = new(dimsizes(p_rate),integer,-1)
; mcol0 = toint(abs(pr_res)/maxpr*ncolormap)+2
 mcol0f = toint(p_rate_f/maxpr*ncolormap)+2
 mcol = toint(p_rate/maxpr*ncolormap)+2
 mcol = where(mcol .gt. ncolormap+1, 1, mcol)  ; black for p_rate > maxpr
 mcol = where(mcol .lt. 2, colind, mcol)       ; grey for p_rate < 0

 mcol = where(abs(pr_res) .lt. 0.01, mcol0f, mcol)

;--------------------------------------- contour
 nl = 1

 wi = 60
 he = 40
 xf = 10
 yf = 10 + (/0,1,2/)*(he+6)

 xstr = (/"phase [cyc]"/)
 ystr = (/"month"/)
 lbopt = (/"[cyc mon:S:-1:N:]"/)

 ymin = (/0/)
 ymax = (/12/)

do il=0, nl-1

; plot = xylineplot(wks,phi,tm,"D",0,"T",1,"C",-1)
;; xyaddmarker(plot,16,"S",0.4,"T",1.0,"C",1)
;; xyaddline(plot,0,(/-999,999/),"D",0,"T",0.5,"C",1)

 plot = cnshadeplot(wks,phib,tmyb*12.,p_rate_f_b,cnlev,cnfil,lbopt(il))
 setvalues plot
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.15
  "lbLabelFontHeightF" : 0.012
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.016
  "lbBoxLinesOn"  : False
  "lbLabelStride" : ncolormap/10
 end setvalues

 gsres0 = True
 gsres0@gsMarkerColor = 0
 dum0 = new(dimsizes(p_rate),graphic)
 do n=0, dimsizes(p_rate)-1
   if ( .not. ismissing(p_rate(n)) ) then
     gsres0@gsMarkerIndex = 1  ; 1 / 16
     gsres0@gsMarkerSizeF = max((/abs(pr_res(n)),0.005/))*1.0
     dum0(n) = gsn_add_polymarker(wks,plot,phi(n),tm(n),gsres0)
   end if
 end do

 gsres = True
 dum = new(dimsizes(p_rate),graphic)
 do n=0, dimsizes(p_rate)-1
   if ( .not. ismissing(p_rate(n)) ) then
     gsres@gsMarkerIndex = 1  ; 1 / 16
     gsres@gsMarkerColor = mcol(n)
     gsres@gsMarkerSizeF = max((/abs(pr_res(n)),0.005/))*1.0 - 0.003
     dum(n) = gsn_add_polymarker(wks,plot,phi(n),tm(n),gsres)
   end if
 end do

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,-999,999,ymin(il),ymax(il),False,False)
 axis_str(plot,"",xstr(il),ystr(il))
 if (il .eq. 0) then
   setvalues plot
    "tmYLMode" : "Explicit"
    "tmYLValues" : ispan(1,12,1) - 0.5
    "tmYLLabels"  : ispan(1,12,1)
   end setvalues
 end if

 draw(plot)

 delete(dum0)
 delete(dum)
 delete(plot)

end do

; lg = lgplot(wks,(/"a","b","c"/),"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
; draw(lg)

 frame(wks)

end

