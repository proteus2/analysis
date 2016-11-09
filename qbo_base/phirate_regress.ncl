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

 wks = cr_wks("fig/phirate_regress_"+t1+"-"+t2t,600,0)
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
 p_rate = f1->p_rate
 t = f1->t

;--------------------------------------- process 0
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

;--------------------------------------- dump

;--------------------------------------- cn level
; lcol = (/1,2,4,3/)
; ltck = (/1.5,0.7,0.7,0.7/)

;--------------------------------------- contour
 nl = 1

 wi = 2 * (toint(t2+1-1.e-4)-t1)
 he = 15
 xf = 10
 yf = 10 + (/0,1,2/)*(he+4)

 xstr = (/"year"/)
 ystr = (/"prog. rate [cyc mon:S:-1:N:]"/)

 vplot = (/p_rate,p_rate_f/)

do il=0, nl-1

 plot = xylineplot(wks,t,vplot,"D",0,"T",(/1.0,0.6/),"C",(/1,2/))
; xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
 xyaddline(plot,t,0,"D",0,"T",0.3,"C",1)

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,t1,toint(t2+1-1.e-4),-999,0.09,False,False)
 axis_str(plot,"",xstr(il),ystr(il))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickStartF" : 1950
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues

 draw(plot)

end do

; lg = lgplot(wks,(/"a","b","c"/),"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
; draw(lg)

 frame(wks)

end

