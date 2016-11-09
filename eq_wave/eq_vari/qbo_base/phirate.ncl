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

 wks = cr_wks("fig/phirate_"+t1+"-"+t2t,600,0)
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

;--------------------------------------- dump

;--------------------------------------- cn level
; lcol = (/1,2,4,3/)
; ltck = (/1.5,0.7,0.7,0.7/)

;--------------------------------------- contour
 nl = 2

 wi = 60
 he = 40
 xf = 10
 yf = 10 + (/0,1,2/)*(he+6)

 xstr = (/"","progression rate [cyc mon:S:-1:N:]"/)
 ystr = (/"phase [cyc]","month"/)

 vplot = (/phi,tm/)
 ymin = (/-0.54,0/)
 ymax = (/0.54,12/)

do il=0, nl-1

 plot = xylineplot(wks,p_rate,vplot(il,:),"D",0,"T",1,"C",-1)
 xyaddmarker(plot,16,"S",0.4,"T",1.0,"C",1)
 xyaddline(plot,(/-999,999/),0,"D",0,"T",0.5,"C",1)
 xyaddline(plot,0,(/-999,999/),"D",0,"T",0.5,"C",1)

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,-999,0.09,ymin(il),ymax(il),False,False)
 axis_str(plot,"",xstr(il),ystr(il))
 if (il .eq. 1) then
   setvalues plot
    "tmYLMode" : "Explicit"
    "tmYLValues" : ispan(1,12,1) - 0.5
    "tmYLLabels"  : ispan(1,12,1)
   end setvalues
 end if

 draw(plot)

 delete(plot)

end do

; lg = lgplot(wks,(/"a","b","c"/),"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
; draw(lg)

 frame(wks)

end

