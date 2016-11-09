load "$UTIL/ncl/header.ncl"

xspace = 80  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 load "./ftn_seasonal.ncl"
 load "./ftn_eof_qbo.ncl"
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
   if (i_des_stations .ne. 0) then
     t1_1 = param0_t1_1
     t2_1 = param0_t2_1
     t1_2 = param0_t1_2
     t2_2 = param0_t2_2
     t1_3 = param0_t1_3
     t2_3 = param0_t2_3
     if ( t2_1 .ne. t1_2 .or. t2_2 .ne. t1_3 ) then
       print("Check param0_t?_?.")
       exit
     end if
   end if
 end if

 if ( tofloat(toint(t1)) .ne. t1 ) then
   print("t1 should be integer.")
   exit
 end if

 ; for the title
 t2t = t2
 if ( tofloat(toint(t2)) .eq. t2 ) then
   l_todec = True
   t2t = t2t - 1
 end if

 wks = cr_wks("fig/qbo_eof_"+t1+"-"+t2t,600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 f1 = addfile(fname,"r")
 u0 = f1->$vname$
 tmpstr = f1->$vname$!0
 p0 = f1->$tmpstr$
 tmpstr = f1->$vname$!1
 t0 = f1->$tmpstr$

 u0@_FillValue = -999.
 t0@_FillValue = -999.

 p = p0
 nz = dimsizes(p)

;--------------------------------------- process 0
 ud0 = u0
 if ( l_deseasonalize ) then
   if (i_des_stations .eq. 0) then
;     ud0 = deseasonalize(u0,i_des_harmonics,0)
     ud0 = deseasonalize(u0,i_des_harmonics,(/t1,t2/))
   else
     ud0 = deseasonalize(u0,i_des_harmonics,(/t1_1,t1_2,t1_3,t2_3/))
   end if
 end if
 if (i_runavg .ne. 0) then
   ud0 = runave_n_Wrap(ud0,i_runavg,0,1)
 end if
 ud = ud0(:,{t1:t2})
 t = t0({t1:t2})

 if (i_runavg .gt. 1) then
   ud(:,0:i_runavg/2-1) = u0@_FillValue
   ud(:,dimsizes(t)-i_runavg/2:) = u0@_FillValue
 end if

 tsn = eofqbo(ud,eofopt,etsopt,True)
 v_eof = tsn@eof

 tsn_s = tsn
 if (shift_minvard .eq. 1) then
   tsn_s = shift_eof(tsn,0.01,3)
;   xy_shift = tsn_s@xy_shift
;   xy_var_d = tsn_s@var_of_d
 end if

;--------------------------------------- cn level
; lcol = (/1,2,4,3/)
; ltck = (/1.5,0.7,0.7,0.7/)

;--------------------------------------- contour
 wi = 25
 he = 30
 xf = 10 + (/0,1/)*(wi+13)
 yf = 10

 xstr = (/"","EOF1"/)
 ystr = (/"p [hPa]","EOF2"/)
 if (shift_minvard .eq. 1) then
   xstr(1) = xstr(1)+" shifted"
   ystr(1) = ystr(1)+" shifted"
 end if

 plot = xylineplot(wks,v_eof,p,"D",(/0,14/),"T",1,"C",1)
 xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",1)
 xyaddline(plot,0,p,"D",0,"T",0.5,"C",1)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,-0.5,0.7,999,-999,False,True)
 axis_str(plot,"",xstr(0),ystr(0))
 setvalues plot
  "tmYLMode" : "Explicit"
  "tmYLValues" : p
  "tmYLLabels" : p
 end setvalues

 draw(plot)

 plot = xylineplot(wks,tsn_s(0,:),tsn_s(1,:),"D",0,"T",0.4,"C",1)
 xyaddmarker(plot,16,"S",0.6,"T",0.6,"C",1)
 xyaddline(plot,(/-100,100/),0,"D",0,"T",0.4,"C",1)
 xyaddline(plot,0,(/-100,100/),"D",0,"T",0.4,"C",1)

 set_size(plot,xf(1),yf,he,he)
 axis_range(plot,-999,999,-999,999,False,False)
 axis_str(plot,"",xstr(1),ystr(1))

 draw(plot)

; lg = lgplot(wks,(/"a","b","c"/),"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
; draw(lg)

 frame(wks)

end

