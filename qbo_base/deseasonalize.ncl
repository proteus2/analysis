load "$UTIL/ncl/header.ncl"

xspace = 160  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 load "./ftn_seasonal.ncl"
 load "./param0.ncl"

 fname = param0_fname
 vname = param0_vname

 t1 = param0_t1
 t2 = param0_t2

 i_des_harmonics = param0_i_des_harmonics
 i_des_stations = param0_i_des_stations

 i_runavg = param0_i_runavg

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

 wks = cr_wks("fig/deseasonalized_"+t1+"-"+t2t,600,0)
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
 if (i_des_stations .eq. 0) then
;   ud0 = deseasonalize(u0,i_des_harmonics,0)
   ud0 = deseasonalize(u0,i_des_harmonics,(/t1,t2/))
 else
   ud0 = deseasonalize(u0,i_des_harmonics,(/t1_1,t1_2,t1_3,t2_3/))
 end if
 if (i_runavg .ne. 0) then
   ud0 = runave_n_Wrap(ud0,i_runavg,0,1)
 end if
 ud = ud0(:,{t1:t2})
 t = t0({t1:t2})

;--------------------------------------- cn level
; lcol = (/1,2,4,3/)
; ltck = (/1.5,0.7,0.7,0.7/)

;--------------------------------------- contour
 nl = nz

 wi = 150
 he = 10
 xf = 10 + (/0,1,2,3/)*(wi+10)
 yf = 10 + ispan(nl-1,0,1)*(he+5)

 mstr = new(nl,string)
 mstr = ""
 mstr(nl-1) = "Deseasonalized"
 xstr = new(nl,string)
 xstr = ""
 xstr(0) = "year"
 ystr = "U [m s:S:-1:N:]"

do il=0, nl-1

 plot = xylineplot(wks,t,ud(il,:),"D",0,"T",1,"C",1)
; xyaddmarker(plot,16,"S",0.3,"T",1.0,"C",1)
 xyaddline(plot,t,0,"D",0,"T",0.5,"C",1)

 set_size(plot,xf(0),yf(il),wi,he)
 axis_range(plot,t1,t2,-999,999,False,False)
 axis_str(plot,mstr(il),xstr(il),ystr)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickStartF" : 1950
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues

 dr_txt(wks,p(il)+" hPa",xf(0)+5,yf(il)+he-1.5,22,5,0.)

 draw(plot)

; lg = lgplot(wks,(/"a","b","c"/),"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
; draw(lg)

end do  ; il

 frame(wks)

end

