load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 0.8
begin

;--------------------------------------- parameters
; expname = "f2000_cntl"
; expname = "f2000_qbo_cyclic"
; expcode = "exp0"
 expname = "f2000_qbo_uni_test"
 expcode = "exp1"
 lata = 10.
 levs = (/25.,55.,70.,90.,100./)
 

 wks = cr_wks("fig/qbo-u_t-"+expcode,600,0)
 gsn_define_colormap(wks,"default")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 fdir = "/grace/s2/kyh/cesm_scr/"+expname+"/run"
if (expcode .eq. "exp0") then
 fdir = "~/tmp/"+expname+"/run"
end if
 files = systemfunc("ls "+fdir+"/"+expname+".cam.h0.*.nc")

 nz = dimsizes(levs)

 f1 = addfiles(files,"r")
 t = f1[:]->time
; date = cd_calendar(t,-3)
 mon = ispan(1,dimsizes(t),1)-0.5
 p_m = round(f1[0]->lev({levs})*10.,3)/10.
;print( f1[0]->lat({-lata:lata}) )
;print( p_m )
 
 u  = transpose( dim_avg_Wrap( dim_avg_Wrap( f1[:]->U      (:,{levs},{-lata:lata},:) ) ) )
 ur = transpose( dim_avg_Wrap( dim_avg_Wrap( f1[:]->QBO_U0 (:,{levs},{-lata:lata},:) ) ) )
 fu = transpose( dim_avg_Wrap( dim_avg_Wrap( f1[:]->QBOTEND(:,{levs},{-lata:lata},:) ) ) )
;printVarSummary(u)

 fu = fu*86400.*10.
if (expcode .eq. "exp0") then
 fu = fu*3.
end if

;--------------------------------------- cn level
 lcol = (/1,2,3,6,4,7/)
 ltck = 0.8

;--------------------------------------- plot
 wi = 60
 he = 12
 xf = 15 + ispan(0,nz-1,1)*(wi+5)*0
 yf = 15 + ispan(0,nz-1,1)*(he+4)

 xstr0 = "month"
 mstr = ""
 xstr = where(yf .eq. max(yf), xstr0, "")
 ystr = "U [m s:S:-1:N:]"

do k=0, nz-1

 plot1 = xylineplot(wks,mon,(/ u(k,:), ur(k,:) /),"D",0,"T",ltck,"C",lcol(:1))
 xyaddline(plot1,mon,fu(k,:),"D",0,"T",0.5,"C",colind)
 xyaddline(plot1,mon,0,"D",0,"T",0.5,"C",colind)

 set_size(plot1,xf(k),yf(k),wi,he)
 axis_range(plot1,-999,120,-999,999,False,False)
 axis_str(plot1,mstr,xstr(k),ystr)
 setvalues plot1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF" : 12
 end setvalues

 draw(plot1)

 dr_txt(wks,p_m(k)+" hPa",xf(k)+5,yf(k)+1.5,21,4,0.)

end do  ; k

 frame(wks)

end

