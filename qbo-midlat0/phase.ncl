load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 casecode = "ra"
; casecode = "hg"

 if (casecode .eq. "ra") then
   load "./param-ra-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ; (/"CGW","CTL"/)
 expname = param0_expname  ;(/"tc" ,"tk" /)
 fdir_u0 = param0_fdir_u0  ; "$DATD/AOL60"+cexp
 fname_u0 = param0_fname_u0  ; "$DATD/AOL60"+cexp
 u0varname = param0_u0varname  ;"u_p_uv"

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 xmin = param0_umin_p1
 xmax = param0_umax_p1
 ymin = param0_umin_p2
 ymax = param0_umax_p2

 if (casecode .eq. "ra") then
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 if (casecode .eq. "hg") then
   wks = cr_wks("fig/phase_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
 else
   if (expname(0) .eq. "era-int_f") then
     wks = cr_wks("fig/phase_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2,600,0)
   else
     wks = cr_wks("fig/phase_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
   end if
 end if
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 u = new((/ne,nt/),float,1.e20)
 ush = new((/ne,nt/),float,1.e20)
 do ie=0, ne-1
   f1 = addfile(fname_u0(ie),"r")
   u(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p1},{-lat0:lat0}) ) )
   ush(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) )
   if (dimsizes(p2) .eq. 2) then
     ush(ie,:) = dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) -  \
                 dim_avg( f1->$u0varname$(it0:it9,{p2(1)},{-lat0:lat0}) )
   end if
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 uc = new((/ne,12/),float,1.e20)
 ushc = new((/ne,12/),float,1.e20)
 if (casecode .eq. "hg") then
   do im=0, 11
     ; use only 50 yrs
     uc(:,im) = dim_avg_n(u(:,im:50*12-1:12), 1)
     ushc(:,im) = dim_avg_n(ush(:,im:50*12-1:12), 1)
   end do
 else
   do im=0, 11
     uc(:,im) = dim_avg_n(u(:,im::12), 1)
     ushc(:,im) = dim_avg_n(ush(:,im::12), 1)
   end do
 end if

 ua = u
 usha = ush
 do n=0, nt-1, 12
 do im=0, 11
   ua(:,n+im) = u(:,n+im) - uc(:,im)
   usha(:,n+im) = ush(:,n+im) - ushc(:,im)
 end do
 end do

;--------------------------------------- cn level
 lcol = (/2,1,3/)

;--------------------------------------- contour
 nl = ne

 mstr = cexp

 wi = 40
 he = 40
 xf = 15 + (/0,1/)*(wi+10)
 yf = 10


do il=0, nl-1

 xstr = "U"+p1+" [m s:S:-1:N:]"
 ystr = (/"U"+p2(0)+" [m s:S:-1:N:]","",""/)
 if (dimsizes(p2) .eq. 2) then
   ystr = (/"U"+p2(0)+" - U"+p2(1)+" [m s:S:-1:N:]","",""/)
 end if

 if (casecode .eq. "ra") then
 if (expname(il) .eq. "era-int_f") then
   xstr = "U"+pp(p1-ppi0)+" [m s:S:-1:N:]"
   ystr = "U"+pp(p2(0)-ppi0)+" [m s:S:-1:N:]"
   if (dimsizes(p2) .eq. 2) then
     ystr = (/"U"+pp(p2(0)-ppi0)+" - U"+pp(p2(1)-ppi0)+" [m s:S:-1:N:]","",""/)
   end if
 end if
 end if

 plot = xylineplot(wks,ua(il,:),usha(il,:),"D",0,"T",0.5,"C",-1)
; do k=0, nz-1
;   xyaddline(plot,t,u(0,k,:),"D",0,"T",1.0,"C",lcol(k))
; end do
 xyaddmarker(plot,16,"S",0.2,"T",1.0,"C",1)

 set_size(plot,xf(il),yf,wi,he)
 axis_range(plot,xmin,xmax,ymin,ymax,False,False)
 axis_str(plot,mstr(il),xstr,ystr(il))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

end do

 frame(wks)

end

