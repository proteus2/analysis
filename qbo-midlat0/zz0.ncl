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
;   load "./param-ra2-0.ncl"
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

 nma = param0_nma
 nph = param0_nph

 xmin = -28 ;-20 ;param0_umin_p1
 xmax = 20 ;14 ;param0_umax_p1
 ymin = -30 ;-32 ;param0_umin_p2
 ymax = 34 ;30 ;param0_umax_p2

 if (casecode .eq. "hg") then
   wks = cr_wks("fig/phase-"+nma+"m"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
 else
   if (expname(0) .eq. "era-int_f") then
     wks = cr_wks("fig/phase-"+nma+"m"+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2,600,0)
   else
;     wks = cr_wks("fig/phase-"+nma+"m"+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
     wks = cr_wks("zz0",600,0)
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
   if (cexp .eq. "qbo_obs") then
;     tmp0 = f1->$u0varname$
     tmp0 = wgt_runave_Wrap( f1->$u0varname$, (/1.,2.,1./), 0 )
     u(ie,:) = tmp0({p1},it0:it9)
     ush(ie,:) = tmp0({p2(0)},it0:it9)
     if (dimsizes(p2) .eq. 2) then
       ush(ie,:) = tmp0({p2(0)},it0:it9) - tmp0({p2(1)},it0:it9)
     end if
   else
     u(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p1},{-lat0:lat0}) ) )
     ush(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) )
     if (dimsizes(p2) .eq. 2) then
       ush(ie,:) = dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) -  \
                   dim_avg( f1->$u0varname$(it0:it9,{p2(1)},{-lat0:lat0}) )
     end if
   end if
 end do
 
 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmp = u
   u(:,0) = tmp(:,nt-1)
   u(:,1:) = tmp(:,:nt-2)
   tmp = ush
   ush(:,0) = tmp(:,nt-1)
   ush(:,1:) = tmp(:,:nt-2)
   delete(tmp)
 end if

 uc = new((/ne,12/),float,1.e20)
 ushc = new((/ne,12/),float,1.e20)
 if (casecode .eq. "hg") then
   do im=0, 11
     ; use only 50 yrs
     uc(:,im) = dim_avg_n(u(:,im:50*12-1:12), 1)
     ushc(:,im) = dim_avg_n(ush(:,im:50*12-1:12), 1)
   end do
   if ( nma .eq. 3 .and. nt .gt. 50*12 ) then
     uc(:,0) = dim_avg_n(u(:,12:50*12:12), 1)
     ushc(:,0) = dim_avg_n(ush(:,12:50*12:12), 1)
   end if
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
 
print(min(ua))
print(max(ua))
print(min(usha))
print(max(usha))


 x = ua - 0.5*tofloat(xmin+xmax)
 y = (usha - 0.5*(ymin+ymax))*(tofloat(xmax-xmin)/tofloat(ymax-ymin))

 deg = where(x .ne. 0., atan(y/x)*180./3.141592, 90.)
 deg = where(abs(y) .gt. abs(x) .and. x*y .gt. 0.,  \
             90. - atan(x/y)*180./3.141592, deg)
 deg = where(abs(y) .gt. abs(x) .and. x*y .lt. 0.,  \
             -(90. - (atan(abs(x/y))*180./3.141592)), deg)
 deg = where(x .eq. 0. .and. y .lt. 0., -90., deg)
 deg = where(x .lt. 0., deg+180., deg)      ; -90 <= d < 270
 deg = 90. - deg                            ; -180 < d <= 180 (clockwise)
 deg = where(deg .lt. 0., deg+360., deg)    ;   0 <= d < 360

 deg00 = 45.
 deg0 = fspan(deg00,deg00+360.-(360./tofloat(nph)),nph)

 phs = where(deg .lt. deg0(0) .or. deg .ge. deg0(nph-1), 1, 999)
 do i=0, nph-2
   phs = where(deg .ge. deg0(i) .and. deg .lt. deg0(i+1), i+2, phs)
 end do

;--------------------------------------- cn level
 lcol = (/2,1,3/)

;--------------------------------------- contour
 nl = ne

 xstr = new(2,string)
 ystr = new(12,string)
 xstr(:) = ""
 ystr(:) = ""

 mstr = (/"",""/)

 wi = 30
 he = 30
 xf = 15 + ispan(0,11,1)*tofloat(wi+8)/nma
 yf = 10 + ispan(0,1,1)*(he+15)
 col = (/2,3,4,6,7,8,13,5/)


do im=0, 11, nma
do il=0, nl-1

 xstr(il) = "U"+p1+" [m s:S:-1:N:]"
 ystr(0) = "U"+p2(0)+" [m s:S:-1:N:]"
 if (dimsizes(p2) .eq. 2) then
   ystr(0) = "U"+p2(0)+" - U"+p2(1)+" [m s:S:-1:N:]"
 end if

 if (casecode .eq. "ra") then
 if (expname(il) .eq. "era-int_f") then
   xstr(il) = "U"+pp(p1-ppi0)+" [m s:S:-1:N:]"
   if (im .eq. 0) then
     ystr(0) = "U"+pp(p2(0)-ppi0)+" [m s:S:-1:N:]"
     if (dimsizes(p2) .eq. 2) then
       ystr(0) = (/"U"+pp(p2(0)-ppi0)+" - U"+pp(p2(1)-ppi0)+" [m s:S:-1:N:]","",""/)
     end if
   end if
 end if
 end if

 ua_pl   = new((/nma,(nt/12)/),float)
 usha_pl = new((/nma,(nt/12)/),float)
 do nn=0, nma-1
   ua_pl  (nn,:) = ua  (il,im+nn::12)
   usha_pl(nn,:) = usha(il,im+nn::12)
 end do
 plot = xylineplot(wks,ua_pl,usha_pl,"D",0,"T",0.5,"C",-1)
; do k=0, nz-1
;   xyaddline(plot,t,u(0,k,:),"D",0,"T",1.0,"C",lcol(k))
; end do
 xyaddmarker(plot,16,"S",0.2,"T",1.0,"C",1)

 dumm = new(nt/(12/nma),graphic)
 plres = True
 plres@gsMarkerIndex = 16
 plres@gsMarkerSizeF = 0.001
 cnt = 0
 do n=im, nt-1, 12
   do nn=0, nma-1
     plres@gsMarkerColor = col(phs(il,n+nn)-1)
     dumm(cnt+nn) = gsn_add_polymarker(wks,plot,ua(il,n+nn),usha(il,n+nn),plres)
   end do
   cnt = cnt + nma
 end do

 set_size(plot,xf(im),yf(il),wi,he)
 axis_range(plot,xmin,xmax,ymin,ymax,False,False)
 axis_str(plot,mstr(il),xstr(il),ystr(im))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
  "tmYLMode" : "Manual"
  "tmYLTickSpacingF"  : 10
  "tmYLMinorPerMajor" : 4
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

delete(plot)
delete(dumm)
end do
end do

 frame(wks)

end

