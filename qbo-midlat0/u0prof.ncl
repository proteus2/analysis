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

 nma = 12
 nph = param0_nph

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
   wks = cr_wks("fig/u0prof-"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
 else
   if (expname(0) .eq. "era-int_f") then
     wks = cr_wks("fig/u0prof-"+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2,600,0)
   else
     wks = cr_wks("fig/u0prof-"+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2,600,0)
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

 axisname = f1->$u0varname$!1
 pver = tofloat( f1->$axisname$ )
 if (expname(0) .eq. "era-int_f") then
   pver_all = asciiread("/data11/data-arch/ERA-I-nr/plev/plev_ps1000_ERA-I",-1,"float")
   pver = pver_all(toint(pver)-1)
 end if
 nz = dimsizes( pver )
 up = new((/ne,nt,nz/),float,1.e20)
 do ie=0, ne-1
   f1 = addfile(fname_u0(ie),"r")
   up(ie,:,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,:,{-lat0:lat0}) ) )
 end do

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmp = u
   u(:,0) = tmp(:,nt-1)
   u(:,1:) = tmp(:,:nt-2)
   tmp = ush
   ush(:,0) = tmp(:,nt-1)
   ush(:,1:) = tmp(:,:nt-2)
   tmpv = up
   up(:,0,:) = tmpv(:,nt-1,:)
   up(:,1:,:) = tmpv(:,:nt-2,:)
   delete(tmp)
   delete(tmpv)
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


 var8m = new((/ne,12/nma,nph,nz/),float)
 navg = new((/ne,12/nma,nph/),integer)
 navg = 0
 tmp = new((/nz,nt/12*nma/),float)
 do ie=0, ne-1
 do im=0, 12/nma-1
 do ip=0, nph-1
   tmp = 0.
   do ma=0, nma-1
   do iy=0, nt/12-1
     imon = iy*12+im*nma+ma
     if (phs(ie,imon) .eq. ip+1) then
       tmp(:,iy*nma+ma) = up(ie,imon,:)
       navg(ie,im,ip) = navg(ie,im,ip) + 1
     end if
   end do
   end do
   if (navg(ie,im,ip) .ne. 0) then
     var8m(ie,im,ip,:) = dim_sum(tmp)/tofloat(navg(ie,im,ip))
   end if
 end do
 end do
 end do

;--------------------------------------- cn level
 lcol = (/2,1,3/)

;--------------------------------------- contour
 nl = ne

 xstr = "U [m s:S:-1:N:]"
 ystr = "p [hPa]"

 wi = 30
 he = 30
 xf = 15 + ispan(0,11,1)*tofloat(wi+8)/nma
 yf = 10 + ispan(0,1,1)*(he+15)
 col = (/2,3,4,6,7,8,13,5/)


do im=0, 11, nma
do il=0, nl-1

 plot = xylineplot(wks,var8m(il,im,:,:),pver,"D",0,"T",0.5,"C",col(:nph-1))
; do k=0, nz-1
;   xyaddline(plot,t,u(0,k,:),"D",0,"T",1.0,"C",lcol(k))
; end do

 set_size(plot,xf(im),yf(il),wi,he)
 axis_range(plot,-40,20,200,5,False,True)
 axis_str(plot,"",xstr,ystr)

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

delete(plot)
end do
end do

 frame(wks)

end

