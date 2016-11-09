load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin
;--------------------------------------- parameters
 cexp    = (/"CGW","MERRA"/)
 expname = "tc"
 ppid    = "j"
 varname = "u"

 y1 = 1980 ;2058 ;1960
 y2 = 2097 ;2097 ;1999
 lat0 = 10.
 zdot = 21.0  ; p0 = 50.

 opt_ra = 0
 y1r = 1979
 y2r = y2

 wks = cr_wks("fig/u_t-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

; HadGEM

 fdir  = "$DATD/AOL60"+cexp(0)
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$
 z = tmp&$tmp!1$ / 1.e3
 nz = dimsizes(z)
 t = tmp&$tmp!0$
 t_char = tochar(t@units)
 if ( t_char(0:10) .ne. "days since " ) then
   print(t@units)
   print("Check this. It should be 'days since ...'.")
   exit
 end if
 y0i = toint(tostring(t_char(11:14)))
 m0i = toint(tostring(t_char(16:17)))
 d0i = toint(tostring(t_char(19:20)))
 print("date0 : "+y0i+" "+m0i+" "+d0i)

 t = y0i + (m0i-1)/12. + (d0i-1)/360. + t/360.
 nt = dimsizes(t)
 print("data : "+t(0)+" - "+t(nt-1))

 u = transpose( dim_avg( tmp(:,:,{-lat0:lat0}) ) )

 delete(tmp)

; t = y1 + (fspan(1,nt,nt)-0.5)/12.
; u = u - conform_dims( dimsizes(u), dim_avg(u), (/0,1/) )

if (opt_ra .ne. 0) then

; MERRA

 print("MERRA")
 f0 = addfile("$DATD/clim/MERRA/merra.u_ypt."+y1r+"-2012.nc","r")
 latnam = f0->$varname$!2
 lat = f0->$latnam$({-lat0-1.0:lat0+1.0})
 print(lat(0)+" - "+lat(dimsizes(lat)-1))
 wgt = tofloat(lat)*0. + 1.
 wgt(0) = 0.5
 wgt(dimsizes(wgt)-1) = 0.5
 pnam = f0->$varname$!1
 p = tofloat( f0->$pnam$({1:100}) )
 nt0 = dimsizes( f0->$varname$&$(f0->$varname$!0)$ )
 u0 = transpose( dim_avg_wgt( \
      tofloat( f0->$varname$(:,{1:100},{-lat0-1.0:lat0+1.0}) ), wgt, 0 ) )
 t0 = tofloat(y1r) + fspan(1./24.,nt0/12.-1./24.,nt0)
 zp = 18.6 + 6.6*log(70./p)

 fh = addfile("$DATD/clim/MERRA/merra.h_p."+y1r+"-2012.nc","r")
 zm = fh->h({1:100}) / 1.e3

; z0 = zp
 z0 = zm
 nz0 = dimsizes(z0)

end if

;--------------------------------------- process 0
 do k=0, nz-1
   if (abs(z(k)-zdot)/zdot .le. 0.02) then
     kk1 = k
   end if
 end do

 ts = new(300,float,1.e20)
 zz = new(300,float,1.e20)
 cnt = 0
 do n=0, nt-2
   if ( ismissing(u(kk1,n+1)) ) then
     break
   end if
   if ( u(kk1,n)*u(kk1,n+1) .lt. 0. ) then
     x = ( abs(u(kk1,n))*t(n+1) + abs(u(kk1,n+1))*t(n) )/ \
         (abs(u(kk1,n))+abs(u(kk1,n+1)))
     ts(cnt) = x
     zz(cnt) = z(kk1)
     cnt = cnt + 1
   end if
 end do

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnlev = fspan(-40,28,18)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 nfig = 1
 if (opt_ra .ne. 0) then
   nfig = 2
 end if

do ifig=0, nfig-1

 wi = 64
 he = 15
 xf = 20
 yf = 10 + ispan(0,5,1)*(he+5)

 if (ifig .eq. 0) then
   varf = runave(u,5,0)
   varx = t
   vary = z
 else
   varf = runave(u0,5,0)
   varx = t0
   vary = z0
 end if

 mstr = ":F25:"+(/"(a)","(b)"/)
 mstr(:) = mstr(:)+"   :F22:"+cexp
 xstr = "year"
 ystr = "height (km)"

print(vary)

 plot1 = xylineplot(wks,varx,varf(28,:),"D",0,"T",0.6,"C",1)
 xyaddline(plot1,varx,varf(29,:),"D",0,"T",0.6,"C",2)

 set_size(plot1,xf,yf,wi,he)
 axis_range(plot1,y1,y2,-999,999,False,False)
 axis_str(plot1,"","","")

 draw(plot1)

 if (nfig .gt. 1) then
   dr_txt2(wks,mstr(ifig),xf-wi-5,yf(0)-4,22,6,0.,"CenterLeft",1)
 end if

 delete(varx)
 delete(vary)
 delete(varf)

end do  ; ifig

 frame(wks)

end

