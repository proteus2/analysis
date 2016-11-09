load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin
;--------------------------------------- parameters
 cexp    = (/"HadGEM2-CGW","MERRA"/)
 expname = "uj"
 ppid    = "j"
 varname = "u"

 y1 = 1953
 y2 = 2006
 lat0 = 5.
; zs = (/20.95,23.31,26.89,31.8/)
; ps = (/50,30,20,10/)
 ymin = -45
 ymax = 30
 zs = (/16.25,18.59,20.95/)
 ps = (/100,70,50/)
 ymin = -20
 ymax = 15

 opt_ra = 1
 y1r = 1979
 y2r = y2

 wks = cr_wks("fig/u_t4",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

; HadGEM

 fdir  = "$DATD/L60CGW-t"
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$(:,{zs*1.e3},:)
 z = tmp&$tmp!1$ / 1.e3
print(z)
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
 f0 = addfile("$DATD/clim/MERRA/merra.u_ypt."+y1r+"-2006.nc","r")
 latnam = f0->$varname$!2
 lat = f0->$latnam$({-lat0-1.0:lat0+1.0})
 print(lat(0)+" - "+lat(dimsizes(lat)-1))
 wgt = tofloat(lat)*0. + 1.
 wgt(0) = 0.5
 wgt(dimsizes(wgt)-1) = 0.5
 pnam = f0->$varname$!1
; p = tofloat( f0->$pnam$({1:100}) )
 p = tofloat( f0->$pnam$({ps}) )
 nt0 = dimsizes( f0->$varname$&$(f0->$varname$!0)$ )
 u0 = transpose( dim_avg_wgt( \
      tofloat( f0->$varname$(:,{ps},{-lat0-1.0:lat0+1.0}) ), wgt, 0 ) )
 t0 = tofloat(y1r) + fspan(1./24.,nt0/12.-1./24.,nt0)
 zp = 18.6 + 6.6*log(70./p)

 fh = addfile("$DATD/clim/MERRA/merra.h_p."+y1r+"-2012.nc","r")
 zm = fh->h({ps}) / 1.e3

; z0 = zp
 z0 = zm
print(z0)
 nz0 = dimsizes(z0)

end if

;--------------------------------------- process 0
 uss = new((/12,nz/),float)
 do n=0, 11
   uss(n,:) = dim_avg(u(:,ispan(n,nt-1,12)))
 end do
 uss = uss - conform(uss,dim_avg_n(uss,0),1)

 udes = u
 do n=0, nt-1
   udes(:,n) = u(:,n) - uss(mod(n,12),:)
 end do
print(dim_stddev(udes))
;print(dim_stddev(udes(:,nt-nt0:nt-1)))

if (opt_ra .ne. 0) then

; MERRA

 uss0 = new((/12,nz0/),float)
 do n=0, 11
   uss0(n,:) = dim_avg(u0(:,ispan(n,nt0-1,12)))
 end do
 uss0 = uss0 - conform(uss0,dim_avg_n(uss0,0),1)

 udes0 = u0
 do n=0, nt0-1
   udes0(:,n) = u0(:,n) - uss0(mod(n,12),:)
 end do
print(dim_stddev(udes0))

end if

;--------------------------------------- cn level
 lcol = (/1,2,3,4,5,6,7,8,9,10,11,12/)

;--------------------------------------- contour
 yl = 13
 yls = 2

 nfig = 1
 if (opt_ra .ne. 0) then
   nfig = 2
 end if

do ifig=0, nfig-1

 wi = 64
 he = 15
 xf = 17
 yf = 10 + ispan(0,5,1)*(he+5)

 if (ifig .eq. 0) then
   varf = udes
   varx = t
   ys = ispan(y1+yls,y2,yl)
 else
   varf = udes0
   varx = t0
   ys = ispan(y1r+yls,y2r,yl)
   yf = yf + nl*(he+5) + 7
 end if
 ys(0) = ys(0) - yls

 nl = dimsizes(ys)

 wi2 = tofloat(wi)/yl*yls

 xf0 = xf - wi2
 wi0 = wi + wi2

 xstr  = new(nl,string)

 mstr = ":F25:"+(/"(a)","(b)"/)
 mstr(:) = mstr(:)+"   :F22:"+cexp
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = "U (m s:S:-1:N:)"


do il=0, nl-1

 plot = xylineplot(wks,varx,varf,"D",0,"T",1.0,"C",lcol)
 xyaddline(plot,varx,0.,"D",0,"T",0.5,"C",1)

 ye = ys(il)+yl
 xf1 = tofloat(xf)
 wi1 = tofloat(wi)
 if (il .eq. 0) then
   ye = ye + yls
   xf1 = xf0
   wi1 = wi0
 end if
 if (ye .gt. y2+1) then
   ye = y2+1
   wi1 = wi*(tofloat(y2+1-ys(il))/yl)
 end if

 set_size(plot,xf1,yf(il),wi1,he)
 axis_range(plot,ys(il),ye,ymin,ymax,False,False)
 axis_str(plot,"",xstr(il),ystr)

 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues


 if (il .gt. 0) then

  setvalues plot
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  plot2 = xylineplot(wks,varx,varf,"D",0,"T",1.0,"C",lcol)
  xyaddline(plot2,varx,0.,"D",0,"T",0.5,"C",1)

  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),ymin,ymax,False,False)
  axis_str(plot2,"","",ystr)

  setvalues plot2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues

  draw(plot2)

 end if


 draw(plot)

end do  ; il

 if (nfig .gt. 1) then
   dr_txt2(wks,mstr(ifig),xf-wi2-5,yf(0)-4,22,6,0.,"CenterLeft",1)
 end if

 delete(varx)
 delete(varf)
 delete(ys)
 delete(mstr)
 delete(xstr)

end do  ; ifig

 frame(wks)

end

