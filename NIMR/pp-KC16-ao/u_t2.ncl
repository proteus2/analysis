load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin
;--------------------------------------- parameters
 cexp    = "CGW"
 expname = "tc"
 ppid    = "j"
 varname = "u"

 y1 = 1955
 y2 = y1 - 1 + 14 ;40
 ydif = 100
 lat0 = 5.
 zs = (/22,28,33/)

 wks = cr_wks("fig/u_t2-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir  = "$DATD/AOL60"+cexp
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$
 z0 = tmp&$tmp!1$ / 1.e3
 iz = ind_nearest_coord(zs,z0,0)
 z = z0(iz)
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

 u = transpose( dim_avg( tmp(:,iz,{-lat0:lat0}) ) )

 delete(tmp)

;--------------------------------------- process 0
 u!1 = "t"
 u&t = t

 us = new((/2,nz,12/),float)
 do n=0, 11
   us(0,:,n) = dim_avg( u(:,{y1+tofloat(n)/12.:y2:12}) )
   us(1,:,n) = dim_avg( u(:,{y1+ydif+tofloat(n)/12.:y2+ydif:12}) )
 end do
 usa = us - conform_dims(dimsizes(us),dim_avg(us),(/0,1/))

print(usa(0,2,:)+"  "+usa(1,2,:))

;--------------------------------------- line/shading att.
 lcol = ispan(1,nz,1)

;--------------------------------------- plot
 yl = 13
 yls = 1;2

 wi = 64
 he = 15
 xf = 20
 yf = 10 + ispan(0,5,1)*(he+5)

 nfig = 2

do ifig=0, nfig-1

 if (ifig .gt. 0) then
   y1 = y1 + ydif
   y2 = y2 + ydif
   yf = yf + (he+5)*3 + 4
yf = yf - (he+5)*2
 end if

 ys = ispan(y1+yls,y2,yl)
 ys(0) = ys(0) - yls

 nl = dimsizes(ys)

 wi2 = tofloat(wi)/yl*yls

 xf0 = xf - wi2
 wi0 = wi + wi2

 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = ":F25:"+(/"(a)","(b)"/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = "U [m s:S:-1:N:]"
 lbopt(:) = ""
 if (ifig .eq. nfig-1) then
   lbopt(nl-1) = "[m s:S:-1:N:]"
 end if


do il=0, nl-1

 plot1 = xylineplot(wks,t,u,"D",0,"T",0.8,"C",lcol)
 xyaddline(plot1,t,0,"D",0,"T",0.5,"C",colind)

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

 set_size(plot1,xf1,yf(il),wi1,he)
 axis_range(plot1,ys(il),ye,-999,999,False,False)
 axis_str(plot1,"",xstr(il),ystr)

 setvalues plot1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues


 if (il .gt. 0) then

  setvalues plot1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  plot2 = xylineplot(wks,t,u,"D",0,"T",0.8,"C",lcol)
 xyaddline(plot2,t,0,"D",0,"T",0.5,"C",colind)

  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),-999,999,False,False)
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


 draw(plot1)

end do  ; il

 if (nfig .gt. 1) then
   dr_txt2(wks,mstr(ifig),xf-wi2-11,yf(0)-1,22,6,0.,"CenterLeft",1)
 end if

end do  ; ifig

 frame(wks)

end

