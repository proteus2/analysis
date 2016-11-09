load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"

 y1 = 1952
 y2 = 2015
 lat0 = 5.
 z0 = 21.0  ; p0 = 50.

 wks = cr_wks("fig/u_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new((/ne,nt/),string)
 fdir  = "../dat/L60CGW-t/xm"
 do ie=0, ne-1
   n = 0
   do iy=y1, y2
   do im=1, 12
     fname(ie,n) = fdir+"/xm.uan"+expname(ie)+"a.p"+ppid(ie)+"_"+iy+sprinti("%0.2i",im)+".nc"
     n = n + 1
   end do
   end do
 end do
 fexist = isfilepresent(fname)

 n1 = new(ne,integer)
 n2 = new(ne,integer)
 n1 = -1
 n2 = nt-1
 do ie=0, ne-1
   do n=0, nt-1
     if ( fexist(ie,n) ) then
       if (n1(ie) .eq. -1) then
         n1(ie) = n
       end if
     else
       if ( n1(ie) .ne. -1 .and. n2(ie) .eq. nt-1 ) then
         n2(ie) = n - 1
       end if
     end if
   end do
   if (n1(ie) .eq. -1) then
     n2(ie) = -1
   end if
 end do
print(n1/12.+"  "+n2/12.)

 ; get z
 ie0 = maxind(n1)
 f0 = addfile(fname(ie0,n1(ie0)),"r")
 tmp = f0->$varname(ie0)$
 zname = tmp!1
 z = tmp&$zname$ / 1.e3
 nz = dimsizes(z)
 delete(tmp)

 u = new((/ne,nz,nt/),float,1.e20)
 do ie=0, ne-1
   if (n1(ie) .ne. -1) then
print(expname(ie)+"")
     f1 = addfiles(fname(ie,n1(ie):n2(ie)),"r")
     ListSetType(f1,"cat")
     tmp = dim_avg( f1[:]->$varname(ie)$(:,:,{-lat0:lat0}) )
     u(ie,:,n1(ie):n2(ie)) = transpose( tmp )
     delete(tmp)
   end if
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

;--------------------------------------- process 0
 do k=0, nz-1
   if (abs(z(k)-z0)/z0 .le. 0.02) then
     kk1 = k
   end if
 end do

 ts = new((/ne,100/),float,1.e20)
 zz = new((/ne,100/),float,1.e20)
 do ie=0, ne-1
   cnt = 0
   do n=0, nt-2
     if ( ismissing(u(ie,kk1,n+1)) ) then
       break
     end if
     if ( u(ie,kk1,n)*u(ie,kk1,n+1) .lt. 0. ) then
       x = ( abs(u(ie,kk1,n))*t(n+1) + abs(u(ie,kk1,n+1))*t(n) )/ \
           (abs(u(ie,kk1,n))+abs(u(ie,kk1,n+1)))
       ts(ie,cnt) = x
       zz(ie,cnt) = z(kk1)
       cnt = cnt + 1
     end if
   end do
 end do

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 nl = 3
 ys = (/1953,1971,1989/)
 yl = 18
 yls = 2

 wi = 126
 he = 20
 xf = 34
 yf = 10 + (/0,1,2/)*28

 wi2 = tofloat(wi)/yl*yls

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:) = ""
; mstr(0) = cexp(0)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "z (km)"
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N:)"

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,z,u(0,:,:),cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,t,z,u(0,:,:),(/0./),0,"D",0,"T",0.8,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50+0.05
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.4
  "pmLabelBarHeightF"        : 0.15
  "lbLabelFontHeightF" : 0.009
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.012
 end setvalues

 set_size(contour1,xf,yf(il),wi,he)
 axis_range(contour1,ys(il),ys(il)+yl,16,45,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

if (il .gt. 0) then
  setvalues contour1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contour2 = cnshadeplot(wks,t,z,u(0,:,:),cnlev,cnfil,"")
  cnaddline(contour2,wks,t,z,u(0,:,:),(/0./),0,"D",0,"T",0.8,"C",1)
  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),16,45,False,False)
  axis_str(contour2,"","",ystr(il))
  setvalues contour2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  cnfinalize(contour2,wks)
  draw(contour2)
end if

 cnfinalize(contour1,wks)

 draw(contour1)


 plot = xylineplot(wks,ts(0,:),zz(0,:),"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot,16,"S",0.8,"T",1.0,"C",1)

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,ys(il),ys(il)+yl,16,45,False,False)
 axis_str(plot,"","","")
 setvalues plot
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues

if (il .gt. 0) then
  plot2 = xylineplot(wks,ts(0,:),zz(0,:),"D",0,"T",1.0,"C",-1)
  xyaddmarker(plot2,16,"S",0.8,"T",1.0,"C",1)
  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),16,45,False,False)
  axis_str(plot2,"","","")
  setvalues plot2
   "tmXBOn" : False
   "tmXTOn" : False
   "tmYLOn" : False
   "tmYROn" : False
  end setvalues
  getvalues plot
   "xyCoordDataSpec" : spec
  end getvalues
  getvalues spec
   "xyMarkerSizeF" : msize
  end getvalues
  getvalues plot2
   "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
   "xyMarkerSizeF" : msize
  end setvalues
;  draw(plot2)
end if

; draw(plot)

end do


 frame(wks)

end



