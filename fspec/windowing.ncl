load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 cmax = 25
 ctck = 5
 
 wks = cr_wks("windowing",600,0)
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- process 0

 nmon_patch = 1
 nhour = 4

 nmon = 1 + 2*nmon_patch
 nt = nmon*30*nhour

 wgt = new(nt,float)

 npad = nmon_patch*30*nhour
 pi05 = 2.*atan(1.)

 wgt = 1.
 do i=1, npad
   wgt(i-1) = sin(pi05*tofloat(i)/tofloat(30*nhour+1))
 end do
 wgt(nt-1:nt-npad) = wgt(0:npad-1)

 nday = (/31,31,28,31,30,31/)
 nn1 = sum(nday(0:2))*nhour
 it1 = (nn1 - nt)/2 + 0
 nn2 = sum(nday(1:3))*nhour
 it2 = (nn2 - nt)/2 + nday(0)*nhour
 nn3 = sum(nday(2:4))*nhour
 it3 = (nn3 - nt)/2 + sum(nday(0:1))*nhour

 t = ispan(0,sum(nday)*nhour-1,1)*(1./tofloat(nhour))

;--------------------------------------- cn level
 xtck = (/0,15,nday(0),nday(0)+15,sum(nday(0:1)),sum(nday(0:1))+15,  \
          sum(nday(0:2)),sum(nday(0:2))+15,  \
          sum(nday(0:3)),sum(nday(0:3))+15,sum(nday(0:4))/)
 xlab = (/"Dec 1","16","Jan 1","16","Feb 1","16","Mar 1","16","Apr 1","16","May 1"/)

;--------------------------------------- contour
 wi = 60
 he = 15
 xf = 10
 yf = 10

 plt = xylineplot(wks,t(it1:it1+nt-1),wgt,"D",0,"T",0.8,"C",1)
 xyaddline(plt,t(it2:it2+nt-1),wgt,"D",0,"T",0.8,"C",2)
 xyaddline(plt,t(it3:it3+nt-1),wgt,"D",0,"T",0.8,"C",3)

 set_size(plt,xf,yf,wi,he)
 axis_range(plt,-1,sum(nday(0:4))+1,0,1.1,False,False)
 axis_str(plt,"","","")
 setvalues plt
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtck
  "tmXBLabels" : xlab
 end setvalues
 draw(plt)

 frame(wks)

end

