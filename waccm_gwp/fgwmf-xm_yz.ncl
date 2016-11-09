load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 4

begin

;--------------------------------------- parameters
 hdname  = "test1_0.9"
 expname = "2"
 varname = "TAU"+(/"E","W","N","S"/)
 titname = (/"East","West","North","South"/)+"ward-MF"

 y1 = 1
 y2 = 5
; monselect = 1  ; 0,1,2,3,4 - Annual, DJF, MAM, JJA, SON
 monselect = loopctl

 zmin = 0
 zmax = 120

;--------------------------------------- workstation
 if (monselect .eq. 0) then
   mon = ispan(1,12,1)
   tstr = "Annual"
 else
   monlist = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
   tstrlist = (/"DJF","MAM","JJA","SON"/)
   mon = monlist(monselect-1,:)
   tstr = tstrlist(monselect-1)
 end if

 wks = cr_wks("fig/fgwmf-xm_yz-"+tstr+"-"+expname,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)

 fname = new(nt,string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(n) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
              sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0),"r")
 p = f0->ilev
 y = f0->lat
 nz = dimsizes(p)
 ny = dimsizes(y)
 z = 7.0*log(1.e3/p)

 ; get variables
 var = new((/nv,nz,ny/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n),"r")
   tmp = abs( dim_avg_n( dim_avg( f1->$varname(iv)$(:,:,:,:) ), 0 ) )
   var(iv,:,:) = var(iv,:,:) + tmp(:,:)
 end do
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)*1.e3

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )

;--------------------------------------- cn level
; cnlev = (/0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56/)
 cnlev = ispan(-8,2,1)*0.2
 cnfil = (/5,6,7,8,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 15 + (/0,1,2,3/)*(wi+8)
 yf = 10 + (/0,0,0,0/)*(he+8)

 xstr  = new(nv,string,"")
 ystr  = new(nv,string,"")
 lbopt = new(nv,string)

 xstr(:) = "latitude (deg)"
 ystr(0) = "height (km)"
 lbopt(:) = ""
 lbopt(nv/2-1) = "(log:B:10:N:M [mPa])"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,y,z,log10(var(iv,:,:)),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.2
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6
  "pmLabelBarHeightF"        : 0.15
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,-999,999,zmin,zmax,False,False)
 axis_str(contour1,titname(iv),xstr(iv),ystr(iv))

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end

end do  ; nnn

