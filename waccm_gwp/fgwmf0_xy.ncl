load "$UTIL/ncl/header.ncl"

xspace = 83;170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 4

begin

;--------------------------------------- parameters
 hdname  = "test1_0.9";_500"
 expname = "2";"500"
 varname = "TAU"+(/"E","N"/)
 titname = (/"E/W","N/S"/)+"-MF"

 y1 = 1
 y2 = 5
 p = 350;350 ;470
; monselect = 1  ; 0,1,2,3,4 - Annual, DJF, MAM, JJA, SON
 monselect = loopctl

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

 wks = cr_wks("fig/fgwmf0_xy-"+tstr+"-"+expname,600,0)
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
 x = f0->lon
 y = f0->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 print( f0->ilev({p}) )

 ; get variables
 var = new((/nv,ny,nx/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n),"r")
   tmp = abs( dim_avg_n( f1->$varname(iv)$(:,{p},:,:), 0 ) )
   var(iv,:,:) = var(iv,:,:) + tmp(:,:)
 end do
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)*1.e3

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )
 xa = array_append_record(x,(/360.d0/),0)
 vara = new((/nv,ny,nx+1/),float,1.e20)
 vara(:,:,:nx-1) = var(:,:,:)
 vara(:,:,nx   ) = var(:,:,0)

;--------------------------------------- cn level
; cnlev = fspan(-4,4,21)
; cnlev = (/0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56/)
 cnlev = ispan(0,11,1)*0.3
 cnlev(dimsizes(cnlev)-1) = 6
 cnfil = (/0,5,6,7,8,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 15 + (/0,0/)*(wi+10)
 yf = 10 + (/0,1/)*(he+12)

 lbopt = new(nv,string)
 lbopt(:) = ""
 lbopt(nv-1) = "(mPa)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,xa,y,vara(iv,:,:),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.64
  "pmLabelBarHeightF"        : 0.08
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,-180,180,-999,999,False,False)
 axis_str(contour1,titname(iv),"","")

 cnfinalize_map(contour1,wks,1)

 draw(contour1)

end do


 frame(wks)

end

end do  ; nnn

