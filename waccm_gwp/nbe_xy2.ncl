load "$UTIL/ncl/header.ncl"

xspace = 83;170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "merra"
 varname = (/"nbe","curv","adv"/)
 titname = (/"NBE (major)","NBE (major) + Curvature","ADV"/)

 y1 = 2011
 y2 = 2011
 p = 200
; monselect = 1  ; 0,1,2,3,4 - Annual, DJF, MAM, JJA, SON
 mon = (/1/)

;--------------------------------------- workstation
; if (monselect .eq. 0) then
;   mon = ispan(1,12,1)
;   tstr = "Annual"
; else
;   monlist = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
;   tstrlist = (/"DJF","MAM","JJA","SON"/)
;   mon = monlist(monselect-1,:)
;   tstr = tstrlist(monselect-1)
; end if

 wks = cr_wks("fig/nbe_xy2",600,0)
 gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)

 fname = new(nt,string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(n) = "$ANAL/nbe/"+expname+".nbe_xyt."+sprinti("%0.4i",iy)+"."+ \
              sprinti("%0.2i",mon(im))+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0),"r")
 x = f0->lon
 y = f0->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 print( f0->p({p}) )

 ; get variables
 var = new((/nv,ny,nx/),float,1.e20)
 var = 0.
 nt0 = 248
 tmp = new((/nv,nt0,ny,nx/),float,1.e20)
 do n=0, nt-1
   f1 = addfile(fname(n),"r")
   do iv=0, nv-1
     tmp(iv,:,:,:) = f1->$varname(iv)$(:,{p},:,:)
   end do
   var(:,:,:) = var(:,:,:) + dim_avg_n( tmp(:,120:120,:,:), 1 )
;   var(:,:,:) = var(:,:,:) + dim_avg_n( tmp(:,:,:,:), 1 )
;   var(:,:,:) = var(:,:,:) + dim_avg_n( abs(tmp(:,:,:,:)), 1 )
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )
 xa = array_append_record(x,(/360./),0)
 vara = new((/nv,ny,nx+1/),float,1.e20)
 vara(:,:,:nx-1) = var(:,:,:)
 vara(:,:,nx   ) = var(:,:,0)

;--------------------------------------- cn level
 cnlev = ispan(-8,8,1)*1.0;*0.4
 cnfil = (/17,19,21,22,23,24,25,26,0,0,6,5,4,3,2,11,9,7/)
; cnlev = ispan(0,11,1)*0.3
; cnfil = (/0,5,6,7,8,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 15 + (/0,0,0/)*(wi+10)
 yf = 10 + (/0,1,2/)*(he+12)

 lbopt = new(nv,string)
 lbopt(:) = ""
 lbopt(nv-1) = "(10:S:-9:N: s:S:-2:N:)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,xa,y,vara(iv,:,:)*1.e9,cnlev,cnfil,lbopt(iv))
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
 axis_range(contour1,-180,180,-90,90,False,False)
 axis_str(contour1,titname(iv),"","")

 cnfinalize_map(contour1,wks,1)

 draw(contour1)

end do


 frame(wks)

end

