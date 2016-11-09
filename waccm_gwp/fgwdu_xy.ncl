load "$UTIL/ncl/header.ncl"

xspace = 83;170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1,1 ;3,3

begin

;--------------------------------------- parameters
 hdname  = "test0";_500"
 expname = "0";"500"
 varname = (/"UTGWSPEC"/)
 titname = (/"GWDx'"/)

 y1 = 1
 y2 = 5
 p = 0.02 ;470
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

 wks = cr_wks("fig/fgwdu_xy-p"+p+"-"+tstr+"-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
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
   tmp = dim_avg_n( f1->$varname(iv)$(:,{p},:,:) - \
                    conform_dims((/nv,ny,nx/),dim_avg(f1->$varname(iv)$(:,{p},:,:)),(/0,1/)), 0 )
   var(iv,:,:) = var(iv,:,:) + tmp(:,:)
 end do
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)*86400.

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )
 xa = array_append_record(x,(/360.d0/),0)
 vara = new((/nv,ny,nx+1/),float,1.e20)
 vara(:,:,:nx-1) = var(:,:,:)
 vara(:,:,nx   ) = var(:,:,0)

;--------------------------------------- cn level
 cnlev = fspan(-8,8,17)
 cnfil = (/17,19,21,22,23,24,25,26,0,0,6,5,4,3,2,11,9,7/)

;--------------------------------------- contour
 wi = 60
 he = 15
 xf = 15 + (/0,1,0,1/)*(wi+10)
 yf = 10 + (/0,0,1,1/)*(he+12)

 lbopt = new(nv,string)
 lbopt(:) = ""
 lbopt(nv-1) = "(10 m s:S:-1:N: day:S:-1:N:)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,xa,y,vara(iv,:,:)*1.e-1,cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12*2
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
 axis_range(contour1,-180,180,0,999,False,False)
 axis_str(contour1,titname(iv),"","")

 cnfinalize_map(contour1,wks,1)

 draw(contour1)

end do


 frame(wks)

end

end do  ; nnn

