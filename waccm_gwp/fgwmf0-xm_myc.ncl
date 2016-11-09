load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 hdname  = "test0" ;1_0.9"
 expname = "0"
 varname = "TAU"+(/"E","N"/)
 titname = (/"E/W","N/S"/)+"-MF"

 y1 = 1
 y2 = 5
 p = 350 ;350 ;470

;--------------------------------------- workstation
 wks = cr_wks("fig/fgwmf0-xm_myc-"+expname,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 nv = dimsizes(varname)

 fname = new(nt,string)
 n = 0
 do iy=y1, y2
 do im=1, 12
   fname(n) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
              sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",im)+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0),"r")
 y = f0->lat
 ny = dimsizes(y)
 print( f0->ilev({p}) )

 ; get variables
 var = new((/nv,ny,nt/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n),"r")
   tmp = abs( dim_avg_n( dim_avg( f1->$varname(iv)$(:,{p},:,:) ), 0 ) )
   var(iv,:,n) = var(iv,:,n) + tmp(:)*1.e3
 end do
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )

 varss = new((/nv,ny,14/),float,1.e20)
 do n=1, 12
   varss(:,:,n) = dim_avg( var(:,:,ispan(n-1,nt-1,12)) )
 end do
 varss(:,:,0 ) = varss(:,:,12)
 varss(:,:,13) = varss(:,:,1 )

 cosphi = cos(tofloat(y)*3.141592/180.)
 do iv=0, nv-1
 do n=0, 13
   varss(iv,:,n) = varss(iv,:,n) * cosphi(:)
 end do
 end do

;--------------------------------------- cn level
 cnlev = ispan(1,12,1)*0.25
; cnlev(dimsizes(cnlev)-1) = 6
 cnfil = (/5,6,7,8,10,12,14,15,16,17,19,20,22/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 15 + (/0,1/)*(wi+8)
 yf = 10 + (/0,0/)*(he+8)

 xstr  = new(nv,string,"")
 ystr  = new(nv,string,"")
 lbopt = new(nv,string)

 xstr(:) = "month"
 ystr(0) = "latitude (deg)"
 lbopt(:) = ""
 lbopt(0) = "(mPa)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,ispan(0,13,1),y,varss(iv,:,:),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.4
  "pmLabelBarHeightF"        : 0.15
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.03
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,-999,999,False,False)
 axis_str(contour1,titname(iv),xstr(iv),ystr(iv))
 setvalues contour1
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end

