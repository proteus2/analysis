load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 4

begin

;--------------------------------------- parameters
 hdname  = (/"test1_0.9","test0"/)
 expname = (/"2"        ,"0"    /)
 varname = (/"V"/)
 titname = (/"COR"/)

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
   monlist = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/), (/4,5,6/) /)
   tstrlist = (/"DJF","MAM","JJA","SON","AMJ"/)
   mon = monlist(monselect-1,:)
   tstr = tstrlist(monselect-1)
 end if

 wks = cr_wks("fig/corudiff-xm_yz-"+tstr+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)

 fname = new((/nt,2/),string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(n,:) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
                sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0,0),"r")
 p = f0->lev
 y = f0->lat
 nz = dimsizes(p)
 ny = dimsizes(y)
 z = 7.0*log(1.e3/p)

 ; get variables
 var = new((/nv,nz,ny/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n,0),"r")
   f2 = addfile(fname(n,1),"r")
   tmp = dim_avg_n( dim_avg( f1->$varname(iv)$(:,:,:,:) ), 0 ) - \
         dim_avg_n( dim_avg( f2->$varname(iv)$(:,:,:,:) ), 0 )
   var(iv,:,:) = var(iv,:,:) + tmp(:,:)
 end do
 end do
 var(:,:,:) = var(:,:,:)/tofloat(nt)*86400.

 fcor = 2.*7.292116e-5*sin(tofloat(y)*3.141592/180.)
 do iv=0, nv-1
 do k=0, nz-1
   var(iv,k,:) = var(iv,k,:) * fcor(:)
 end do
 end do

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = fspan(-4,4,21)
; cnlev = (/-80,-60,-40,-20,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,20,40,60,80/)
 cnlev = (/-1000,-500,-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200,500,1000/)*0.1
 cnfil = (/17,19,21,22,23,24,25,26,0,0,6,5,4,3,2,11,9,7/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 15 + (/0,1,2/)*(wi+8)
 yf = 10 + (/0,0,0/)*(he+8)

 xstr  = new(nv,string,"")
 ystr  = new(nv,string,"")
 lbopt = new(nv,string)

 xstr(:) = "latitude (deg)"
 ystr(0) = "height (km)"
 lbopt(:) = ""
 lbopt(nv/2) = "(m s:S:-1:N: day:S:-1:N:)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,y,z,var(iv,:,:),cnlev,cnfil,lbopt(iv))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.0
  "pmLabelBarHeightF"        : 0.15
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.04
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

