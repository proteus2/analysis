load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 4

begin

;--------------------------------------- parameters
 hdname  = "test"+(/"1","0"/)
 expname = (/"1","0"/)
 varname = "TAU"+(/"E","W","N","S"/)
 titname = (/"East","West","North","South"/)+"ward-MF"

 y1 = 1
 y2 = 2
; monselect = 1  ; 0,1,2,3,4 - Annual, DJF, MAM, JJA, SON
 monselect = loopctl

 p = (/350.,1./)

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

 wks = cr_wks("fig/fgwmf-xm_y2-"+tstr+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)
 ne = dimsizes(expname)

 fname = new((/ne,nt/),string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(:,n) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
                sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 ; get axis
 f0 = addfile(fname(0,0),"r")
 y = f0->lat
 nz = dimsizes(p)
 ny = dimsizes(y)
 z = 7.0*log(1.e3/p)

 ; get variables
 var = new((/ne,nv,nz,ny/),float,1.e20)
 var = 0.
 do ie=0, ne-1
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(ie,n),"r")
   tmp = abs( dim_avg_n( dim_avg( f1->$varname(iv)$(:,{p},:,:) ), 0 ) )
   var(ie,iv,:,:) = var(ie,iv,:,:) + tmp(:,:)
 end do
 end do
 end do
 var(:,:,:,:) = var(:,:,:,:)/tofloat(nt)*1.e3

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 var = where( var .eq. 0., var@_FillValue, var )

 varp = new((/nv,ne*nz,ny/),float,1.e20)
 cnt = 0
 do ie=0, ne-1
 do k=0, nz-1
   varp(:,cnt,:) = var(ie,:,k,:)
   cnt = cnt + 1
 end do
 end do

;--------------------------------------- cn level
 ldsh = (/14,0,14,0/)
 ltck = (/1.6,1.6,0.5,0.5/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 15 + (/0,1,2,3/)*(wi+8)
 yf = 10 + (/0,0,0,0/)*(he+8)

 xstr  = new(nv,string,"")
 ystr  = new(nv,string,"")

 xstr(:) = "latitude (deg)"
 ystr(0) = "M (mPa)"

do iv=0, nv-1

 plot = xylineplot(wks,y,varp(iv,:,:),"D",ldsh,"T",ltck,"C",1)

 set_size(plot,xf(iv),yf(iv),wi,he)
 axis_range(plot,-999,999,-999,999,False,False)
 axis_str(plot,titname(iv),xstr(iv),ystr(iv))

 draw(plot)

end do

 lg = lgplot(wks,p+" hPa","D",ldsh,"T",1.0,"C",1,2)
 set_size_lg(lg,xf(0)+wi-13.0,yf(0)+0.1,wi*0.25,8)
 setvalues lg
;  "lgLabelOffsetF" : 0.2
  "lgLineDashSegLenF"  : 0.026
  "lgLabelFontHeightF" : 0.005
 end setvalues
 draw(lg)


 frame(wks)

end

end do  ; nnn

