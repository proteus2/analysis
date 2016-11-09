load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "/data17/kyh/dat/gwjet_wrf"
 var1 = "U"
 lev = 300
 var1name = var1+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_kl2-"+var1name,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/fcoef_klo/fft_"+var1name+"_klo.nc"
 f1 = addfile(fname1,"r")
 nx = 432
 ny = 648
 nt = 432
 dx = 1./10.8
 dy = 1./10.8
 dt = 300.
 v1 = ( f1->$("fc_"+var1name+"_r")$(:,:,:) )^2 + ( f1->$("fc_"+var1name+"_i")$(:,:,:) )^2
 v1(:,:,1:(nx-1)/2) = v1(:,:,1:(nx-1)/2)*2.
 v1 = v1*(nx*dx)*(ny*dy)*(nt*dt)
 nd1 = dimsizes(v1)

;--------------------------------------- process 0
 nk0 = nd1(2)
 nk = nx
 nl = ny
 vs = new((/nl,nk/),float,1.e32)

 no = nt
 kwn = ispan(-nk/2,(nk-1)/2,1)/tofloat(nx*dx)
 lwn = ispan(-nl/2,(nl-1)/2,1)/tofloat(ny*dy)

 vs = 0.
 vs(nl/2:  ,nk/2+1:) = dim_sum_n( v1((no+1)/2:,:(nl-1)/2,1:nk0-2), 0 )
 vs(:nl/2-1,nk/2+1:) = dim_sum_n( v1((no+1)/2:,(nl+1)/2:,1:nk0-2), 0 )
 vs(nl/2:  ,:nk/2-1) = dim_sum_n( v1((no-1)/2:0,nl-1:(nl+1)/2,nk0-1:1), 0 )
 vs(:nl/2-1,:nk/2-1) = dim_sum_n( v1((no-1)/2:0,(nl-1)/2:0   ,nk0-1:1), 0 )
 vs(nl/2:  ,nk/2) = dim_sum_n( v1((no+1)/2:,:(nl-1)/2,0), 0 )
 vs(:nl/2-1,nk/2) = dim_sum_n( v1((no+1)/2:,(nl+1)/2:,0), 0 )
 vs(nl/2:  ,nk/2) = vs(nl/2:  ,nk/2) + dim_sum_n( v1((no-1)/2:0,nl-1:(nl+1)/2,0), 0 )
 vs(:nl/2-1,nk/2) = vs(:nl/2-1,nk/2) + dim_sum_n( v1((no-1)/2:0,(nl-1)/2:0   ,0), 0 )

 vs = vs/(nt*dt)

; xc = ispan(0,nx,1)*dx
; nvc = nv
; nvc(ndim-1) = nvc(ndim-1) + 1
; vc = new(nvc,float)
; vc(:,:,:nx-1) = v0(:,:,:nx-1)
; vc(:,:,nx   ) = v0(:,:,0    )

;--------------------------------------- cn level
 if (var1 .eq. "W") then
   cnlev = fspan(-9.0,-4.0,11)*1.
   kl9 = 5
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-5.0,0.0,11)*1.
   kl9 = 2
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 50
 he = wi
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

 mstr = ""
 xstr = "k [cyc deg:S:-1:N:]"
 ystr = "l [cyc deg:S:-1:N:]"
 lbopt = "log:B:10:N:(PSD[m:S:2:N: s:S:-2:N: / cyc:S:2:N: deg:S:-2:N:])"

do il=0, 0

 contour1 = cnshadeplot(wks,kwn,lwn,log10(vs),cnlev,cnfil,lbopt(il))
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.5
  "pmLabelBarHeightF"        : 0.08
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.015
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-kl9,kl9,-kl9,kl9,False,False)
 axis_str(contour1,mstr,xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

end do

 frame(wks)

end



