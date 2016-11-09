load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "/data18/GW_jet_wrf/dat"
 var1 = "W"
 var2 = "U"
 lev = 300
 nave = 12*12
 n1s = ispan(nave/2-1,1151,nave)
; n1s = ispan(nave-1,1151,nave)
 var1name = var1+"_"+lev
 var2name = var2+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/field2_xy-"+var1name+"-"+var2name,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/"+var1+"/"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname2 = fdir+"/x1/"+var2+"/"+var2name+"__x1__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 v1  = f1->$var1name$(n1s,53:52+ny,:)
 v2i = f2->$var2name$(n1s,53:52+ny,:)

;--------------------------------------- process 0
 nv = dimsizes(v1)
 nt = nv(0)
 ndim = dimsizes(nv)

 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

 if ( var2 .eq. "U" ) then
   v2 = 0.5*(v2i(:,:,0:nx-1) + v2i(:,:,1:nx))
 end if
 if ( var2 .eq. "V" ) then
   v2 = 0.5*(v2i(:,0:ny-1,:) + v2i(:,1:ny,:))
 end if

; xc = ispan(0,nx,1)*dx
; nvc = nv
; nvc(ndim-1) = nvc(ndim-1) + 1
; vc = new(nvc,float)
; vc(:,:,:nx-1) = v0(:,:,:nx-1)
; vc(:,:,nx   ) = v0(:,:,0    )

;--------------------------------------- cn level
 cnlev1 = fspan(-10,10,21)*3.
 cnlev2 = fspan(-10,10,21)*4.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

 mstr  = new(nt,string)
 xstr  = new(nt,string)
 ystr  = new(nt,string)
 lbopt = new(nt,string)

 mstr(:) = ""
 xstr(:) = ""
 xstr(nt/2:) = "lon [deg]"
 ystr(:) = ""
 ystr(ispan(0,nt-1,nt/2)) = "lat [deg]"
 lbopt(:) = ""
 lbopt(nt-1) = "w [mm s:S:-1:N:]"
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

do il=0, nt-1

 contour1 = cnshadeplot(wks,x,y,v1(il,:,:)*1.e3,cnlev1,cnfil,lbopt(il))
 cnaddline(contour1,wks,x,y,v2(il,:,:),cnlev2,1,"D",0,"T",0.5,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,40,25,75,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; dr_txt(wks,subl(il),xf(0)-8,yf(il),25,6,0.)

end do

 frame(wks)

end



