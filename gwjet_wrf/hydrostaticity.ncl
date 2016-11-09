load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "/data18/GW_jet_wrf/dat"
 var1 = "PRS"
 var2 = "P_HYD"
 lev = 300
do it=0, 16
; it = 16
 var1name = "P"
 var2name = "P_HYD"
;	Time = UNLIMITED ; // (1152 currently)
 
;--------------------------------------- read data
 fname1 = fdir+"/x6/"+var1+"/"+var1+"__x6__08_00_00__12_00_00__6hrly.nc"
 fname2 = fdir+"/x6/"+var2+"/"+var2+"__x6__08_00_00__12_00_00__6hrly.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 p  = f1->$var1name$(it,:,53:52+ny,:)
 ph = f2->$var2name$(it,:,53:52+ny,:)

;--------------------------------------- process 0
 nd = dimsizes(p)
 nz = nd(0)

 dp = p
 dph = ph
 do k=0, nz-2
   dp (k,:,:) = p (k,:,:) - p (k+1,:,:)
   dph(k,:,:) = ph(k,:,:) - ph(k+1,:,:)
 end do
 dp (nz-1,:,:) = dp (nz-2,:,:)
 dph(nz-1,:,:) = dph(nz-2,:,:)

 hyd = abs(dp - dph)/abs(dph)

 print("Maximum nonhydrostaticity: "+max(hyd))
end do

 exit

 wks = cr_wks("fig/hydrostaticity",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*0.01
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 12 + (/0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1,2,2,2,2/)*(he+6)

 mstr  = new(nz,string)
 xstr  = new(nz,string)
 ystr  = new(nz,string)
 lbopt = new(nz,string)

 mstr(:) = ""
 xstr(:) = ""
 xstr(nz/2:) = "lon [deg]"
 ystr(:) = ""
 ystr(ispan(0,nz-1,nz/2)) = "lat [deg]"
 lbopt(:) = ""
 lbopt(nz-1) = " "

do il=0, nz/4-2

 contour1 = cnshadeplot(wks,x,y,hyd(il*4,:,:),cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,x,y,hyd(il*4,:,:),cnlev,1,"D",0,"T",0.5,"C",1)
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

end do

 frame(wks)

end



