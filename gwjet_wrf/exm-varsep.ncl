load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9
sz_mstr = p2n(5)

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf/x6_5min"
; fdname = "domain_x6.nc"
 ftail = "__x6__08_00_00__12_00_00__5min.nc"
 varn = "W"
 vard = varn
; vard = "PRS"
 lev = 8000 ;5000
 tims = (/6.5,7.0/)
 ymin = 40
 ymax = 70
 xmin = 10
 xmax = 50
 
 wks = cr_wks("fig/exm-varsep",600,0)
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.,1.,0.)

;--------------------------------------- main body
 nv = 3
 nt = dimsizes(tims)

 fnames = fdir+"/"+vard+"_z/prt_d300km/"+ \
          (/"mean","prt"/)+"_"+vard+"__z"+sprinti("%5.5i",lev)+ftail
; fnames(0) = fnames(1)  ; for 1-hrly data
 
 fd = addfile(fnames(0),"r")
 x = fd->lon - 40
 y = fd->lat
 nx = dimsizes(x)
 ny = dimsizes(y)
 dx = x(1) - x(0)
 dy = y(1) - y(0)
 nt0 = dimsizes( fd->t )

 tims = round(tims*1.e4,0)*1.e-4
 time = round(fspan(4,8,nt0)*1.e4,0)*1.e-4  ; [day]
 it = get1Dindex(time,tims)
print("No. of lon, lat, t :  "+nx+", "+ny+", "+nt+"/"+nt0)

 nxp = nx*2

do n=0, nt-1  ; global loop
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; read -----------------------------------------------------
 fm = addfile(fnames(0),"r")
 fp = addfile(fnames(1),"r")

 vm = fm->$("mean_"+varn)$(it(n),:,:)
 vp = fp->$("prt_"+varn)$(it(n),:,:)
 vs = (/vm+vp,vm,vp/)

 if (varn .eq. "W") then
   vs = vs * 1.e3
 end if 
 if (varn .eq. "P") then
   vs = vs / 100.
 end if
 if (varn .eq. "T") then
   vs(0:1,:,:) = vs(0:1,:,:) + 300.
 end if

; process --------------------------------------------------
 xpl = new(nxp,float,1.e20)
 xpl(:nx-1) = x(:)
 xpl(nx:) = xpl(nx-1) + (x(1)-x(0))*ispan(1,nxp-nx,1)

 vpl = new((/nv,ny,nxp/),float,1.e20)
 do j=0, ny-1
   vpl(:,j,0:nx-1) = vs(:,j,:)
   vpl(:,j,0+nx:nx+nx-1) = vs(:,j,:)
 end do

; plot setting ---------------------------------------------
 mstr = ""
 xstr = ""
 ystr = ""
 tday = toint(tims(n)+1.e-4)
 tmp = (tims(n)+1.e-4 - tday)*24.
 thour = toint(tmp)
 tmin = toint((tmp - thour)*60.)
 mstr = "Day "+tday+",  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)
 if ((nt - n) .le. 4) then
   xstr = "longitude [deg]"
 end if
 if (mod(n,4) .eq. 0) then
   ystr = "latitude [deg]"
 end if

 lbopt = ""
 if (n .eq. 0) then
   lbopt = "[mm s:S:-1:N:]"
 end if

; U
; cnlev0 = ispan(-5,5,1)*5 + 25
; cnlev2 = ispan(-5,5,1)*0.5
; T
; cnlev0 = ispan(-5,5,1)*1 + 330
; cnlev2 = ispan(-5,5,1)*0.2
; cnlev1 = cnlev0
; W
 cnlev0 = ispan(-5,5,1)*3
 cnlev2 = ispan(-5,5,1)*3
 cnlev1 = cnlev0
 
; cnfil = ispan(6,61,5) + (/-4,-4,-4,-3,-2,-1,1,2,3,4,4,4/)
 cnfil = array_append_record( ispan(2,32,6), ispan(35,65,6), 0 )
 nfil = dimsizes(cnfil)
 cnfil(nfil/2-1:nfil/2) = 0

; plot -----------------------------------------------------
 pscl = 0.8 ;0.45
 rlonlat = 0.82  ; ~55 deg  ; 0.71  ; ~45 deg  ; 0.64  ; ~50 deg
 wi = (xmax - xmin)*rlonlat*pscl
 he = (ymax - ymin)*pscl
 xf = 10 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+5)
 yf = 20 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+9)

 contour1 = cnshadeplot(wks,xpl,y,vpl(2,:,:),cnlev2,cnfil,lbopt)
; cnaddline(contour1,wks,xpl,y,vpl(0,:,:),cnlev0,0, \
;           "D",0,"T",0.4,"C",1)
 cnaddline(contour1,wks,xpl,y,vpl(1,:,:),cnlev1,0, \
           "D",(/0,0,16/),"T",0.8,"C",(/1,-1,1/))
;           "D",0,"T",1.0,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.1
  "lbTitlePosition"  : "Right"
  "lbTitleDirection" : "Across"
  "lbLabelFontHeightF" : 0.04
  "lbTitleOffsetF"     : 0.04
 end setvalues
 set_size(contour1,xf(n),yf(n),wi,he)
 axis_range(contour1,xmin,xmax,ymin,ymax,False,False)
 axis_str(contour1,mstr,xstr,ystr)
 setvalues contour1
  "tmYLMode" : "Manual"
  "tmYLTickSpacingF" : 10
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end do  ; global loop

frame(wks)

end

