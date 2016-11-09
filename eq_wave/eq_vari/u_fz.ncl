load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = "era-int"
 dirname = "ERA-I"

 lat0 = 20.
 ymax0 = 0.5

 x_wmrgw = 1 ;3

 wks = cr_wks("fig/u_fz",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "$DATD/clim/"+dirname(0)+"/"+expname(0)+".u_ypt.nc"

   f0 = addfile(fnam0,"r")
   lat = f0->g0_lat_2({-lat0:lat0})
   p = f0->lv_ISBL1

   ny = dimsizes(lat)
   nz = dimsizes(p)
   nt = 408

   u = tofloat( f0->U_GDS0_ISBL(:nt-1,:,{-lat0:lat0}) )

;--------------------------------------- process 0
   uss = new((/12,nz,ny/),float,1.e20)
   uss!2 = "lat"
   uss&lat = lat
   do n=0, 11
     uss(n,:,:) = dim_avg_n(u(ispan(n,nt-1,12),:,:),0)
   end do

   uds = new((/nt,nz,ny/),float,1.e20)
   uds!2 = "lat"
   uds&lat = lat
   do n=0, 11
   do nn=n, nt-1, 12
     uds(nn,:,:) = u(nn,:,:) - uss(n,:,:)
   end do
   end do

   uds_h = new((/nt,nz,ny/),float,1.e20)
   uds_h!2 = "lat"
   uds_h&lat = lat
   do k=0, nz-1
   do j=0, ny-1
     cf = ezfftf(uds(:,k,j))
     cf(:,:nt/12-1) = 0.
     uds_h(:,k,j) = ezfftb(cf,0.)
   end do
   end do

   uds_l = uds - uds_h

 fuss = dim_stddev_n( uss(:,:,:), 0 )
 fuds = dim_stddev_n( uds(:,:,:), 0 )
 fuds_h = dim_stddev_n( uds_h(:,:,:), 0 )
 fuds_l = dim_stddev_n( uds_l(:,:,:), 0 )

 ff = (/ fuss, fuds, fuds_h, fuds_l /)

;--------------------------------------- cn level
 cnlev = ispan(1,24,1)*0.5
 cnfil = ispan(26,2,1)

;--------------------------------------- contour
 wi = 40
 he = 40
 xf = 20+(/0,1,2,3/)*(wi+8)
 yf = 10

 mstr = (/"U_s","U_des","U_des_h","U_des_l"/)
 xstr = "latitude (deg)"
 ystr = (/"pressure (hPa)","","",""/)
 lbopt = "(m s:S:-1:N:)"
; txtstr = "U"

do ic=0, 3

 contour1 = cnshadeplot(wks,lat,p,ff(ic,:,:),cnlev,cnfil,lbopt)
; contour1 = cnlineplot(wks,t,zu,u0,cnlev,0)
; cnline_att(contour1,"D",(/0,0,1/),"T",0.8,"C",1)
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.12 ;+ 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.60
  "pmLabelBarHeightF"        : 0.10
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.016
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(ic),yf,wi,he)
 axis_range(contour1,-999,999,500,50,False,True)
 axis_str(contour1,mstr(ic),xstr,ystr(ic))

 cnfinalize(contour1,wks)

 draw(contour1)

; dr_txt2(wks,txtstr,xf+1,yf+he-2,22,6,0.,"CenterLeft",1)

end do  ; ic

 frame(wks)

end



