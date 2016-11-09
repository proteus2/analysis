load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 lat0 = 25.
 z0    = 18.59  ; 18.59  ; 17.41
 ymax0 = 0.4    ; 0.4    ; 0.6

 wks = cr_wks("fig/fz-rw_y-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")
print(fnam2+"")

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nt = 647 + 1
   nw = 8

   rw_i = new((/nt-1,nw,ny/),float)
   rw_i(:,:,:) = f1->f_z(:,:nw-1,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     f2->f_z(:,8:9,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(:,2,{z0*1.e3},{-lat0:lat0}) + \
                      f2->f_z(:,3,{z0*1.e3},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(:,13,{z0*1.e3},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(:,16:17,{z0*1.e3},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,nw,ny/),float)
   rw0(:nt-2,:,:) = rw_i(:,:,:)
   rw0(nt-1,:,:) = rw0(nt-2,:,:)

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   yp = dim_avg_n( rw(:,:,:,:), 1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 ldsh = (/0,16,14/)

;--------------------------------------- contour
 wname = (/ (/"K","E-MRG","E-IG"/), (/"R","W-MRG","W-IG"/) /)
 wi = 30
 he = 22
 xf = 12 + (/0,1/)*38
 yf = 10

 xstr = "lat (deg)"
 ystr = (/":F10:F:B:z:NF22: (mPa)",""/)
 ymax = ymax0*(/-1,1/)

do idr=0, 1

 plot = xylineplot(wks,lat1,yp(:,idr,:),"D",ldsh,"T",1.5,"C",1)
 set_size(plot,xf(idr),yf,wi,he)
 axis_range(plot,-lat0,lat0,0,ymax(idr),False,False)
 axis_str(plot,"",xstr,ystr(idr))

 draw(plot)

 lg = lgplot(wks,wname(idr,:),"D",ldsh,"T",1.5,"C",1,2)
 set_size_lg(lg,xf(idr),yf+he+8,wi*0.6,12)
 draw(lg)

end do  ; idr


 frame(wks)

end



