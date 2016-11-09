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
 z0    = 10.1 ;18.59

 ; wave type
 idr = 0
 iv = 0

 wks = cr_wks("fig/fz-k_ty-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"
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


   fnam0 = "$DATD/L60CGW/xm.uan"+expname(0)+"a.pj.nc"
   f0 = addfile(fnam0,"r")
   u = f0->u(:,{10.e3:17.41e3},{-lat0:lat0})
;   u = f0->u(:,{fspan(3,16,10)*1.e3},{-lat0:lat0})
   nz = dimsizes(u&$u!1$)
print(u&$u!1$/1.e3)

   t = fspan(0,13,14)

   wss = new((/ny,14/),float)
   do n=1, 12
     wss(:,n) = dim_avg_n( rw(iv,ispan(n-1,nt-1,12),idr,:), 0 )
   end do
   wss(:,0 ) = wss(:,12)
   wss(:,13) = wss(:,1 )
   wss!0 = "lat"
   wss&lat = lat1

   uss = new((/nz,ny,14/),float)
   do n=1, 12
     uss(:,:,n) = dim_avg_n( u(ispan(n-1,nt-1,12),:,:), 0 )
   end do
   uss(:,:,0 ) = uss(:,:,12)
   uss(:,:,13) = uss(:,:,1 )
   uss!1 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0
   wss = wss*(-1.)

;--------------------------------------- cn level
 cnlev = (/0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10/)*0.1
 cnfil = (/5,6,7,8,10,12,14,15,16,17,18,19,20,21/)

 cnlevu = ispan(-30,30,2)
 ldsh = (/0,16,14/)

;--------------------------------------- contour
 wi = 30
 he = 22
 xf = 12 + (/0,1,0,1,0,1,0,1,0,1,0,1/)*38
 yf = 10 + (/0,0,1,1,2,2,3,3,4,4,5,5/)*28

 xstr  = new(nz,string,"")
 ystr  = new(nz,string,"")
 lbopt = new(nz,string)
 xstr(max((/nz-2,0/)):) = "month"
 ystr(ispan(0,nz-1,2)) = "lat (deg)"
 lbopt = ""
 lbopt((nz-1)/2*2) = "(mPa)"

 do k=0, nz-1

 contour1 = cnshadeplot(wks,t,lat1,wss,cnlev,cnfil,lbopt(k))
 cnaddline(contour1,wks,t,lat1,uss(k,:,:),cnlevu,0,"D",(/0,0,1/),"T",(/0.5,1.5,0.5/),"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.01
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(k),yf(k),wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,"",xstr(k),ystr(k))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

; dr_txt2(wks,txtstr(k),xf(k)+0.6,yf(k)+he-2,22,6,0.,"CenterLeft",1)

end do

;  txtstrb = ":F22:"+lat0avg+":S:o:N:N:NF29:-:F22:"+lat0avg+":S:o:N:S"
;  dr_txt(wks,txtstrb,xf(3)+wi-8,yf(3)+2,22,6,0.)

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,6,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(3)-2,25,6,0.)


 frame(wks)

end



