load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 nmon = 54*12
 lat0 = 25.
 lat0avg = 5
 z0 = 18.59e3

 wks = cr_wks("fig/rwmfss-abs_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW/epf/uan"+expname(0)+".epfsa0_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW/epf/uan"+expname(0)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})
   z  = f1->z / 1.e3

   ny = dimsizes(lat1)
   nz = dimsizes(z)
   nt = 648
   nw = dimsizes(f1->wg) + 2

   rw_i = new((/647,nw,ny/),float)
   rw_i(:,:nw-3,:) = f1->f_z(:,:,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw_i(:,nw-2:,:) = f2->f_z(:,8:9,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(:,2,{z0},{-lat0:lat0}) + \
                      f2->f_z(:,3,{z0},{-lat0:lat0}) ) / 6370.e3 * 1.e3
   rw0 = new((/nt,nw,ny/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
   rw0(647,:,:) = rw0(646,:,:)

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,7:8,:)/)

   rw(:,:,0,:) = rw(:,:,0,:)*(-1.)

   u = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = rw(:,n,:,:)
   end do

   t = fspan(0,13,14)

   uss = new((/3,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )
   uss!2 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
; cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)*0.5
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 lcol = (/ (/1,3,23/), (/12,3,23/) /)

 cnlev2 = (/-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,-0.02,0, \
            0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)

;--------------------------------------- contour
 wname = (/ (/"Kelvin wave","E-MRGW","E-IGW"/), (/"Rossby wave","W-MRGW","W-IGW"/) /)
 wi = 30
 he = 24
 xf = 12 + (/0,1,2,3/)*38 + (/0,0,1,1/)*8
 yf = 10 + (/0,1,2,3,4/)*30
 ymax = (/0.55,1.10/)

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-F_z","Westward-F_z",lat0avg+"N:S:o:NF29:-:F22:"+lat0avg+":S:o:N:S"/)
 xstr = (/"","","month"/)
 ystr = (/"lat (deg)","","F_z (mPa)",""/)
 txtstr = (/ (/"Kelvin wave","E-MRGW","E-IGW"/), \
             (/"Rossby wave","W-MRGW","W-IGW"/) /)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(mPa)",""/)


do iv=0, 2

do idr=0, 1

 contour1 = cnshadeplot(wks,t,lat1,abs(uss(iv,idr,:,:)),cnlev,cnfil,lbopt(iv,idr))
 if (idr .eq. 0) then
   cnaddline(contour1,wks,t,lat1,uss(iv,idr,:,:),cnlev2,0,"D",0,"T",(/0.5,2.0,0.5/),"C",(/-1,1,1/))
 end if
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.18
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

  dr_txt2(wks,txtstr(idr,iv),xf(idr)+0.6,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

 plot = xylineplot(wks,t,dim_avg_n(uss(iv,idr,{-lat0avg*3:lat0avg*3},:),0),"D",0,"T",1.5,"C",lcol(idr,iv))
 set_size(plot,xf(2+idr),yf(0)+he/3,wi,he)
 axis_range(plot,0.5,12.5,0,ymax(0),False,False)
 axis_str(plot,mstr(0,idr),"",ystr(2+idr))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

 plot = xylineplot(wks,t,dim_avg_n(uss(iv,idr,{-lat0avg:lat0avg},:),0),"D",0,"T",1.5,"C",lcol(idr,iv))
 set_size(plot,xf(2+idr),yf(1)+he/3,wi,he)
 axis_range(plot,0.5,12.5,0,ymax(1),False,False)
 axis_str(plot,"",xstr(iv),ystr(2+idr))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

 lg = lgplot(wks,wname(idr,:),"D",0,"T",1.5,"C",lcol(idr,:),2)
 set_size_lg(lg,xf(2+idr),yf(2)+10,wi-6,20)
 draw(lg)

end do  ; idr

end do  ; iv

  dr_txt(wks,lat0avg*3+"N:S:o:NF29:-:F22:"+lat0avg*3+":S:o:N:S",xf(3)+wi-8,yf(0)+he/3+2,22,6,0.)
  dr_txt(wks,lat0avg+"N:S:o:NF29:-:F22:"+lat0avg+":S:o:N:S",xf(3)+wi-8,yf(1)+he/3+2,22,6,0.)
;  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end



