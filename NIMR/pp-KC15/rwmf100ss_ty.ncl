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
 z0 = 17.5e3

 wks = cr_wks("fig/rwmf100ss_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/epf/uanuj.epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   lat1 = f1->lat({-lat0:lat0})
   z  = f1->z / 1.e3

   ny = dimsizes(lat1)
   nz = dimsizes(z)
   nt = 648
   nw = dimsizes(f1->wg)

   rw_i = f1->f_z(:,:,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw0 = new((/nt,nw,ny/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
   rw0(647,:,:) = rw0(646,:,:)

   rw0a = rw0
   rw0a(1:646,:,:) = 0.5*(rw0(0:645,:,:) + rw0(2:647,:,:))

   rw0 = where(abs(rw0) .ge. 1.e10, rw0a, rw0)

   rw = (/ rw0(:,4:5,:)+rw0(:,6:7,:), rw0(:,8:9,:), rw0(:,0:1,:) /)

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
 cnlev = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.5/)
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38 + (/0,0,1/)*8
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-F_z","Westward-F_z",lat0avg+"N:S:o:NF29:-:F22:"+lat0avg+":S:o:N:S"/)
 xstr = (/"","","month"/)
 ystr = (/"lat (deg)","","F_z (hPa)"/)
 txtstr = (/"Large-scale waves","Inertia-GWs","Total resolved waves"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(mPa)",""/)


do iv=0, 2

do idr=0, 1

 contour1 = cnshadeplot(wks,t,lat1,uss(iv,idr,:,:),cnlev,cnfil,lbopt(iv,idr))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.18
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.4
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

  dr_txt2(wks,txtstr(iv),xf(0)+0.6,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; idr

 plot = xylineplot(wks,t,dim_avg_n(uss(iv,:,{-lat0avg:lat0avg},:),1),"D",0,"T",1.2,"C",(/3,23/))
 set_size(plot,xf(2),yf(iv),wi,he)
 axis_range(plot,0.5,12.5,0,0.9,False,False)
 axis_str(plot,mstr(iv,2),xstr(iv),ystr(2))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

end do  ; iv

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,8,0.)
;  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end



