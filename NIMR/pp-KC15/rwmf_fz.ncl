load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 lat0 = 5.

 wks = cr_wks("fig/rwmf_fz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/epf/uanuj.epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   lat = f1->lat({-lat0:lat0})
   z  = f1->z / 1.e3

   ny = dimsizes(lat)
   nz = dimsizes(z)
   nt = 648/3
   nw = dimsizes(f1->wg)

   rw_i = dim_avg( f1->f_z(:,:,:,{-lat0:lat0}) ) / 6370.e3 * 1.e3
   rw0 = new((/nt*3,nw,nz/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
   rw0(647,:,:) = rw0(646,:,:)

   rw0a = rw0
   rw0a(1:646,:,:) = 0.5*(rw0(0:645,:,:) + rw0(2:647,:,:))

   rw0 = where(abs(rw0) .ge. 1.e10, rw0a, rw0)

   rw0t = new((/nw,3,nz,nt/),float)
   do iw=0, nw-1
   do k=0, nz-1
     rw0t(iw,0,k,:) = rw0(0:215,iw,k)
     rw0t(iw,1,k,:) = rw0(216:431,iw,k)
     rw0t(iw,2,k,:) = rw0(432:,iw,k)
   end do
   end do

   rw = (/ rw0t(4:5,:,:,:)+rw0t(6:7,:,:,:), rw0t(8:9,:,:,:), rw0t(0:1,:,:,:) /)

   rw(:,0,:,:,:) = rw(:,0,:,:,:)*(-1.)

;--------------------------------------- process 0
 u = dtrend(rw ,False)

 win = sin(3.141592*fspan(0,1,nt))        ; cosine window
win = 1.

 mswin = avg((win)^2)

 do n=0, nt-1
   u(:,:,:,:,n) = u(:,:,:,:,n)*win(n)
 end do

   nf  = nt/2
   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
;print(fr(0)+"  "+fr3(0))

   fftu  = ezfftf(u )

   pws0 = ( fftu (0,:,:,:,:,:)^2 + fftu (1,:,:,:,:,:)^2 )*0.5*nt /12. / mswin

   pws = dim_avg_n(pws0, 2)

 pws = sqrt(pws*fr(0))/fr(0)

;--------------------------------------- cn level
 cnlev = (/0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6/)
 cnfil = (/0,6,4,2,26,24,22,21,19,17/)

;--------------------------------------- contour
 wi = 60
 he = 24
 xf = 12 + (/0,1,2/)*68 + (/0,0,0/)*6.6
 yf = 10 + (/0,1,2,3,4/)*30

 varn = (/"Large-scale waves","Inertia-GWs","Total resolved waves"/)
 mstr = (/ (/"Eastward-F_z","Westward-F_z"/), (/"",""/), (/"",""/) /)
 xstr = (/"","","frequency (cycle yr:S:-1:N:)"/)
 ystr = (/"height (km)","",""/) ;"log-p height (km)"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(mPa / cycle yr:S:-1:N:)",""/) ;(/"(m s:S:-1:N:)","",""/)


do ie=0, 2
do iv=0, 1

 contour1 = cnshadeplot(wks,fr,z,pws(ie,iv,:,:),cnlev,cnfil,lbopt(ie,iv))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.2 ;1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
  "pmLabelBarHeightF"        : 0.2
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.020
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.25
  "lbTitleFontHeightF" : 0.024
 end setvalues

 set_size(contour1,xf(iv),yf(ie),wi,he)
 axis_range(contour1,-999,2.5,17.5,40,False,False)
 axis_str(contour1,mstr(ie,iv),xstr(ie),ystr(iv))
 setvalues contour1
  "tmXBMaxTicks" : 5
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt2(wks,varn(ie),xf(iv)+wi-1,yf(ie)+2,22,6,0.,"CenterRight",1)

end do
end do

 frame(wks)

end



