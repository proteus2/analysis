load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
; expname = "era-int_f"
; dirname = "ERA-I-nr"
; uname   = "u"
 expname = "jra55c_f"
 dirname = "JRA55C-nr"
 uname   = "U_GRD_GDS4_HYBL"

 lat0 = 25.
 p0    = 100.  ; 54, 66, 80, 96
if (expname .eq. "era-int_f") then
 p0u   = 23   ; 22, 23, 24, 25
else
 p0u   = 38   ; 39, 38, 37, 36
end if

 wks = cr_wks("fig/epf-wig_ty-p"+numeric2int(p0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
   fnam0 = "$DATD/clim/"+dirname(0)+"/"+expname(0)+".u_ypt.nc"
   fnam1 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")

   f0 = addfile(fnam0,"r")
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nt = 1 + 407 + 36
if (expname .eq. "jra55c_f") then
   nt = 1 + 407 - 12
end if
   nw = 8

   u = f0->$uname(0)$(:nt-1,{p0u},{-lat0:lat0})

   rw_i = new((/nt-1,nw,ny/),float)
   rw_i(:,:,:) = f1->f_z(:nt-2,:nw-1,{p0},{-lat0:lat0})
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     f2->f_z(:nt-2,8:9,{p0},{-lat0:lat0})
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(:nt-2,2,{p0},{-lat0:lat0}) + \
                      f2->f_z(:nt-2,3,{p0},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(:nt-2,13,{p0},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(:nt-2,16:17,{p0},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,nw,ny/),float,1.e20)
   rw0!2 = "lat"
   rw0&lat = lat1
   rw0(1:,:,:) = rw_i(:,:,:)
;   rw0(0,:,:) = rw0(1,:,:)
   rw0(0,:,:) = 1.e20

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   fz = new((/3,2,ny,nt/),float,1.e20)
   fz!2 = "lat"
   fz&lat = lat1
   do n=0, nt-1
     fz(:,:,:,n) = (/ rw(:,n,:,:) /)
   end do

   t = 1979 + ( ispan(0,nt-1,1)+0.5 )/12.

   fzss = new((/3,2,ny,14/),float,1.e20)
   fzss!2 = "lat"
   fzss&lat = lat1
   do n=1, 12
     fzss(:,:,:,n) = dim_avg(fz(:,:,:,ispan(n-1,nt-1,12)))
   end do
   fzss(:,:,:,0 ) = fzss(:,:,:,12)
   fzss(:,:,:,13) = fzss(:,:,:,1 )

;--------------------------------------- process 0
 fzds = fz
 fzssp = fz
 do n=1, 12
   nn = ispan(n-1,nt-1,12)
   fzds(:,:,:,nn) = fz(:,:,:,nn) -  \
       conform_dims((/3,2,ny,dimsizes(nn)/),fzss(:,:,:,n),(/0,1,2/))
   fzssp(:,:,:,nn) = conform_dims((/3,2,ny,dimsizes(nn)/),fzss(:,:,:,n),(/0,1,2/))
   delete(nn)
 end do

;--------------------------------------- cn level
 cnlev = (/0.5,1,2,3,4,5,6,7,8,9,10,11,12/)*0.1
 cnfil = (/5,6,7,8,10,12,14,15,16,17,18,19,20,21/)

;--------------------------------------- contour
 wi = 50
 he = 40
 xf = 15
 yf = 10 + ispan(0,3,1)*(he+10)

 ystr = "latitude [deg]"
 xstr = "month"
 lbtxt = "F:B:z:N: [mPa]"

ne = 1

fzssr = -fzss(0,0,:,:)/conform_dims(dimsizes(fzss(0,0,:,:)),dim_max_n(-fzss(0,0,:,:),0),1)

do ie=0, ne-1

 contour1 = cnshadeplot(wks,ispan(0,13,1),lat1,fzss(2,1,:,:),cnlev,cnfil,lbtxt)
; cnaddline(contour1,wks,ispan(0,13,1),lat1,fzssr,exp(-1.),0,"D",0,"T",1.0,"C",0)
 cnaddline(contour1,wks,ispan(0,13,1),lat1,fzss(2,1,:,:),cnlev,0,"D",0,"T",0.5,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
;  "pmLabelBarParallelPosF"   : 0.50-0.40
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.60
  "pmLabelBarHeightF"        : 0.06
  "lbLabelFontHeightF" : 0.012  ; 0.011
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.015
 end setvalues
 set_size(contour1,xf,yf(ie),wi,he)
 axis_range(contour1,0.5,12.5,-999,999,False,False)
 axis_str(contour1,"",xstr,ystr)
 setvalues contour1
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : ispan(1,12,1)
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

end do

 frame(wks)

end

