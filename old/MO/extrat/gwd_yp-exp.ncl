load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 expn = (/"i","0"/)

 month = (/1,7/)

;--------------------------------------- parameters
 mstr = (/"CGWD + NGWD","NGWD (CTL)","CGWD","NGWD"/)

;--------------------------------------- read data
 dir_dat   = "$PORTAL/MO/dat"
 file_dat1 = dir_dat+"/ng_"+expn(1)+".nc"
 file_dat2 = dir_dat+"/ng_"+expn(0)+".nc"
 file_dat3 = dir_dat+"/cg_"+expn(0)+".nc"

 f1 = addfile(file_dat1,"r")
 f2 = addfile(file_dat2,"r")
 f3 = addfile(file_dat3,"r")

 lat = f1->latitude
 z   = f1->hybrid_ht/1.e3

 nmon = dimsizes(f3->t(4:))

 u1i = f1->field424(4:,:,:,0)*86400.
 u2i = f2->field424(4:,:,:,0)*86400.
 u3i = f3->unspecified(4:,:,:,0)*86400.

 ny = dimsizes(lat)
 nz = dimsizes(z)

 nyr = new(12,integer)
 um0 = new((/12,nz,ny/),float,1.e32)
 um1 = new((/12,nz,ny/),float,1.e32)
 um2 = new((/12,nz,ny/),float,1.e32)
 do n=0, 11
   nyr(n) = dimsizes(ispan(n,nmon-1,12))
   if (nyr(n) .gt. 1) then
     um0(n,:,:) = dim_avg_n( u1i(ispan(n,nmon-1,12),:,:),0 )
     um1(n,:,:) = dim_avg_n( u2i(ispan(n,nmon-1,12),:,:),0 )
     um2(n,:,:) = dim_avg_n( u3i(ispan(n,nmon-1,12),:,:),0 )
   else
     um0(n,:,:) = u1i(n,:,:)
     um1(n,:,:) = u2i(n,:,:)
     um2(n,:,:) = u3i(n,:,:)
   end if
 end do

;--------------------------------------- process 0
 do im=0, dimsizes(month)-1

 wks = cr_wks("gwd_yp-"+month(im)+"-"+expn(0),600,0)
 gsn_define_colormap(wks,"posneg_2")

 print(nyr(month(im)-1)+"")

 uf1 = um0(month(im)-1,:,:)
 uf2 = um2(month(im)-1,:,:)
 uf3 = um1(month(im)-1,:,:)
 uf0 = uf2 + uf3

;--------------------------------------- line attribute
 cnlev = (/-5,-3,-2,-1,0,1,2,3,5/)*0.1
 cnfil = (/6,8,9,10,0,0,13,14,15,17/)

 cnlev2 = (/-50,-20,-10,-5,-2,-1,1,2,5,10,20,50/)*1.  ;fspan(-20,20,41)

;--------------------------------------- plot
 wi = 58
 he = 32
 xf = 20 + (/0,1/)*72
 yf = 10 + (/0,1/)*56

 contour1 = cnshadeplot(wks,lat,z,uf0,cnlev,cnfil,1)
 cnaddline(contour1,wks,lat,z,uf0,cnlev2,0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour1,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,-90,90,10,64,False,False)
 axis_str(contour1,mstr(0),"lat","model height (km)")

 cnfinalize(contour1,wks)

 draw(contour1)


 contour2 = cnshadeplot(wks,lat,z,uf1,cnlev,cnfil,1)
 cnaddline(contour2,wks,lat,z,uf1,cnlev2,0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour2,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour2,xf(1),yf(0),wi,he)
 axis_range(contour2,-90,90,10,64,False,False)
 axis_str(contour2,mstr(1),"lat","model height (km)")

 cnfinalize(contour2,wks)

 draw(contour2)


 contour3 = cnshadeplot(wks,lat,z,uf2,cnlev,cnfil,1)
 cnaddline(contour3,wks,lat,z,uf2,cnlev2,0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour3,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour3,xf(0),yf(1),wi,he)
 axis_range(contour3,-90,90,10,64,False,False)
 axis_str(contour3,mstr(2),"lat","model height (km)")

 cnfinalize(contour3,wks)

 draw(contour3)


 contour4 = cnshadeplot(wks,lat,z,uf3,cnlev,cnfil,1)
 cnaddline(contour4,wks,lat,z,uf3,cnlev2,0, \
           "D",(/0,0,0/),"T",(/0.8,0.5,0.8/),"C",(/1,0,1/))
; cnline_att(contour4,"D",(/0,0,1/),"T",(/0.8,0.5,0.8/),"C",1)

 set_size(contour4,xf(1),yf(1),wi,he)
 axis_range(contour4,-90,90,10,64,False,False)
 axis_str(contour4,mstr(3),"lat","model height (km)")

 cnfinalize(contour4,wks)

 draw(contour4)


 frame(wks)


 end do

end

