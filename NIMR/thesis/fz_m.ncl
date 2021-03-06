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
 lat0avg = 5
 z0    = 18.59  ; 18.59  ; 17.41
 ymax0 = 1.0    ; 0.5    ; 0.8

 wks = cr_wks("fig/fz_m-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"
   fnam3 = "$DATD/L60CGW-t/mmxm.uan"+expname(0)+"a.pd.nc"
   fnam4 = "$DATD/L60CGW-t/mmxm.uan"+expname(0)+"a.pg.nc"
print(fnam1+"")
print(fnam2+"")
print(fnam3+"")
print(fnam4+"")

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")
   f4 = addfile(fnam4,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nt = 647 + 1
   nwr = 8
   nwp = 3

   rw_i = new((/nt-1,nwr,ny/),float)
   rw_i(:,:,:) = f1->f_z(:,:nwr-1,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nwr-2:,:) = rw_i(:,nwr-2:,:) + \
                     f2->f_z(:,8:9,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nwr-1,:) = rw_i(:,nwr-1,:) + \
                    ( f2->f_z(:,2,{z0*1.e3},{-lat0:lat0}) + \
                      f2->f_z(:,3,{z0*1.e3},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(:,13,{z0*1.e3},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(:,16:17,{z0*1.e3},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,nwr,ny/),float)
   rw0(:nt-2,:,:) = rw_i(:,:,:)
   rw0(nt-1,:,:) = rw0(nt-2,:,:)

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   pw = new((/nwp,2,nt,ny/),float)
   pw(0,0,:,:) = f3->F_e_gwconv (:,{z0*1.e3},{-lat0:lat0})
   pw(0,1,:,:) = f3->F_w_gwconv (:,{z0*1.e3},{-lat0:lat0})
   pw(1,0,:,:) = f4->Fp_e_gwussp(:,{z0*1.e3},{-lat0:lat0})
   pw(1,1,:,:) = f4->Fp_w_gwussp(:,{z0*1.e3},{-lat0:lat0})*(-1.)
   pw(2,0,:,:) = pw(0,0,:,:) + pw(1,0,:,:)
   pw(2,1,:,:) = pw(0,1,:,:) + pw(1,1,:,:)

   pw(:,:,:,:) = pw(:,:,:,:)*(-1.)*1.e3

   u = new((/5,2,ny,nt/),float)
   do j=0, ny-1
   do n=0, nt-1
     u(0:2,:,j,n) = (/ rw(0:2,n,:,j) /)
;     u(3:5,:,:,n) = (/ pw(0:2,:,n,j) /)
     u(3:4,:,j,n) = (/ pw((/0,2/),:,n,j) /)*cos(lat1(j)*3.141592/180.)
   end do
   end do

   t = fspan(0,13,14)

   uss = new((/5,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )
   uss!2 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10/)*0.1
 cnfil = (/5,6,7,8,10,12,14,15,16,17,18,19,20,21/)

 ldsh  = (/0,  0,0,14,0/)
 ldsh2 = (/0,0,0,0,14,0/)
 lcol  = (/ (/1,7,3,2,2/), \
            (/4,7,3,2,2/) /)
 lcol2 =  (/1,4,7,3,2,2/)

;--------------------------------------- contour
 wname = (/"Kelvin wave","Rossby wave","MRG wave","IG wave", \
           "CGW","CGW+BGW"/)
 wi = 30
 he = 22
 xf = 12 + (/0,1/)*38
 yf = 10

 mstr = (/"Eastward waves","Westward waves"/)
 ystr = (/":F10:F:B:z:NF22: (mPa)",""/)
 ymax = ymax0*(/-1,1/)

do idr=0, 1

 plot = xylineplot(wks,t,dim_avg_n(uss(:,idr,{-lat0avg:lat0avg},:),1),"D",ldsh,"T",1.5,"C",lcol(idr,:))
 set_size(plot,xf(idr),yf,wi,he)
 axis_range(plot,0.5,12.5,-ymax(idr)*0.01,ymax(idr),False,False)
 axis_str(plot,mstr(idr),"month",ystr(idr))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

end do  ; idr

 lg = lgplot(wks,wname(:2),"D",ldsh2(:2),"T",1.5,"C",lcol2(:2),2)
 set_size_lg(lg,xf(1)-20,yf+he+8,wi*0.6,12)
 draw(lg)
 lg2 = lgplot(wks,wname(3:),"D",ldsh2(3:),"T",1.5,"C",lcol2(3:),2)
 set_size_lg(lg2,xf(1)+5,yf+he+8,wi*0.6,12)
 draw(lg2)


 frame(wks)

end



