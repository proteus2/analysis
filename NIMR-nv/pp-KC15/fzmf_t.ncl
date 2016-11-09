load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 nmon = 54*12
 lat0 = 20.
 lat0avg = 5
 z0 = 26.6e3 ;26.89e3 ;17.41e3

 wks = cr_wks("fig/fzmf_t-z"+numeric2int(z0/1.e3,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW-t/epf/uan"+expname(0)+".epfsa0_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf/uan"+expname(0)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})
print(f1->z({z0}))

   ny = dimsizes(lat1)
   nt = 648
   nw = dimsizes(f1->wg) + 2

   rw_i = new((/647,nw,ny/),float)
   rw_i(:,:nw-3,:) = f1->f_z(:,:,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw_i(:,nw-2:,:) = f2->f_z(:,8:9,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) ;+ \
;     ( f2->f_z(:,2,{z0},{-lat0:lat0}) + \
;       f2->f_z(:,3,{z0},{-lat0:lat0}) ) / 6370.e3 * 1.e3
;     tide
   rw0 = new((/nt,nw,ny/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
;   rw0(647,:,:) = rw0(646,:,:)
   rw0(647,:,:) = rw0(647-12,:,:)

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,7:8,:)/)

;   rw(:,:,0,:) = rw(:,:,0,:)*(-1.)


   fnam1 = "../dat/L60CGW-t/mmxm.uanuja.pd.nc"

   f1 = addfile(fnam1,"r")
   mfe1 = f1->F_e_gwconv(:,{z0},{-lat0:lat0})
   mfw1 = f1->F_w_gwconv(:,{z0},{-lat0:lat0})*(-1.)
;   lat1 = f1->latitude({-lat0:lat0})
print(f1->hybrid_ht_1({z0}))

   fnam2 = "../dat/L60CGW-t/mmgwusspz_xm.uanuja.pg.nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->Fp_e_gwussp(:,{z0},{-lat0:lat0})
   mfw2 = f2->Fp_w_gwussp(:,{z0},{-lat0:lat0})
;   mfe2 = mfe2 + mfe1
;   mfw2 = mfw2 + mfw1

   u = new((/2,3,2,ny,nt/),float)
   do n=0, nt-1
     u(0,:,:,:,n) = rw(:,n,:,:)
     u(1,1:2,:,:,n) = (/ (/mfe1(n,:),mfw1(n,:)/), \
                       (/mfe2(n,:),mfw2(n,:)/) /)*1.e3
   end do
   u(1,0,:,:,:) = u(1,1,:,:,:) + u(1,2,:,:,:)


   t = fspan(0,13,14)

   u_cli = dim_avg(u)
   u_cli!3 = "lat"
   u_cli&lat = lat1

   uss = new((/2,3,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,:,n) = dim_avg(u(:,:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,:,0 ) = uss(:,:,:,:,12)
   uss(:,:,:,:,13) = uss(:,:,:,:,1 )
   uss!3 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0

;--------------------------------------- cn level
; lcol = (/ (/1,3,23/), (/12,3,23/) /)
; ldsh = (/0,0,0/)
 ltck = (/1.2,1.2,1.2,1.8/)
 ldsh = (/0,16,14,0/)

;--------------------------------------- contour
 wname = (/ (/"Kelvin wave","E-MRGW","E-IGW","E-GW"/), \
            (/"Rossby wave","W-MRGW","W-IGW","W-GW"/) /)
 wi = 30
 he = 24
 xf = 12 + (/0,1,2,3/)*38
 yf = 10 + (/0,1,2,3/)*30

 mstr = (/ (/"Eastward waves","Westward waves"/), (/"",""/) /)
 xstr = (/"","latitude (deg)"/)
 ystr = (/ (/"F:B:z:N: (mPa)","M (mPa)"/), (/"",""/) /)


do idr=0, 1

 plot_y = xylineplot(wks,lat1,u_cli(0,:,idr,:),"D",ldsh(:2),"T",ltck(:2),"C",1)
 xyaddline(plot_y,lat1,0.,"D",0,"T",0.8,"C",colind)
 set_size(plot_y,xf(idr),yf(0),wi,he)
 if (idr .eq. 0) then
   axis_range(plot_y,-20,20,0.,-0.15,False,False)
 else
   axis_range(plot_y,-20,20,0.,0.15,False,False)
 end if
 axis_str(plot_y,mstr(0,idr),xstr(0),ystr(idr,0))
 draw(plot_y)

 plot_y = xylineplot(wks,lat1,u_cli(1,0,idr,:),"D",ldsh(3),"T",ltck(3),"C",1)
 set_size(plot_y,xf(idr),yf(1),wi,he)
 axis_range(plot_y,-20,20,0.,0.9,False,False)
 axis_str(plot_y,mstr(1,idr),xstr(1),ystr(idr,1))
 draw(plot_y)

 lg = lgplot(wks,wname(idr,:),"D",ldsh,"T",ltck,"C",1,2)
 set_size_lg(lg,xf(idr),yf(1)+32,wi-8,15)
 draw(lg)

end do

  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,6,0.)
  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,6,0.)


 frame(wks)

end



