load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.4

begin

 wks = cr_wks("t_xz",600,0)
 gsn_define_colormap(wks,"hotcold_18lev")
 icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
t1 = 24
t2 = 24
y0 = 2241

;--------------------------------------- read data

   f0  = addfile("/data4/atmosdata/T_SPECTRUM/p1/rho_ave_profile_p1.nc","r")
   rho = f0->rho / f0->rho({20000})

   f1 = addfile("../res/temp2/field/tp_p1_f10-e1.nc","r")
   time = f1->time({t1:t2:18})
   y    = f1->sout_north({y0})
   z    = f1->height/1.e3
   x    = f1->west_east

   var1 = f1->t_up({t1:t2:18},:,{y0},:)
   var2 = f1->t_dn({t1:t2:18},:,{y0},:)

 nt = dimsizes(time)

 do k=0, dimsizes(z)-1
   var1(:,k,:) = var1(:,k,:)*sqrt(rho(k))
   var2(:,k,:) = var2(:,k,:)*sqrt(rho(k))
 end do

 var0 = var1 + var2

;--------------------------------------- line attribute
 cnlev0 = (/0.5,1,2,3,4,6,8,10/)
 ncnl0 = dimsizes(cnlev0)
 cnfil0 = icol_center + ispan(0,ncnl0,1)*1
 cnfil0 = cnfil0 + 1
 cnfil0(0) = 0

 cnlev = (/-2,-1.5,-1,-0.5,0,0.5,1,1.5,2/)
 ncnl = dimsizes(cnlev)
 cnfil = icol_center + ispan(-ncnl/2,ncnl/2+1,1)*1
 cnfil(:ncnl/2) = cnfil(:ncnl/2) - 2
 cnfil(ncnl/2+1:) = cnfil(ncnl/2+1:) + 1
;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 12 + (/0,1,2/)*53 ;+ (/0,6,6/)
 yf = 10

 mstr = (/"T","T:B:up:N","T:B:down:N"/)

do n=0, nt-1

 contour0 = cnshadeplot(wks,x,z,var0(n,:,:)*1.e1, \
           cnlev,cnfil,0)

 set_size(contour0,xf(0),yf,wi,he)
 axis_range(contour0,-999,999,20,999,False,False)
 axis_str(contour0,mstr(0),"x [km]","z [km]")

 cnfinalize(contour0,wks)

 contour1 = cnshadeplot(wks,x,z,var1(n,:,:)*1.e1, \
           cnlev,cnfil,"[10:S:-1:N: K]")

 set_size(contour1,xf(1),yf,wi,he)
 axis_range(contour1,-999,999,20,999,False,False)
 axis_str(contour1,mstr(1),"x [km]","")

 cnfinalize(contour1,wks)

 contour2 = cnshadeplot(wks,x,z,var2(n,:,:)*1.e1, \
           cnlev,cnfil,0)

; cnaddline(contour2,wks,x,z,var2(n,:,:)*1.e1,(/-sqrt(0.1),sqrt(0.1)/),0, \
;           "D",0,"T",0.6,"C",1)

 set_size(contour2,xf(2),yf,wi,he)
 axis_range(contour2,-999,999,20,999,False,False)
 axis_str(contour2,mstr(2),"x [km]","")

 cnfinalize(contour2,wks)


 draw(contour0)
 draw(contour1)
 draw(contour2)

;dr_txt(wks,"(a)",xf(0)-2,yf-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf-4,25,8,0.)
;dr_txt(wks,"(c)",xf(2)-2,yf-4,25,8,0.)


 frame(wks)

end do  ; n

end

