load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("test2",600,0)
 gsn_define_colormap(wks,"hotcold_18lev")
 icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
t1 = 24
t2 = 24
y0 = 2241

;--------------------------------------- read data

   f0  = addfile("../res/rho_ave_profile_p1.nc","r")
   rho = f0->rho / f0->rho({20000})

   f1o = addfile("../res/temp2/field/tp_p1_f10-e1.nc","r")

   f1 = addfile("../res/temp2/field/wp_p1_f10-e1.nc","r")
   time = f1->time({t1:t2:18})
   y    = f1->sout_north({y0})
   z    = f1->height/1.e3
   x    = f1->west_east

   var1t = f1o->t_up({t1:t2:18},:,{y0},:)
   var2t = f1o->t_dn({t1:t2:18},:,{y0},:)

   var1u = f1->u_up({t1:t2:18},:,{y0},:)
   var2u = f1->u_dn({t1:t2:18},:,{y0},:)
   var1w = f1->w_up({t1:t2:18},:,{y0},:)
   var2w = f1->w_dn({t1:t2:18},:,{y0},:)

 nt = dimsizes(time)

 do k=0, dimsizes(z)-1
   var1t(:,k,:) = var1t(:,k,:)*sqrt(rho(k))
   var2t(:,k,:) = var2t(:,k,:)*sqrt(rho(k))
   var1u(:,k,:) = var1u(:,k,:)*sqrt(rho(k))
   var2u(:,k,:) = var2u(:,k,:)*sqrt(rho(k))
   var1w(:,k,:) = var1w(:,k,:)*sqrt(rho(k))
   var2w(:,k,:) = var2w(:,k,:)*sqrt(rho(k))
 end do

 var0t2 = (var1t+var2t)^2

 var0u = var1u + var2u
 var0w = var1w + var2w
 var0uw = (var1u + var2u)*(var1w + var2w)
 var0ut = var0u*(var1t+var2t)
 utmn = dim_avg(var0ut)

;--------------------------------------- line attribute
 cnlev0 = (/0.5,1,2,3,4,6,8,10/)
 ncnl0 = dimsizes(cnlev0)
 cnfil0 = icol_center + ispan(0,ncnl0,1)*1
 cnfil0 = cnfil0 + 1
 cnfil0(0) = 0

 cnlev = (/-3,-2,-1,-0.5,-0.1,0.1,0.5,1,2,3/)
 ncnl = dimsizes(cnlev)
 cnfil = icol_center + ispan(-ncnl/2,ncnl/2,1)*1
 cnfil(:ncnl/2-1) = cnfil(:ncnl/2-1) - 3
 cnfil(ncnl/2+1:) = cnfil(ncnl/2+1:) + 1

 cnlevu = (/-5,-3,-2,-1,1,2,3,5/)
 cnlevw = (/-10,-6,-4,-2,2,4,6,10/)
 ncnl1 = dimsizes(cnlevu)
 cnfil1 = icol_center + (/-8,-6,-5,-4,0,2,3,4,6/)
;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 10 + (/0,1,2/)*58 ;+ (/0,6,6/)
 yf = 10 + (/0,1/)*68

 mstr = (/"T:S:2:N:","uw","uT","w"/)

do n=0, nt-1

 contour = cnshadeplot(wks,x,z,var0t2(n,:,:)*1.e2, \
           cnlev0,cnfil0,"[10:S:-2:N: K:S:2:N:]")

 set_size(contour,xf(0),yf(0),wi,he)
 axis_range(contour,-999,999,20,999,False,False)
 axis_str(contour,mstr(0),"x [km]","z [km]")

 cnfinalize(contour,wks)

 contour0 = cnshadeplot(wks,x,z,var0uw(n,:,:)*1.e3, \
           cnlev,cnfil,"[10:S:-3:N: m:S:2:N: s:S:-2:N:]")

 set_size(contour0,xf(1),yf(0),wi,he)
 axis_range(contour0,-999,999,20,999,False,False)
 axis_str(contour0,mstr(1),"x [km]","z [km]")

 cnfinalize(contour0,wks)

 contour1 = cnshadeplot(wks,x,z,var0ut(n,:,:)*1.e2, \
           cnlevu,cnfil1,"[10:S:-2:N: m s:S:-1:N:]")

 set_size(contour1,xf(0),yf(1),wi,he)
 axis_range(contour1,-999,999,20,999,False,False)
 axis_str(contour1,mstr(2),"x [km]","z [km]")

 cnfinalize(contour1,wks)

; contour2 = cnshadeplot(wks,x,z,var0w(n,:,:)*1.e3, \
;           cnlevw,cnfil1,"[10:S:-3:N: m s:S:-1:N:]")
;
;; cnaddline(contour2,wks,x,z,var0uw(n,:,:)*1.e2,(/0.1/),0, \
;;           "D",0,"T",0.6,"C",1)
;
; set_size(contour2,xf(1),yf(1),wi,he)
; axis_range(contour2,-999,999,20,999,False,False)
; axis_str(contour2,mstr(3),"x [km]","z [km]")
;
; cnfinalize(contour2,wks)


 draw(contour)
 draw(contour0)
 draw(contour1)
; draw(contour2)

 plot = xylineplot(wks,utmn,z,"D",1,"T",1.2,"C",1)

 set_size(plot,xf(1),yf(1),wi,he)
 axis_range(plot,-999,999,20,45,False,False)
 axis_str(plot,"","","")

 draw(plot)

dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
dr_txt(wks,"(b)",xf(1)-2,yf(0)-4,25,8,0.)
dr_txt(wks,"(c)",xf(0)-2,yf(1)-4,25,8,0.)
dr_txt(wks,"(d)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end do  ; n

end

