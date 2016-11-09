load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.4

begin

 wks = cr_wks("t2e_z",300,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 colind2 = NhlNewColor(wks,0.8,0.8,0.8)

;--------------------------------------- parameters
t1 = 24
t2 = 24
y0 = 2241

 nf = 1
 nl = 3

;--------------------------------------- read data

   f0  = addfile("/data4/atmosdata/T_SPECTRUM/p1/rho_ave_profile_p1.nc","r")
   rho = f0->rho / f0->rho({20000})

   f1 = addfile("../res/temp2/field/tp_p1_f10-e1.nc","r")
   time = f1->time({t1:t2:18})
   y    = f1->sout_north({y0})
   z    = f1->height({20000:})/1.e3
   x    = f1->west_east

   var1 = f1->t_up({t1:t2:18},{20000:},{y0},:)
   var2 = f1->t_dn({t1:t2:18},{20000:},{y0},:)

 nt = dimsizes(time)
 nz = dimsizes(z)

 do k=0, dimsizes(z)-1
   var1(:,k,:) = var1(:,k,:)*sqrt(rho(k))
   var2(:,k,:) = var2(:,k,:)*sqrt(rho(k))
 end do

 var0 = var1 + var2

;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat = 1.e32
 var_dat(0,0,:) = dim_avg(var1^2 - var2^2)
 var_dat(0,1,:) = dim_avg(var0^2)
 kavg = 10
 do k=kavg, nz-1-kavg
   var_dat(0,2,k) = avg(var_dat(0,1,k-kavg:k+kavg))
 end do
 var_dat(0,:,:kavg-1) = 1.e32
 var_dat(0,:,nz-kavg:) = 1.e32
;;;

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
; labels = (/"F:B:W:N:","F:B:T:N:",/)
 lcol = (/1,2,4,2/)
 ltck = (/1.3,1.3,1.3,1.3/)
 ldsh = (/0,0,0,16/)

;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 12 + (/0,1,2/)*(wi+10)
 yf = 10 + (/0,0,0/)*(he+10)

 ystr = (/"z [km]","","","z [km]","",""/)

do ifig=0, nf-1

 plot = xylineplot(wks,var(ifig,:,:),z,"D",ldsh,"T",ltck,"C",lcol)

 xyaddline(plot,0,z,"D",0,"T",0.6,"C",colind2)

 set_size(plot,xf(ifig),yf(ifig),wi,he)
 axis_range(plot,0,0.32,-999,45,False,False)
 axis_str(plot,"","T:S:2:N: [K:S:2:N:]",ystr(ifig))
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)

; lg1 = lgplot(wks,labels,"D",ldsh, \
;                  "T",ltck,"C",lcol,1)
;
; set_size_lg(lg1,xf(ifig)+wi-10,yf,10,8)
;
; if (ifig .eq. 0) then
;   draw(lg1)
; end if

; dr_txt(wks,"date",xf+wi/2,yf+he+10,22,6,0.)

end do

;dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

