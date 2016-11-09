load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 wks = cr_wks("./numdmp",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
; f0 = addfile("../dat/arch/uj/pe/uanuja.pe_195801.nc","r")
; f1 = addfile("../dat/arch/uj/pa/uanuja.pa_1958010100.nc","r")
; f2 = addfile("../dat/arch/uj/pa/uanuja.pa_1958020100.nc","r")
; f3 = addfile("../dat/tem/1958/uanuj.tem_yzt.1958.01.nc","r")
; f4 = addfile("../dat/arch/uj/pj/uanuja.pj_195801.nc","r")
 f0 = addfile("/data7/kyh/arch/up/pe/uanupa.pe_194401.nc","r")
 f1 = addfile("/data7/kyh/arch/up/pa/uanupa.pa_1944010100.nc","r")
 f2 = addfile("/data7/kyh/arch/up/pa/uanupa.pa_1944020100.nc","r")
 f3 = addfile("../dat/tem_up/1944/uanup.tem_yzt.1944.01.nc","r")
 f4 = addfile("/data7/kyh/arch/up/pj/uanupa.pj_194401.nc","r")

 du = (dim_avg(f2->u(23,:,:,:)) - dim_avg(f1->u(23,:,:,:)))/30.

 y = f0->latitude
 z = f0->hybrid_ht_1/1.e3
 ny = dimsizes(y)
 nz = dimsizes(z)

 var = new((/6,nz,ny/),float)

 var(0,:,:) = dim_avg(f0->du_gw    (0,:,:,:))*72.
 var(1,:,:) = dim_avg(f0->du_adv   (0,:,:,:))*72.
 var(2,:,:) = dim_avg(f0->du_bl    (0,:,:,:))*72.
 var(3,:,:) = dim_avg(f0->du_solver(0,:,:,:))*72.
 var(4,:,:) = dim_avg(f0->dudt_conv_p(0,:,:,:))*86400.

 var(5,:,:) = dim_sum_n(var(0:4,:,:),0)

 utend = dim_avg_n( f3->u_tend(0,:,:), 0 )

 var2 = new((/11,nz,ny/),float)
 var2(0,:,:) = dim_avg_n( f3->cor    (:,:,:), 0 )
 var2(1,:,:) = dim_avg_n( f3->uadv_y (:,:,:), 0 )
 var2(2,:,:) = dim_avg_n( f3->uadv_z (:,:,:), 0 )
 var2(3,:,:) = dim_avg_n( f3->epd    (:,:,:), 0 )
 var2(4,:,:) = 0;f3->ru_eddy(0,:,:)
 var2(5,:,:) = dim_avg_n( f3->u_force(:,:,:), 0 ) + var(0,:,:) + var(2,:,:) + var(4,:,:)
 var2(6,:,:) = dim_avg_n( f3->u_tend (:,:,:), 0 )

 var2(7,:,:)  = dim_avg_n( f3->epd_y(:,:,:), 0 )
 var2(8,:,:)  = dim_avg_n( f3->epd_z(:,:,:), 0 )
 var2(9,:,:)  = dim_avg_n( f3->v_res(:,:,:), 0 )
 var2(10,:,:) = dim_avg_n( f3->w_res(:,:,:), 0 )*1.e3

 rho = f1->rho(23,:,:,:)
 u   = f1->u  (23,:,:,:)
 rum1 = dim_avg((rho - conform_dims(dimsizes(rho),dim_avg(rho),(/0,1/)))*(u - conform_dims(dimsizes(u),dim_avg(u),(/0,1/))))

 rho = f2->rho(23,:,:,:)
 u   = f2->u  (23,:,:,:)
 rum2 = dim_avg((rho - conform_dims(dimsizes(rho),dim_avg(rho),(/0,1/)))*(u - conform_dims(dimsizes(u),dim_avg(u),(/0,1/))))

 var2(4,:,:) = -(rum2 - rum1)/30./dim_avg(f4->rho(0,:,:,:))

 var2(5,:,:) = var2(5,:,:) + var2(4,:,:)

 um = dim_avg(f4->u(0,:,:,:))

 dz2u = new((/nz,ny/),float)
 do k=1, nz-2
   dz2u(k,:) = (um(k-1,:)+um(k+1,:)-2.*um(k,:))
 end do

 dy2u = new((/nz,ny/),float)
 dy2u = 0.
 do j=1, ny-2
   dy2u(:,j) = (um(:,j-1)+um(:,j+1)-2.*um(:,j))
 end do

;--------------------------------------- process 0
var2(6,:,:) = du

 fvar = (/ var2(6,:,:), var2(5,:,:), var2(5,:,:)-var2(6,:,:), \
           var2(7,:,:)+var2(8,:,:), dy2u+dz2u   , dz2u        , \
           var2(0,:,:), var2(1,:,:), var2(2,:,:), \
           var (0,:,:), um/10., var2(10,:,:) /)

;fvar(1,:,:) = var(1,:,:)
;fvar(1,:,:) = ( var2(0,:,:) + var2(1,:,:) + var2(2,:,:) + var2(3,:,:) ) - var(1,:,:)

;--------------------------------------- cn level
 cnlev = fspan(-0.5,0.5,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlev2= (/-16,-8,-4,-2,-1,1,2,4,8,16/)

;--------------------------------------- contour
 nl = 12

tit = (/"tend","F_total","Diff","EPD","d2u","d2u_z","COR","ADV_y","ADV_z","GWD","U","W_star"/)

 wi = 60
 he = 40
 xf = 34 + (/0,1,2,0,1,2,0,1,2,0,1,2/)*70
 yf = 10 + (/0,0,0,1,1,1,2,2,2,3,3,3/)*60

do il=0, nl-1

 contour1 = cnshadeplot(wks,y,z,fvar(il,:,:),cnlev,cnfil," ")
 cnaddline(contour1,wks,y,z,fvar(il,:,:),cnlev2,0,"D",0,"T",1,"C",1)

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-90,90,20,60,False,False)
 axis_str(contour1,tit(il),"","")

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end



