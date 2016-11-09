load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 wks = cr_wks("./ft2",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fff0 = systemfunc("ls ../dat/arch/uj/pj/uanuja.pj_195301.nc")
 fff3 = systemfunc("ls ../dat/tem/1953/uanuj.tem_yzt.1953.01.nc")
 fff1 = systemfunc("ls /data4/kyh/epf/1953/uanuj.epf_koyz.1953.01.nc")

 f0 = addfiles(fff0,"r")
 f3 = addfiles(fff3,"r")
; f31= addfile("../dat/tem/1999/uanuj.tem_yzt.1999.12.nc","r")
; f32= addfile("../dat/tem/2001/uanuj.tem_yzt.2001.01.nc","r")
 f1 = addfiles(fff1,"r")
; ListSetType (f1, "join")

 y = f1[0]->lat({-25:25})
 z = f1[0]->z({18590:80000})/1.e3
 ny = dimsizes(y)
 nz = dimsizes(z)

 var = new((/6,nz,ny/),float)
 var = 0.

 var(0,:,:) = dim_avg_n( f3[:]->epd(:,{18590:80000},{-25:25}), 0 )
; var(1,:,:) = 0.5*var(0,:,:) + \
;              0.25*dim_avg_n( f31->epd(:,{18590:80000},{-25:25}), 0 ) + \
;              0.25*dim_avg_n( f32->epd(:,{18590:80000},{-25:25}), 0 )
; var(1,:,:) = ( var(0,:,:)*12. + \
;                dim_avg_n( f31->epd(:,{18590:80000},{-25:25}), 0 )*0.5 + \
;                dim_avg_n( f32->epd(:,{18590:80000},{-25:25}), 0 )*0.5 ) / 13.

 do j=0, dimsizes(fff1)-1
   var(2,:,:) = var(2,:,:) + dim_sum_n( f1[j]->epd({18590:80000},{-25:25},:,:), (/2,3/) )
 end do
 var(2,:,:) = var(2,:,:)/tofloat(dimsizes(fff1))

 var(3,:,:) = var(2,:,:) - var(0,:,:)

; var2(7,:,:)  = dim_avg_n( f3->epd_y(:,:,:), 0 )
; var2(8,:,:)  = dim_avg_n( f3->epd_z(:,:,:), 0 )

; um = dim_avg(f0->u(0,:,:,:))

;--------------------------------------- process 0
 fvar = (/ var(0,:,:), var(1,:,:), var(2,:,:), var(3,:,:) /) ;*0.001

;--------------------------------------- cn level
 cnlev = fspan(-0.5,0.5,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlev2= (/-16,-8,-4,-2,-1,1,2,4,8,16/)

;--------------------------------------- contour
 nl = 4

tit = (/"EPD","EPD3","FC","DIFF"/)

 wi = 60
 he = 40
 xf = 34 + (/0,1,2,0,1,2,0,1,2,0,1,2/)*70
 yf = 10 + (/0,0,0,1,1,1,2,2,2,3,3,3/)*60

do il=0, nl-1

 contour1 = cnshadeplot(wks,y,z,fvar(il,:,:),cnlev,cnfil," ")
 cnaddline(contour1,wks,y,z,fvar(il,:,:),cnlev2,0,"D",0,"T",1,"C",1)

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-25,25,19,60,False,False)
 axis_str(contour1,tit(il),"","")

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end



