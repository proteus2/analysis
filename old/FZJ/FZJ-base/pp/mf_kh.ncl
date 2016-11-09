load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("mf_kh",600,0)
 gsn_define_colormap(wks,"gsdtol")
 icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

 f01 = addfile("/data4/atmosdata/T_SPECTRUM/P1/n_profile_ave_p1.nc","r")
 f02 = addfile("/data4/atmosdata/T_SPECTRUM/P1/u_profile_ave_p1.nc","r")

 y0  = f01->height({20000:}) / 1.e3

 x0      = new((/2,dimsizes(y0)/),float,1.e32)
 x0(0,:) = f01->n_ave({20000:})
 x0(1,:) = f02->u_ave({20000:})


   f1 = addfile("../res/temp2/mf_kh/mfx_kh_p1_f1.nc","r")
   y    = f1->height/1.e3
   x    = f1->kh

   var1 = f1->mfx_up_e*1.e3
   var2 = f1->mfx_dn_e*1.e3

;--------------------------------------- line attribute
 cnlev1 = (/2.5,5,10,15,20,25,30/)*2
 cnlev2 = (/-1.5,-1,-0.5,-0.25/)*2
 cnfil1 = (/0,28,24,20,16,12,8,4/)
 cnfil2 = (/13,18,23,28,0/)
;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 10 + (/0,1,2/)*50 + (/0,6,6/)
 yf = 10

 mstr = (/"","MFx_T_up","MFx_T_down"/)


 plot = xylineplot(wks,x0(1,:),y0,"D",1,"T",1.2,"C",1)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,-999,999,-999,45,False,False)
 axis_str(plot,"","U [m s:S:-1:N:]","z [km]")
 setvalues plot
  "tmXTOn" : False
 end setvalues

 plot2 = xylineplot(wks,x0(0,:)*100.,y0,"D",0,"T",1.2,"C",1)

 set_size(plot2,xf(0),yf,wi,he)
 axis_range(plot2,1.5,3,-999,45,False,False)
 axis_str(plot2,"","N [10:S:-2:N: s:S:-1:N:]","")
 setvalues plot2
  "tmXBOn" : False
  "tmXTLabelsOn" : True
  "tiXAxisSide" : "Top"
 end setvalues


 draw(plot)
 draw(plot2)


 contour1 = cnshadeplot(wks,x*1.e6,y,var1*1.e-4, \
           cnlev1,cnfil1,"[10 mPa / cycle km:S:-1:N:]")

 set_size(contour1,xf(1),yf,wi,he)
 axis_range(contour1,-999,999,20,999,False,False)
 axis_str(contour1,mstr(1),"k:B:h:N: [10:S:-3:N: cycle km:S:-1:N:]","z [km]")

 cnfinalize(contour1,wks)

 contour2 = cnshadeplot(wks,x*1.e6,y,var2*1.e-4, \
           cnlev2,cnfil2,"[10 mPa / cycle km:S:-1:N:]")

 set_size(contour2,xf(2),yf,wi,he)
 axis_range(contour2,-999,999,20,999,False,False)
 axis_str(contour2,mstr(2),"k:B:h:N: [10:S:-3:N: cycle km:S:-1:N:]","")

 cnfinalize(contour2,wks)


 draw(contour1)
 draw(contour2)


 frame(wks)

end

