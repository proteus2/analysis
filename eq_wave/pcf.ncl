load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 dx = 1.25/6.
 lx = 90./6.
 xrng = 20.

 wks = cr_wks("pcf",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 pi = 3.141592

;--------------------------------------- read data

;--------------------------------------- process 0
 x = fspan(-lx,lx,toint(2.*lx/dx)+1)
 nx = dimsizes(x)

 pcf = new((/8,nx/),float)
 pcfsum = new(8,float)

 x2 = x^2
 gau = exp(-0.25*x2)
 pcf(0,:) = 1.
 pcf(1,:) = x
 pcf(2,:) = x2 - 1.
 pcf(3,:) = x*(x2 - 3.)
 pcf(4,:) = x2*(x2 - 6.) + 3.
 pcf(5,:) = x*(x2*(x2 - 10.) + 15.)
 pcf(6,:) = x2*(x2*(x2 - 15.) + 45.) - 15.
 pcf(7,:) = x*(x2*(x2*(x2 - 21.) + 105.) - 105.)

 pcf = pcf*conform_dims((/8,nx/),gau,1)
 pcf!1 = "x"
 pcf&x = x*6.
 pcfsum = (/1.,1.,2.,6.,24.,120.,720.,5040./)
 pcfsum = (pcfsum*sqrt(2.*pi))/dx
 pcfsum1 = dim_sum(pcf(:,{-xrng:xrng})^2)

 pcfco = new((/8,8/),float)
 do j=0, 7
 do i=0, 7
   pcfco(j,i) = sum(pcf(j,{-xrng:xrng})*pcf(i,{-xrng:xrng}))/sqrt(pcfsum1(j)*pcfsum1(i))
 end do
 end do

 print("theoretical  rel. err.")
 print("  "+pcfsum+"   "+(1.-pcfsum1/pcfsum))

;--------------------------------------- cn level
 lcol = (/1,2,3,4,6,13/)

;--------------------------------------- plot
 nl = 1

 wi = 75
 he = 40
 xf = 34
 yf = 10 + (/0,1,2/)*60

 plot = xylineplot(wks,x*6.,0.,"D",0,"T",0.5,"C",1)
 do k=0, 4
   xyaddline(plot,x*6.,pcf(k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot,xf,yf(0),wi,he)
 axis_range(plot,-45,45,-4,4,False,False)
 axis_str(plot,"","","")

 draw(plot)


 plot2 = xylineplot(wks,ispan(-1,8,1),0.,"D",0,"T",0.5,"C",1)
 do k=0, 4
   xyaddline(plot2,ispan(0,7,1),pcfco(k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot2,xf,yf(1),wi,he)
 axis_range(plot2,-0.5,7,-0.6,1.2,False,False)
 axis_str(plot2,"","","")
 setvalues plot2
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tmXBTickStartF"    : 0
 end setvalues

 draw(plot2)


 frame(wks)

end



