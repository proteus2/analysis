load "$UTIL/ncl/header.ncl"
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

beta = 2.*7.292e-5/6370.e3
pi = 3.141592
kwn = 2.*pi/90.  ; [rad/deg]
ome = 2.*pi/2.   ; [rad/day]
koo = kwn/ome*(86400./111.e3)  ; [s/m]

x = fspan(0,360-3,120)
;y = fspan(-15,15,31)
y = fspan(-90,90,181)
nx = dimsizes(x)
ny = dimsizes(y)
ym = y*111.e3
d2xm = 2.*(x(1) - x(0))*111.e3
d2ym = 2.*(y(1) - y(0))*111.e3

u = new((/ny,nx/),float)
u = conform_dims(dimsizes(u),exp(-0.5*beta*koo*ym^2),0)*  \
    conform_dims(dimsizes(u),cos(kwn*x),1)
u = u*10.01

dudy0 = u*conform_dims(dimsizes(u),-beta*koo*ym,0)

phi0 = ilapsF(-dudy0,0)

u_phi0 = phi0
u_phi0(1:ny-2,:) = -(phi0(2:,:) - phi0(:ny-3,:))/d2ym
u_phi0(0,:) = 0.
u_phi0(ny-1,:) = 0.

v_phi0 = phi0
v_phi0(:,1:nx-2) = (phi0(:,2:) - phi0(:,:nx-3))/d2xm
v_phi0(:,0) = (phi0(:,1) - phi0(:,nx-1))/d2xm
v_phi0(:,nx-1) = (phi0(:,0) - phi0(:,nx-2))/d2xm


wks = cr_wks("zrot",600,0)

 contour1 = cnlineplot(wks,x,y,u,ispan(-10,10,1),0)
 cnline_att(contour1,"D",(/0,0,5/),"T",0.8,"C",(/1,-1,1/))
; setvalues contour1
;  "cnLineLabelFontHeightF" : 0.022
; end setvalues

 set_size(contour1,10,10,60,15)
 axis_range(contour1,0,360,-15,15,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 draw(contour1)

 contour2 = cnlineplot(wks,x,y,u_phi0,ispan(-10,10,1),0)
 cnline_att(contour2,"D",(/0,0,5/),"T",0.8,"C",(/1,-1,1/))
 cnaddline(contour2,wks,x,y,v_phi0,ispan(-10,10,1),0,  \
           "D",(/0,0,5/),"T",0.4,"C",(/1,-1,1/))
; setvalues contour2
;  "cnLineLabelFontHeightF" : 0.022
; end setvalues

 set_size(contour2,10,30,60,15)
 axis_range(contour2,0,360,-15,15,False,False)
 axis_str(contour2,"","","")

 cnfinalize(contour2,wks)

 draw(contour2)

 contour3 = cnlineplot(wks,x,y,u - u_phi0,ispan(-10,10,1),0)
 cnline_att(contour3,"D",(/0,0,5/),"T",0.8,"C",(/1,-1,1/))
 cnaddline(contour3,wks,x,y,-v_phi0,ispan(-10,10,1),0,  \
           "D",(/0,0,5/),"T",0.4,"C",(/1,-1,1/))
; setvalues contour3
;  "cnLineLabelFontHeightF" : 0.022
; end setvalues

 set_size(contour3,10,50,60,15)
 axis_range(contour3,0,360,-15,15,False,False)
 axis_str(contour3,"","","")

 cnfinalize(contour3,wks)

 draw(contour3)


 frame(wks)

