load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 load "./param0.ncl"
 load "./param1.ncl"
 load "./param_ra.ncl"
 load "./ftn_get_color_S-RIP.ncl"
 load "./ftn_read_spec.ncl"

 p_max = 100
 p_min = 2

 ra = (/"ERAI","MERRA2","JRA55"/)
 ra_m = ra+"_m"

 p_all = param_ra_plevs

 h_scale = param0_h_scale

 ranames = get_ra_names(ra)
 ranames_m = get_ra_names(ra_m)
 
 wks = cr_wks("fig/dz-p_z",600,0)
 icol_s = get_col_srip(wks)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 ne = dimsizes(ra)
 nfig = 2

 yy = new((/2,ne,12/),float,1.e20)
 yy = yy@_FillValue
 dz = new((/nfig,ne,12/),float,1.e20)
 dz = dz@_FillValue

do ie=0, ne-1

 ir = toint( ranames(ie,0) )
 ir0 = toint( ranames_m(ie,0) )

 p = p_all(ir,ind( p_all(ie,:) .le. p_max .and. p_all(ie,:) .ge. p_min ))
 p0 = p_all(ir0,ind( .not. ismissing(p_all(ir0,:)) ))

 z  = 18.6 + h_scale*log(70./p )
 z0 = 18.6 + h_scale*log(70./p0)

 iz1 = ind_nearest_coord(z,z0,0)
 z1 = z0(iz1)

 z2 = z1  ; dummy for declaration
 do k=0, dimsizes(z)-1
   z0x = where(z0 .eq. z1(k), z0@_FillValue, z0)
   iz2 = ind_nearest_coord(z(k),z0x,0)
   z2(k) = z0(iz2)
   delete(z0x)
 end do

 yy(0,ie,:dimsizes(z1)-1) = p
 yy(1,ie,:dimsizes(z1)-1) = z*1.e-3
 dz(0,ie,:dimsizes(z1)-1) = abs(z1 - z)
 dz(1,ie,:dimsizes(z1)-1) = abs(z1 - z)/abs(z2 - z1)

 delete( [/p,p0,z,z0,iz1,z1,z2/] )

end do  ; ie

;--------------------------------------- cn level
 ltck = (/1,1,1,1,1/)*0.8
 ldsh = (/0,0,0,0,0/)*16
; lcol = (/ icol_s@ei , icol_s@m , icol_s@m2 , icol_s@c , icol_s@j5 /)
 lcol = (/ icol_s@ei , icol_s@m2 , icol_s@j5 /)

;--------------------------------------- ctr
 wi = 30
 he = 30
 xf = 12 + (/0,1/)*(wi+10)
 yf = 10

 mstr = ""
 xstr = (/"distance from the nearest level [m]","distance relative to the grid spacing"/)
 ystr = "p [hPa]"

 ytckl = (/100,70,50,30,20,10,7,5,3,2,1/)
 ytckv = 18.6 + h_scale*1.e-3*log(70./ytckl)
 tmpx = (/40,4/)
 ymntckv = 18.6 + h_scale*1.e-3*log(70./tmpx)
 ymax = 18.6 + h_scale*1.e-3*log(70./p_min)
 ymin = 18.6 + h_scale*1.e-3*log(70./p_max)


do ip=0, nfig-1

 plot1 = xylineplot(wks,dz(ip,:,:),yy(1,:,:),"D",ldsh,"T",ltck,"C",lcol)
 xyaddmarker(plot1,16,"S",1.5,"T",0.6,"C",lcol)

 set_size(plot1,xf(ip),yf,wi,he)
 axis_range(plot1,-999,999,ymin,ymax,False,False)
; axis_range(plot1,999,0,ymin,ymax,False,False)
 axis_str(plot1,"",xstr(ip),ystr)
 setvalues plot1
  "tmYLMode" : "Explicit"
  "tmYLValues" : ytckv
  "tmYLLabels" : ytckl
  "tmYLMinorValues" : ymntckv
  "tmXBMaxTicks" : 7
 end setvalues

 draw(plot1)


; dr_txt2(wks,mstr,xf(ip)+1,yf+2,22,4.5,0., 1,1)

end do

 lgtit = ranames(:,1)
 lgtit = where(lgtit .eq. "MERRA2", lgtit+" (& MERRA)", lgtit)
 lg = lgplot(wks,lgtit,"D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(0),yf+he+7,wi*0.5,2.8*3)
 draw(lg)


 frame(wks)

end

