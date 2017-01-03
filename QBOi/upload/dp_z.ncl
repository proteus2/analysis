load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"L60CTL","L60CGW"/)
 expcode = (/"ui"    ,"uj"    /)

 latb = 15.
 y1 = 1979
 y2 = 2006

 ; 100--0.4 hPa
 p1 = 200
 p2 = 0.2
 z1 = 15
 z2 = 60

 wks = cr_wks("fig/dp_z",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nyr = y2-y1+1
 nt = nyr*12

 dname = "$DATD/"+expname+"/dp"
 fname = new(nt*ne,string)

 i = 0
 do ie=0, ne-1
 do m=1, 12
   mm = sprinti("%0.2i",m)
   fname(i:i+nyr-1) = dname(ie)+"/"+ispan(y1,y2,1)+"/uan"+expcode(ie)+  \
                      ".dp_yzt."+ispan(y1,y2,1)+"."+mm+".nc"
   i = i + nyr
 end do
 end do

 f1 = addfiles(fname,"r")
 ListSetType(f1,"join")
 p_m  = dim_avg_n( dim_avg( f1[:]->p_m (:,{z1*1.e3:z2*1.e3},{-latb:latb}) ), 0 )
 p2_m = dim_avg_n( dim_avg( f1[:]->p2_m(:,{z1*1.e3:z2*1.e3},{-latb:latb}) ), 0 )

 f0 = addfile(fname(0),"r")
 zm = f0->z({z1*1.e3:z2*1.e3}) * 1.e-3
 nz = dimsizes(zm)

 sd = sqrt( p2_m - p_m^2 )

 print(p_m+"  "+sd)


 f2 = addfile("/hippo0/HG2AMIP/"+expname(0)+"/invariant/L60_z.nc","r")
 zt = f2->hybrid_ht({min(zm)*1.e3:max(zm)*1.e3}) * 1.e-3

 coef_rho0 = ( zm(1:nz-1) - zt(0:nz-2) )/  \
             ( zm(1:nz-1) - zm(0:nz-2) )
 coef_rho1 = 1. - coef_rho0
 print(coef_rho0)

 p_th = exp( log(p_m(0:nz-2))*coef_rho0 + log(p_m(1:nz-1))*coef_rho1 )

 print(p_th)

;--------------------------------------- cn level

;--------------------------------------- ctr
 wi = 60
 he = 25
 xf = 20
 yf = 10

 xstr = "p [hPa]"
 ystr = "z [km]"

 xmin = min(p_m) - 2.*sd(minind(p_m))

 yv = zm
 nsd = 2*sd
 tcks = 0.8

 dum1 = new((/3,nz/),graphic)

 plot1 = xylineplot(wks,p_m,yv,"D",0,"T",1,"C",-1)
 xyaddmarker(plot1,1,"S",3.0,"T",0.6,"C",2)

 plres = True
 plres@gsLineThicknessF = 0.8
 plres@gsFillColor      = 1
 do k=0, nz-1
   dum1(0,k) = gsn_add_polyline(wks,plot1,(/p_m(k)-nsd(k),p_m(k)+nsd(k)/),  \
                                (/yv(k),yv(k)/),plres)
   dum1(1,k) = gsn_add_polyline(wks,plot1,(/p_m(k)-nsd(k),p_m(k)-nsd(k)/),  \
                                (/yv(k)-tcks,yv(k)+tcks/),plres)
   dum1(2,k) = gsn_add_polyline(wks,plot1,(/p_m(k)+nsd(k),p_m(k)+nsd(k)/),  \
                                (/yv(k)-tcks,yv(k)+tcks/),plres)
 end do

 set_size(plot1,xf,yf,wi,he)
 axis_range(plot1,p2*0.999,p1*1.001,z1,z2,True,False)
 axis_str(plot1,"Model-level mean (red) and 2SD range (black) of pressure over 15N-15S",xstr,ystr)

 draw(plot1)

 frame(wks)

end

