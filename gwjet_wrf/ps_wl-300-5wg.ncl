load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 lev = 8000
 nw = 5
 dphi = 15.
 phi1 = (/-75,-30/)
 phi2 = (/-15,30/)
 phi3 = (/45, 75/)
 phi4a = (/0, 90/)
 phi4b = (/-180, -90/)
 phi5a = (/105,165/)
 phi5b = (/-75,-15/)
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_wl-300"+var1name+"-"+nw+"wg",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-1.nc"
 fname2 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-2.nc"

 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 nx = 432
 ny = 302
 dx = 1./10.8
 dy = 1./10.8
 dt = 1./12.
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy

 kh = f1->kh / 111.
 nkh = dimsizes(kh)

 vs1 = dim_sum( f1->psd_phi_kh(:,{phi1(0):phi1(1)}) )
 vs2 = dim_sum( f1->psd_phi_kh(:,{phi2(0):phi2(1)}) )
 vs3 = dim_sum( f1->psd_phi_kh(:,{phi3(0):phi3(1)}) )
 vs4 = dim_sum( f2->psd_phi_kh(:,{phi4a(0):phi4a(1)}) ) + \
       dim_sum( f2->psd_phi_kh(:,{phi4a(0):phi4a(1)}) )
 vs5 = dim_sum( f2->psd_phi_kh(:,{phi5a(0):phi5a(1)}) ) + \
       dim_sum( f2->psd_phi_kh(:,{phi5a(0):phi5a(1)}) )
 vs = (/vs1,vs2,vs3,vs4,vs5/) * 111. * dphi
; vs(:,:,1:(nx-1)/2) = vs(:,:,1:(nx-1)/2)*2.
; vs = vs*(nx*dx)*(ny*dy)
; vs(:,:,0) = 1.e32

;--------------------------------------- process 0
 vsh = vs*log(10.)*conform(vs,kh,1)

 vsh = where(vsh .eq. 0., 1.e-30, vsh)

 wl = where(kh .ne. 0., 1./kh, 999)

;--------------------------------------- cn level
 lcol = (/18,16,13,10,1/)

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "horizontal wavelength [km]"
ystr = "PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"

 il = 0

 plot = xylineplot(wks,wl,vsh*1.e5,"D",0,"T",1,"C",lcol)

 set_size(plot,xf(il),yf(il),wi,he)
 axis_range(plot,999,10,-999,999,True,False)
 axis_str(plot,mstr(il),xstr(il),ystr(il))

 draw(plot)

; dr_txt(wks,subl(il),xf(0)-8,yf(il),25,6,0.)

 frame(wks)

end



