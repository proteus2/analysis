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
 dz = 250
 nz = 4*10+1
 dphi = 15.
 var1name = "prt_"+var1
 it1 = ispan(144,431,1)
 it2 = ispan(0,287,1)
 
 wks = cr_wks("fig/ps_wlz2-300"+var1name,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 z = ispan(0,nz-1,1)*dz

do k=0, nz-1

 fname0 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_wn__z"+sprinti("%5.5i",z(k))+"-"+ispan(1,2,1)+"a.nc"

 f01 = addfile(fname0(0),"r")
 f02 = addfile(fname0(1),"r")

 kh = f01->kh / 111.
 nkh = dimsizes(kh)

 vs01 = dim_avg_n( dim_sum( f01->psd_phi_kh(it1,:,:) ), 0 )
 vs02 = dim_avg_n( dim_sum( f02->psd_phi_kh(it2,:,:) ), 0 )

 if (k .eq. 0) then
   vs = new((/2,nz,nkh/),float)
 end if

 vs(:,k,:) = (/vs01,vs02/) * 111. * dphi * exp(-z(k)/7.e3)

end do

;--------------------------------------- process 0
 vsh = vs*log(10.)*conform(vs,kh,2)

 vsh = where(vsh .eq. 0., 1.e-30, vsh)

 wl = where(kh .ne. 0., 1./kh, 999)

;--------------------------------------- cn level
 cnlev = (/ ispan(0,10,1), ispan(0,10,1) /)*0.25 - 6
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 50
 he = 30
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "horizontal wavelength [km]"
ystr = (/"z [km]",""/)
lbopt = "PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"

do il=0, 1

 contour1 = cnshadeplot(wks,wl,z/1.e3,log10(vsh(il,:,:)),cnlev(il,:),cnfil,lbopt)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -0.8;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.5
;  "pmLabelBarHeightF"        : 0.10
;  "lbLabelStride"      : 2
;  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.06
;  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,2001,10,-999,999,True,False)
 axis_str(contour1,mstr,xstr,ystr(il))

 cnfinalize(contour1,wks)

 draw(contour1)

end do

 frame(wks)

end

