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

 wks = cr_wks("fig/ps_wlzt-300"+var1name,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

do nn=0, 0 ;1
;
; it1 = ispan(144,287,1) + nn*144
; it2 = ispan(0,143,1) + nn*144
 it1 = ispan(0  ,11,1)
 it2 = ispan(0  ,11,1)
 
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
; cnlev = (/ (/0.3,0.5,0.7,1,1.5,2,2.5,3.5,5,7,10/), \
;            (/0.3,0.5,0.7,1,1.5,2,2.5,3.5,5,7,10/) /)
; cnfil = (/0,6,7,10,12,14,15,16,17,19,20,22/)
 cnlev = (/ (/0.3,0.5,0.7,1,1.5,2,2.5,3.5,5,7/), \
            (/0.3,0.5,0.7,1,1.5,2,2.5,3.5,5,7/) /)
 cnfil = (/0,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 40
 he = 30
 xf = 12 + (/0,1/)*(wi+6)
 yf = 10 + (/0,1/)*(he+6)

mstr = ""
;xstr = (/"","horizontal wavelength [km]"/)
xstr = (/"horizontal wavelength [km]",""/)
ystr = (/"z [km]",""/)
;lbopt = (/ (/"",""/),(/"","PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"/) /)
lbopt = (/ (/"","PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"/), (/"",""/) /)

do il=0, 1

 contour1 = cnshadeplot(wks,wl,z/1.e3,vsh(il,:,:)*1.e5,cnlev(il,:),cnfil,lbopt(nn,il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 0.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
;  "pmLabelBarHeightF"        : 0.10
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size(contour1,xf(il),yf(nn),wi,he)
 axis_range(contour1,2001,10,-999,999,True,False)
 axis_str(contour1,mstr,xstr(nn),ystr(il))

 cnfinalize(contour1,wks)

 draw(contour1)

end do

end do  ; nn

 dr_txt(wks,"D1",xf(0)+wi/2,yf(0)-4.0,22,7.0,0.)
 dr_txt(wks,"D2",xf(1)+wi/2,yf(0)-4.0,22,7.0,0.)
; dr_txt(wks,"Day 6, 12::00-24::00",xf(0)+12,yf(0)+1.8,22,5.5,0.)
; dr_txt(wks,"Day 7, 00::00-12::00",xf(0)+12,yf(1)+1.8,22,5.5,0.)
; dr_txt(wks,"Day 7, 00::00-12::00",xf(1)+wi-12,yf(0)+1.8,22,5.5,0.)
; dr_txt(wks,"Day 7, 12::00-24::00",xf(1)+wi-12,yf(1)+1.8,22,5.5,0.)

 frame(wks)

end

