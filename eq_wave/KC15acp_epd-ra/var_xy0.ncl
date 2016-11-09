load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"u","v","w","T"/)
 varnames = (/"U_GRD_GDS0_ISBL","V_GRD_GDS0_ISBL","V_VEL_GDS0_ISBL","TMP_GDS0_ISBL"/)

 yy = "2012";"1983"
 mm = "07";"05"
 p = 30
 it = 120 - 4 ; 60
 lat0 = 26

 wks = cr_wks("fig/var_xy0",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

;do it=0, 123, 4
;--------------------------------------- read data
 nv = dimsizes(varname)
 nw = 4

 fdir  = "/data11/data-arch/"+(/"JRA55","JRA55C"/)
 fname0 = fdir+"/"+yy+"/"+mm+"/"+(/"jra55","jra55c"/)+".ugrd.anal.pl."+yy+""+mm+".nc"
 fname1 = fdir+"/"+yy+"/"+mm+"/"+(/"jra55","jra55c"/)+".vgrd.anal.pl."+yy+""+mm+".nc"
 fname2 = fdir+"/"+yy+"/"+mm+"/"+(/"jra55","jra55c"/)+".vvel.anal.pl."+yy+""+mm+".nc"
 fname3 = fdir+"/"+yy+"/"+mm+"/"+(/"jra55","jra55c"/)+".tmp.anal.pl."+yy+""+mm+".nc"
print(fname3+"")

; axis
 f0 = addfile(fname0(0),"r")
 s_v0 = varnames(0)
 lat = f0->$(f0->$s_v0$!2)$({-lat0:lat0})
 lon = f0->$(f0->$s_v0$!3)$
 ny = dimsizes(lat)
 nx = dimsizes(lon)

 vvv = new((/4,ny,nx/),float)
 vvv0 = new((/4,ny,nx/),float)
 vvv1 = new((/4,ny,nx/),float)

 f0 = addfile(fname0(0),"r")
 f1 = addfile(fname1(0),"r")
 f2 = addfile(fname2(0),"r")
 f3 = addfile(fname3(0),"r")
 vvv0(0,:,:) = f0->$varnames(0)$(it,{p},{-lat0:lat0},:)
 vvv0(1,:,:) = f1->$varnames(1)$(it,{p},{-lat0:lat0},:)
 vvv0(2,:,:) = f2->$varnames(2)$(it,{p},{-lat0:lat0},:)
 vvv0(3,:,:) = f3->$varnames(3)$(it,{p},{-lat0:lat0},:)
 f0 = addfile(fname0(1),"r")
 f1 = addfile(fname1(1),"r")
 f2 = addfile(fname2(1),"r")
 f3 = addfile(fname3(1),"r")
 vvv1(0,:,:) = f0->$varnames(0)$(it,{p},{-lat0:lat0},:)
 vvv1(1,:,:) = f1->$varnames(1)$(it,{p},{-lat0:lat0},:)
 vvv1(2,:,:) = f2->$varnames(2)$(it,{p},{-lat0:lat0},:)
 vvv1(3,:,:) = f3->$varnames(3)$(it,{p},{-lat0:lat0},:)

 vvv = vvv0 ;- vvv1

 vm = dim_avg(vvv)

 vvv = vvv - conform(vvv,vm,(/0,1/))

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)
 cnlevu = ispan(-10,10,1)*1.5
 cnlevv = ispan(-10,10,1)*0.5
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
; cnfil = (/22,23,24,25,26,0,0,11,10,9,8,7/)
 cnlevh = ispan(-10,10,1)*0.3
 cnfac = (/1,0.3333,1,0.25/)

;--------------------------------------- contour
 wi = 60
 he = 25
 xf = 8 + (/0,1,0,1/)*52
 yf = 10 + (/0,0,1,1/)*30

 mstr  = new(nw,string)
 xstr  = new(nw,string)
 ystr  = new(nw,string)
 lbopt = new(nw,string)

 mstr = "";"Rossby"
; xstr = (/"","""longitude (deg)","longitude (deg)"/)
 xstr = (/"longitude (deg)"/)
 ystr = (/"latitude (deg)","","latitude (deg)",""/)
 lbtxt = "(K)"

il=0
;do il=0, nw-1

; do iv=0, nv-1
 iv=3

 if (varname(iv) .eq. "T") then
   contour1 = cnshadeplot(wks,lon,lat,vvv(iv,:,:), \
              cnlevh*cnfac(il),cnfil,lbtxt)
 end if

 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.15-0.05
  "pmLabelBarParallelPosF"   : 0.50;-0.40
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.08
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10-0.05
  "lbTitleFontHeightF" : 0.020
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,360,-25,25,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))

 cnfinalize(contour1,wks)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vvv(0:1,::2,::4)
  "vfXArray"    : lon(::4)
  "vfYArray"    : lat(::2)
 end create
 ref = cnfac*10.
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref(il)
  "vcRefLengthF"    : 0.030
  "vcRefAnnoOn"     : True
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(il),yf(il),wi,he)
 axis_range(vec1,0,360,-25,25,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

;end do

 frame(wks)

;end do

end

