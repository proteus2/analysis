load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = "era-int_f"
 varname = (/"epd","f_y","f_z"/)

 yy = 2007
 mm = 6
 y1 = 1979
 lat1 = 25.
 z1p = 17
 z2p = 39
 z1 = z1p - 1.5

 iwvs = (/3,2,1/)
 ; 0 - 6   : K, R, E-MRG, W-MRG, E-IG, W-IG, mig.-tide
 swv = (/"k","r","emrg","wmrg","eig","wig","tid"/)
 scl = (/0.8,0.8,0.8/)

 wks = cr_wks("fig/epf-mrg_yp-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 it = (yy - y1)*12 + (mm - 1)

do iw=0, 2
  iwv = iwvs(iw)
;--------------------------------------- read data
 nv = dimsizes(varname) + 1

 z2 = z2p + 4

 p1 = numeric2int( 70.*exp((18.6-z1)/6.6), 1 )
 p2 = numeric2int( 70.*exp((18.6-z2)/6.6), 1 )

 f2 = addfile("/data18/kyh/dat/tem/ERA-I-nr/era-int_f.tem_ypt.1979-2013.nc","r")
 ut  = f2->u_tend(it,{p1:p2},{-lat1:lat1})
 lat = f2->lat({-lat1:lat1})
 ht0 = tofloat( 18.6e3 + 6.6e3*log( 70./ f2->p ) )
 ht0!0 = "z"
 ht0&z = ht0
 ht = ht0({z1*1.e3:z2*1.e3})
 ny = dimsizes(lat)
 nz = dimsizes(ht )

 fname = "/data18/kyh/dat/clim/ERA-I-nr/era-int_f.u_ypt.1979-2012.nc"
 f1 = addfile(fname,"r")
 tmp = (/ tofloat( f1->u(it,:,{-lat1:lat1}) ) /)
 tmp!0 = "z"
 tmp&z = tofloat(ht0)
 u = tmp({z1*1.e3:z2*1.e3},:)


 f2a = addfile("/data18/kyh/dat/tem/ERA-I-nr/epf/era-int_f.epfsa_yz_recon0.nc","r")
 f2b = addfile("/data18/kyh/dat/tem/ERA-I-nr/epf/era-int_f.epf_yz_recon0.nc","r")

 if (iw .eq. 0) then
   wv  = new((/nv-1,nz,ny/),float)
   wv2 = new((/nv-1,nz,ny/),float)
   wv3 = new((/nv-1,nz,ny/),float)
 end if
 ; 0 - 6   : K, R, E-MRG, W-MRG, E-IG, W-IG, mig.-tide, MRG

 do iv=0, nv-2
   if (iwv .le. 5) then
     wv3(iv,:,:) = f2a->$varname(iv)$(it-1,iwv+2, \
        {p1:p2},{-lat1:lat1})
     if (iwv .eq. 1) then
       wv3(iv,:,:) = wv3(iv,:,:) + f2a->$varname(iv)$(it-1,iwv+12, \
          {p1:p2},{-lat1:lat1})
     end if
     if (iwv .ge. 4) then
       wv3(iv,:,:) = wv3(iv,:,:) + f2b->$varname(iv)$(it-1,iwv+4, \
          {p1:p2},{-lat1:lat1}) \
                                 + f2a->$varname(iv)$(it-1,iwv+12, \
          {p1:p2},{-lat1:lat1})
     end if
   end if
   if (iwv .eq. 6) then
     wv3(iv,:,:) = f2b->$varname(iv)$(it-1,2, \
        {p1:p2},{-lat1:lat1}) + \
                   f2b->$varname(iv)$(it-1,3, \
        {p1:p2},{-lat1:lat1})
   end if
 end do
 if ( iwv .le. 5 .and. mod(iwv,2) .eq. 0 ) then
   wv3(1:2,:,:) = wv3(1:2,:,:)*(-1.)  ; for eastward waves
 end if

 if (iw .eq. 0) then
   wv = wv3
   wv3 = 0.
 end if
 if (iw .eq. 1) then
   wv2 = wv3
   wv3 = 0.
 end if

end do  ; iw

;--------------------------------------- process 0
 vari = new((/nv,nz,ny/),float)
 vari(0,:,:) = (/ u /)
 vari(1:,:,:) = wv

 vari2 = new((/nv,nz,ny/),float)
 vari2(0,:,:) = (/ u /)
 vari2(1:,:,:) = wv2

 vari3 = new((/nv,nz,ny/),float)
 vari3(0,:,:) = (/ u /)
 vari3(1:,:,:) = wv3

 vari(1,:,:) = vari(1,:,:)*30.  ; m/s/month
 vari(2:3,:,:) = vari(2:3,:,:) / 6371.e3 * 1.e3  ; mPa
 vari2(1,:,:) = vari2(1,:,:)*30.  ; m/s/month
 vari2(2:3,:,:) = vari2(2:3,:,:) / 6371.e3 * 1.e3  ; mPa
 vari3(1,:,:) = vari3(1,:,:)*30.  ; m/s/month
 vari3(2:3,:,:) = vari3(2:3,:,:) / 6371.e3 * 1.e3  ; mPa

;--------------------------------------- cn level
 cnlev = fspan(-5,5,21)*0.5
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 nfig = 2
 nu   = 3
 nstrd= 4
 ustt = 1

 wi = 22
 he = 33
 xf = 20 + (/0,1,2,3,4/)*30
; yf = 10 + (/0,1,2,3/)*42 + (/0,0,1,1/)*6
 yf = 10 + (/0,1,2/)*38

 r_yz = 3.141592/180.*6371.
 vari (2,:,:) = vari (2,:,:)/r_yz
 vari2(2,:,:) = vari2(2,:,:)/r_yz
 vari3(2,:,:) = vari3(2,:,:)/r_yz
 do k=0, nz-1
   vari (2:3,k,:) = vari (2:3,k,:)*exp(ht(k)/(2.*6.6e3))
   vari2(2:3,k,:) = vari2(2:3,k,:)*exp(ht(k)/(2.*6.6e3))
   vari3(2:3,k,:) = vari3(2:3,k,:)*exp(ht(k)/(2.*6.6e3))
 end do

do iw=0, 0;1;2
  iwv = iwvs(iw)

  iwv = iwvs(iw)

 mstr  = (/"W","E"/)+"-MRG"
 xstr  = (/"lat (deg)",""/)
 ystr  = (/"z (km)","","","",""/)
 refann = (/False,False,False,False,True/)

  vvv = vari
  if (iw .eq. 1) then
    vvv = vari2
  end if
  if (iw .eq. 2) then
    vvv = vari3
  end if

 contour1 = cnshadeplot(wks,lat,ht/1.e3,vvv(1,:,:),cnlev,cnfil,0)
 contour0 = cnlineplot(wks,lat,ht/1.e3,vvv(0,:,:),ispan(-50,50,5),0)
 cnline_att(contour0,"D",(/0,0,16/),"T",(/0.6,1.8,0.7/),"C",(/1,colind,1/))

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vvv(2:3,:,::4)
  "vfXArray"    : lat(::4)
  "vfYArray"    : ht/1.e3
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
;  "vcLineArrowThicknessF" :  vthc
;  "vcLineArrowColor"      :  vcol
  "vcRefMagnitudeF" : scl(iw)
  "vcRefLengthF"    : 0.05
  "vcGlyphStyle"              : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.12
  "vcFillArrowWidthF"         : 0.1
  "vcPositionMode"            : "ArrowTail"
  "vcRefAnnoOn"          : refann(0)
  "vcRefAnnoArrowAngleF" : 90.
  "vcRefAnnoFontHeightF" : 0.024
  "vcRefAnnoString1On"   : False
;  "vcRefAnnoString1"     : "1 mPa"
;  "vcRefAnnoString2On"   : False
  "vcRefAnnoString2"     : scl(iw)+" mPa"
  "vcRefAnnoPerimSpaceF" : 0.6
  "vcRefAnnoSide" : "Right"
  "vcRefAnnoOrthogonalPosF" : 0.1 ;0.15
  "vcRefAnnoParallelPosF"   : 0.1 ;1. - nn*0.1
 end create

 set_size(contour1,xf(iw),yf(0),wi,he)
 axis_range(contour1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour1,mstr(iw),xstr(0),ystr(iw))

 set_size(contour0,xf(iw),yf(0),wi,he)
 axis_range(contour0,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour0,"","","")

 set_size(vec1,xf(iw),yf(0),wi,he)
 axis_range(vec1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)
 cnfinalize(contour0,wks)

 draw(contour1)

 overlay(contour0,vec1)
 draw(contour0)
 NhlRemoveOverlay(contour0,vec1,False)

end do  ; iw

 label = create "lab" labelBarClass wks
  "lbOrientation" : "horizontal"
  "lbBoxCount" : dimsizes(cnfil)
  "lbFillColors" : cnfil
  "lbMonoFillPattern" : True
  "lbLabelsOn" : True
  "lbLabelStrings" : cnlev
  "lbAutoManage" : False
  "lbLabelAutoStride" : True
;  "lbLabelPosition" : "Bottom" ;"Right"
  "lbLabelAlignment" : "InteriorEdges" ;"BoxCenters" ;"ExternalEdges"
;  "lbLabelAngleF" : 0.
  "lbPerimOn" : False
  "lbTitleString" : "(m s:S:-1:N: month:S:-1:N:)"
  "lbTitlePosition" : "Bottom" ;"Right" ;"Top"
  "lbTitleDirection" : "Across"
 end create
 set_size_lb(label,xf(0)-9,yf(0)+he+10,40,4)

 draw(label)

 frame(wks)

end

