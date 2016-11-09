load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "uj"
 varname = (/"epd","f_y","f_z"/)

 y0 = 1953
 y1 = 1980
 y2 = 2005
 lat1 = 25.
 z1p = 17
 z2p = 39
 z1 = z1p - 1.5

 iwvs = (/3,1/)
 ; 0 - 5   : K, R, E-MRG, W-MRG, E-IG, W-IG
 swv = (/"k","r","emrg","wmrg","eig","wig"/)
 scl = (/1.,2.,0.1,0.2,1.,1./)

 wks = cr_wks("fig/epf-mrg_yz-"+expname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
do iw=0, 0
 iwv = iwvs(iw)

 n1 = (y1 - y0)*12
 n2 = (y2 - y0 + 1)*12 - 1
 nt = (y2 - y1 + 1)*12
 nv = dimsizes(varname) + 1

 z2 = z2p + 4

 fname = "$DATD/L60CGW-t/xm.uan"+expname+"a.pj.nc"
 f1 = addfile(fname,"r")
 lat = f1->latitude({-lat1:lat1})
 ht = f1->hybrid_ht({z1*1.e3:z2*1.e3})
 ny = dimsizes(lat)
 nz = dimsizes(ht )
 u  = f1->u(n1:n2,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("$DATD/L60CGW-t/tem/uan"+expname+".tem-wc_yzt.nc","r")
 ut  = f2->u_tend(n1:n2,{z1*1.e3:z2*1.e3},{-lat1:lat1})

 f2a = addfile("$DATD/L60CGW-t/epf-wc/uan"+expname+".epfsa_yz_recon0.nc","r")
 f2b = addfile("$DATD/L60CGW-t/epf-wc/uan"+expname+".epf_yz_recon0.nc","r")

 if (iw .eq. 0) then
   wv  = new((/nv-1,nt,nz,ny/),float)
   wv2 = new((/nv-1,nt,nz,ny/),float)
 end if
 ; 0 - 5   : K, R, E-MRG, W-MRG, E-IG, W-IG

 do iv=0, nv-2
   wv2(iv,:,:,:) = f2a->$varname(iv)$(n1:n2,iwv+2, \
      {z1*1.e3:z2*1.e3},{-lat1:lat1})
   if (iwv .eq. 1) then  ; Rossby wave
     wv2(iv,:,:,:) = wv2(iv,:,:,:) + f2a->$varname(iv)$(n1:n2,iwv+12, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1})
   end if
   if (iwv .ge. 4) then  ; IG wave
     wv2(iv,:,:,:) = wv2(iv,:,:,:) + f2b->$varname(iv)$(n1:n2,iwv+4, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1}) \
                                   + f2a->$varname(iv)$(n1:n2,iwv+12, \
        {z1*1.e3:z2*1.e3},{-lat1:lat1})
   end if
 end do
 if ( mod(iwv,2) .eq. 0 ) then
   wv2(1:2,:,:,:) = wv2(1:2,:,:,:)*(-1.)  ; for eastward waves
 end if

 if (iw .eq. 0) then
   wv = wv2
   wv2 = 0.
 end if

end do  ; iw

;--------------------------------------- process 0
 vari  = new((/nv,nt,nz,ny/),float)
 vari2 = new((/nv,nt,nz,ny/),float)
 vari (0,:,:,:) = u
 vari2(0,:,:,:) = u
 vari (1:,:,:,:) = (/ wv /)
 vari2(1:,:,:,:) = (/ wv2 /)

 vari(1,:,:,:) = vari(1,:,:,:)*30.  ; m/s/month
 vari(2:3,:,:,:) = vari(2:3,:,:,:) / 6371.e3 * 1.e3  ; mPa
 vari2(1,:,:,:) = vari2(1,:,:,:)*30.  ; m/s/month
 vari2(2:3,:,:,:) = vari2(2:3,:,:,:) / 6371.e3 * 1.e3  ; mPa

;--------------------------------------- cn level
 cnlev = fspan(-5,5,21)*scl(iwvs(0))
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 21
 he = 33
 xf = 12 + (/0,1,2,3,4,5,0,1,2,3,4,5/)*(wi+6)
; yf = 10 + (/0,1,2,3/)*42 + (/0,0,1,1/)*6
 yf = 10 + (/0,0,0,0,0,0,1,1,1,1,1,1/)*(he+10)

 r_yz = 3.141592/180.*6371.
 vari (2,:,:,:) = vari (2,:,:,:)/r_yz
 vari2(2,:,:,:) = vari2(2,:,:,:)/r_yz
 do k=0, nz-1
   vari (2:3,:,k,:) = vari (2:3,:,k,:)*exp(ht(k)/(2.*6.6e3))
   vari2(2:3,:,k,:) = vari2(2:3,:,k,:)*exp(ht(k)/(2.*6.6e3))
 end do

do iw=0, 0
  iwv = iwvs(iw)

 xstr  = (/"","","","","","","lat [deg]","lat [deg]","lat [deg]","lat [deg]","lat [deg]","lat [deg]"/)
 ystr  = (/"z [km]","","","","","","z [km]","","","","",""/)
 refann = new(nt,logical,False)
 if (iw .eq. 1) then
   refann(nt-1) = True
 end if

  vvv = vari
  if (iw .eq. 1) then
    vvv = vari2
  end if

do n=0, 11 ;nt-1

 contour1 = cnshadeplot(wks,lat,ht/1.e3,vvv(1,n,:,:),cnlev,cnfil,0)
 contour0 = cnlineplot(wks,lat,ht/1.e3,vvv(0,n,:,:),ispan(-50,50,5),0)
 cnline_att(contour0,"D",(/0,0,16/),"T",(/0.6,1.8,0.7/),"C",(/1,colind,1/))

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vvv(2:3,n,::2,::4)
  "vfXArray"    : lat(::4)
  "vfYArray"    : ht(::2)/1.e3
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
;  "vcLineArrowThicknessF" :  vthc
;  "vcLineArrowColor"      :  vcol
  "vcRefMagnitudeF" : scl(iwvs(0))
  "vcRefLengthF"    : 0.05
  "vcGlyphStyle"              : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.12
  "vcFillArrowWidthF"         : 0.1
  "vcPositionMode"            : "ArrowTail"
  "vcRefAnnoOn"          : refann(n)
  "vcRefAnnoArrowAngleF" : 90.
  "vcRefAnnoFontHeightF" : 0.024
  "vcRefAnnoString1On"   : False
;  "vcRefAnnoString1"     : "1 mPa"
;  "vcRefAnnoString2On"   : False
  "vcRefAnnoString2"     : scl(iwvs(0))+" mPa"
  "vcRefAnnoPerimSpaceF" : 0.6
  "vcRefAnnoOrthogonalPosF" : 0.15
;  "vcRefAnnoParallelPosF"   :  1. - nn*0.1
 end create

 set_size(contour1,xf(n),yf(n),wi,he)
 axis_range(contour1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour1,"",xstr(n),ystr(n))

 set_size(contour0,xf(n),yf(n),wi,he)
 axis_range(contour0,-lat1,lat1,z1p,z2p,False,False)
 axis_str(contour0,"","","")

 set_size(vec1,xf(n),yf(n),wi,he)
 axis_range(vec1,-lat1,lat1,z1p,z2p,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)
 cnfinalize(contour0,wks)

 draw(contour1)

 overlay(contour0,vec1)
 draw(contour0)
 NhlRemoveOverlay(contour0,vec1,False)

end do

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
  "lbTitleString" : "[m s:S:-1:N: month:S:-1:N:]"
  "lbTitlePosition" : "Right" ;"Top"
  "lbTitleDirection" : "Across"
 end create
; set_size_lb(label,xf(0),yf(3)+he+10,90,4)
 set_size_lb(label,xf(0),max(yf)+he+10,90,4)

 draw(label)

 frame(wks)

end

