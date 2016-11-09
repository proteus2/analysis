load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"ERA-I"/)
; expname = (/"era-int"/)
; varname = (/"U_GDS0_ISBL"/)
 expname = (/"era-int_f"/)
 varname = (/"u"/)
 ppid    = (/"j" /)

 y0 = 1953
 y1 = 1980+8
 y2 = 1991+8
 lat1 = 1.5

 wks = cr_wks("fig/fu-req1_tz0-eraf1.5-obs",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = 336 - 96
 t = 1987. + (fspan(1,nt,nt)-0.5)/12.
 za1 = 15
 if (expname .eq. "era-int_f") then
   za1 = 2
 end if

 fnam0 = "/data18/kyh/qbo_obs/qbo_sonde_hr.nc"
 fnam1 = "$DATD/clim/"+cexp(0)+"/"+expname(0)+".u_ypt.1979-2006.nc"
 fnam2 = "$DATD/tem/"+cexp(0)+"/"+expname(0)+".tem_ypt.1979-2006.nc"
 f0 = addfile(fnam0,"r")
 f1 = addfile(fnam1,"r")
 f2 = addfile(fnam2,"r")

 zo = 18.6 + 6.6*log( 70./ f0->p )
 if (expname .eq. "era-int_f") then
   fp = addfile("/data7/ERA-I-nr/ERA-Interim_coordvars.nc","r")
   p00 = (fp->a_model_alt)/100. + (fp->b_model_alt)*1.e3
   zu = 18.6 + 6.6*log( 70./ p00( f1->$f1->$varname(0)$!1$ - 1 ) )
 else
   zu = 18.6 + 6.6*log( 70./ f1->$f1->$varname(0)$!1$ )
 end if
 za = 18.6 + 6.6*log( 70./ f2->p )

 uo = f0->u(:,:nt-1)
 u0 = transpose( ( f1->$varname(0)$(96:335,:,{lat1}) ) )

 ut   = transpose( ( f2->u_tend(96:335,:,{lat1}) ) )
; advz = transpose( ( f2->uadv_z(96:335,:,{lat1}) ) )
 advy = transpose( ( f2->uadv_y(96:335,:,{lat1}) ) )
 advy = advy + transpose( ( f2->cor(96:335,:,{lat1}) ) )

 ws = transpose( ( f2->w_res(96:335,:,{lat1}) ) )

;--------------------------------------- process 0
 nzo = dimsizes(zo)
 nzu = dimsizes(zu)
 nza = dimsizes(za)

 us = u0*0.
 uso = uo*0.
 do k=1, nzu-2
   us(k,:) = (u0(k+1,:) - u0(k-1,:))/(zu(k+1) - zu(k-1))/1.e3
 end do
 do k=1, nzo-2
   uso(k,:) = (uo(k+1,:) - uo(k-1,:))/(zo(k+1) - zo(k-1))/1.e3
 end do
 uso(nzo-1,:) = uso(nzo-2,:)

 uto = uo*0.
 do n=1, nt-2
   uto(:,n) = (uo(:,n+1) - uo(:,n-1))/60.
 end do

 wsi = uo*0.
 advyi = uo*0.
 uti = uo*0.
 do n=0, nt-1
   wsi  (:,n) = ftcurv(za(za1:),ws  (za1:,n),zo)
   advyi(:,n) = ftcurv(za(za1:),advy(za1:,n),zo)
   uti  (:,n) = ftcurv(za(za1:),ut  (za1:,n),zo)
 end do

 advz = -wsi*uso*86400.

; tend = (/uti,advz,uti-advz-advyi/)
 tend = (/uto,advz,uto-advz-advyi/)

 tend = tend*30.

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*2.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 3

 wi = 64 + (64./13.*2.)
 he = 15
 xf = 17 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,7,1)*(he+5)

 mstr  = new((/dimsizes(cexp),nl/),string)
 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:,:) = ""
 mstr(:,0) = "OBS-ERA combined"
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"z (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz","REQ"/)

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,zo,tend(il,:,:),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues

 set_size(contour1,xf(0),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,mstr(0,il),xstr(il),ystr(0))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

; contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 contouru = cnlineplot(wks,t,zo,uo,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(0),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,19,33,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(il),xf(0)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)

 delete(contour1)
 delete(contouru)

end do


 frame(wks)

end



