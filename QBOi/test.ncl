load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 0.8

load "read_dims.ncl"
load "read_var.ncl"
load "hintp.ncl"

begin

;--------------------------------------- parameters
 models = (/"ERA-I-nr","MRI","ERA-I","CMAM","MIROC-ESM","MIROC-AGCM-LL"/)
; models = (/"ERA-I-nr","MRI","CMAM","ERA-I"/)
; yrs = (/1998,2003/)
 yrs = 1994  ;  96 / 98 / 01 / 03 / 05 / 07
; mmi = "05"
 mmi = "11"
 hfcst = 24  ; 6
 latb = 15.
 zp1 = 21.5  ;24.5  ; 30.2 hPa
 vnheads = (/"D","u","v"/)
 wn_min = 1
 wn_max = 8
 nv = dimsizes(vnheads)


; wks = cr_wks("fig/div-v_xy-"+yrs+mmi+"-h"+sprinti("%2.2i",hfcst),600,0)
 wks = cr_wks("test",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 dd = sprinti("%2.2i",1 + hfcst/24)
 hh = sprinti("%2.2i",hfcst - (hfcst/24)*24)

 nm = dimsizes(models)
 nt = dimsizes(yrs)

do im=0, nm-1

;--------------------------------------- read data
 tag = 0
 do iv=0, nv-1
 do iyr=0, nt-1
   yr = yrs(iyr)

   files = getvarfile_q(models(im),vnheads(iv),yr,mmi,hh)
   vname = files@vname
   if ( dimsizes(vname) .eq. 1 ) then
     f0 = addfile(files(0),"r")
   else
     f0 = addfile(files(0,0),"r")
     f1 = addfile(files(1,0),"r")
   end if
   delete( files )

   ; time: n1
   tname = f0->$vname(0)$!0
   t = f0->$tname$
   date = cd_calendar(t,-3)
   n1 = ind( date+"" .eq. yr+mmi+dd+hh )
   if ( ismissing(n1) ) then
     print(" MISSING: "+models(im)+" for "+(yr+mmi+dd+hh))
     exit
   end if

   if (iyr .eq. 0) then
     tag = 1

     lon0 = getlon_q(models(im),f0,vname(0))
     nx0 = dimsizes(lon0)

     lat0 = getlat_q(models(im),f0,vname(0))
     ny0 = dimsizes(lat0)

     if (models(im) .eq. "ERA-I-nr") then
       lon = lon0(::2)
       lat = lat0(0::2)
       nx = dimsizes(lon)
       ny = dimsizes(lat)
     else
       lon = lon0
       lat = lat0
       nx = nx0
       ny = ny0
     end if

     if ( iv .gt. 0 ) then
       lev2 = lev
     end if
     lev = getlev_q(models(im),f0,vname(0))
     if ( iv .gt. 0 ) then
       if ( any(lev .ne. lev2) ) then
         print("Check differences in levels between variables.")
         exit
       end if
     end if

     pscale = 1.
     levp = tofloat( lev )
     if ( isatt(lev, "units") ) then
     if (lev@units .eq. "Pa") then
       pscale = 100.
     end if
     end if
     levp = levp/pscale
     nz = dimsizes(levp)
     zp = 7.*log(1.e3/levp)
     k1 = minind( abs(zp - zp1) )
print(models(im)+"  "+zp(k1)+"  "+tofloat(lat(1)-lat(0)))

     if (iv .eq. 0) then
       u = new((/nv,nt,ny,nx/),float)
     end if

     if (vnheads(iv) .eq. "uw") then
       uw_u0 = new((/nt,ny,nx/),float)
       uw_w0 = new((/nt,ny,nx/),float)
     end if
     if (vnheads(iv) .eq. "D") then
       div_u0 = new((/ny0,nx0/),float)
       div_v0 = new((/ny0,nx0/),float)
       div = new((/ny0,nx0/),float)
     end if

   end if

   if ( all(vnheads(iv) .ne. (/"uw","D"/)) ) then

     if (models(im) .eq. "ERA-I-nr") then
       u0 = f0->$vname(0)$(n1,k1,0::2,::2)
     else
       u0 = f0->$vname(0)$(n1,k1,:,:)
     end if

     if (vnheads(iv) .eq. "w") then
       if ( models(im) .eq. "ERA-I-nr" .or. models(im) .eq. "ERA-I" ) then
         u0 = u0*(-7.e3/(levp(k1)*100.))
       end if
       u0 = u0*1.e3  ; [mm/s]
     end if

     u(iv,iyr,:,:) = u0

   end if

   if (vnheads(iv) .eq. "uw") then

     if (models(im) .eq. "ERA-I-nr") then
       uw_u0(iyr,:,:) = f0->$vname(0)$(n1,k1,0::2,::2)
       uw_w0(iyr,:,:) = f1->$vname(1)$(n1,k1,0::2,::2)
     else
       uw_u0(iyr,:,:) = f0->$vname(0)$(n1,k1,:,:)
       uw_w0(iyr,:,:) = f1->$vname(1)$(n1,k1,:,:)
     end if

     if ( models(im) .eq. "ERA-I-nr" .or. models(im) .eq. "ERA-I" ) then
       uw_w0(iyr,:,:) = uw_w0(iyr,:,:)*(-7.e3/(levp(k1)*100.))
     end if

     u(iv,iyr,:,:) = 0.

   end if

   if (vnheads(iv) .eq. "D") then

     deg2rad = 3.141592/180.
     idlon2 = tofloat(1./((lon0(2) - lon0(0))*deg2rad*6370.e3*cos(lat0*deg2rad)))
     idlat2 = tofloat(1./((lat0(2) - lat0(0))*deg2rad*6370.e3))
     div_u0 = f0->$vname(0)$(n1,k1,:,:)
     div_v0 = f1->$vname(1)$(n1,k1,:,:)*  \
              conform_dims(dimsizes(div_v0),tofloat(cos(lat0*deg2rad)),0)
     div(:,1:nx0-2) = (div_u0(:,2:nx0-1) - div_u0(:,0:nx0-3))*  \
                      conform_dims((/ny0,nx0-2/),idlon2,0)
     div(:,0) = (div_u0(:,1) - div_u0(:,nx0-1))*idlon2
     div(:,nx0-1) = (div_u0(:,0) - div_u0(:,nx0-2))*idlon2
     div(1:ny0-2,:) = div(1:ny0-2,:) +  \
                      (div_v0(2:ny0-1,:) - div_v0(0:ny0-3,:))*idlat2/  \
              conform_dims((/ny0-2,nx0/),tofloat(cos(lat0(1:ny0-2)*deg2rad)),0)
     div(0,:) = div(1,:)
     div(ny0-1,:) = div(ny0-2,:)
 
     div = div*86400.  ; [/day]

     if (models(im) .eq. "ERA-I-nr") then
       u(iv,iyr,:,:) = div(0::2,::2)
     else
       u(iv,iyr,:,:) = div
     end if

   end if

 end do  ; iyr
   delete( vname )
 end do  ; iv

;--------------------------------------- calculation
 ul = xrescale_q(u, lon(0), wn_min, wn_max, 1.5)

 lono = ul&lon

 do iv=0, nv-1
   if (vnheads(iv) .ne. "uw") then
     if (vnheads(iv) .ne. "v") then
       ul(iv,:,:,:) = 0.5*(ul(iv,:,:,:) + ul(iv,:,ny-1:0,:))
     else
       ul(iv,:,:,:) = 0.5*(ul(iv,:,:,:) - ul(iv,:,ny-1:0,:))
     end if
   end if

   if (vnheads(iv) .eq. "uw") then
     uw_u0 = 0.5*(uw_u0 + uw_u0(:,ny-1:0,:))
     uw_w0 = 0.5*(uw_w0 + uw_w0(:,ny-1:0,:))

     uw_ul0 = xrescale_q(uw_u0, lon(0), wn_min, wn_max, 1.5)
     uw_wl0 = xrescale_q(uw_w0, lon(0), wn_min, wn_max, 1.5)

     ul(iv,:,:,:) = uw_ul0*uw_wl0
     ul(iv,:,:,:) = ul(iv,:,:,:)*(levp(k1)/1.e3)
     ul(iv,:,:,:) = ul(iv,:,:,:)*1.e3  ; [mPa]
   end if
 end do

;--------------------------------------- cn level
 cnlev = new((/nv,21/),float)
 ltck  = new((/nv,21/),float)
 do iv=0, nv-1
   cnlev(iv,:) = tofloat(ispan(-10,10,1))
   ltck (iv,:) = where(mod(cnlev(iv,:),5) .eq. 0, 1.8, 0.8)
   if (vnheads(iv) .eq. "t") then
     cnlev(iv,:) = cnlev(iv,:)*0.5
   else if (vnheads(iv) .eq. "w") then
     cnlev(iv,:) = cnlev(iv,:)*1.5
   else if (vnheads(iv) .eq. "uw") then
;     cnlev(iv,:) = cnlev(iv,:)*0.1
     cnlev(iv,:) = cnlev(iv,:)*0.2
   else if (vnheads(iv) .eq. "D") then
     cnlev(iv,:) = cnlev(iv,:)*0.1
   end if
   end if
   end if
   end if
 end do
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; lcol = (/2,3,4,5,6/)*0+1

;--------------------------------------- contour1
 wi = 60
 he = 15
 xf = 10 + (/0,0,1,1,1,1/)*(wi+10)
 yf = 15 + (/0,1,0,1,2,3/)*(he+7)
 nm_x0 = 2

 if (mmi .eq. "11") then
   mmis = "Nov"
 else
   mmis = "May"
 end if
 mstr = where(vnheads(0) .eq. "t", "T", vnheads(0))+"'"
 do iv=1, nv-1
   mstr = mstr+", "+where(vnheads(iv) .eq. "t", "T", vnheads(iv))+"'"
 end do
 mstr = mstr+" at "+zp1+" km, +"+sprinti("%2.2i",hfcst)+" h ("+mmis+" "+yrs(0)+")"
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 if ( .not. ( im .eq. nm-1 .or. im .eq. nm_x0-1 ) ) then
   xstr = ""
 end if
 lbopt = new(nv,string)
 lbopt = ""
 if (im .eq. nm-1) then
   lbopt = where(vnheads .eq. "u", "[m s:S:-1:N:]", lbopt)
   lbopt = where(vnheads .eq. "t", "[K]", lbopt)
   lbopt = where(vnheads .eq. "w", "[mm s:S:-1:N:]", lbopt)
   lbopt = where(vnheads .eq. "uw", "[mPa]", lbopt)
   lbopt = where(vnheads .eq. "D", "[day:S:-1:N:]", lbopt)
 end if
 

;do ipl=0, nt-1
ipl = im

 iv = 0
 contour1 = cnshadeplot(wks,lono,lat,ul(iv,0,:,:),cnlev(iv,:),cnfil,lbopt(iv))
; iv = 1
; cnaddline(contour1,wks,lono,lat,ul(iv,0,:,:),cnlev(iv,:),0,  \
;           "D",(/0,0,5/),"T",ltck(iv,:),"C",(/1,-1,1/))

 set_size(contour1,xf(ipl),yf(ipl),wi,he)
 axis_range(contour1,0,360,-latb,latb,False,False)
 axis_str(contour1,"",xstr,ystr)
 setvalues contour1
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(0,360,90)
  "tmXBLabels" : ispan(0,360,90)
  "tmXBMinorValues" : ispan(0,360,15)
 end setvalues
; if (im .ne. 0) then
;   setvalues contour1
;    "tiMainOn"  : False
;    "tiXAxisOn" : False
;    "tiYAxisOn" : False
;   end setvalues
; end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.3
  "pmLabelBarHeightF"        : 0.04
;  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.006
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.008
 end setvalues

 cnfinalize(contour1,wks)

 itvvec = 10
 uvvec = yintp_q(ul(1:2,0,:,::itvvec), lat, 3., 15.)
 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : uvvec
  "vfXArray" : lono(::itvvec)
  "vfYArray" : uvvec&lat
 end create


 ref = 10. ;/(ir+1)
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(ipl),yf(ipl),wi,he)
 axis_range(vec1,0,360,-latb,latb,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)

 draw(contour1)

;end do  ; ipl

 delete(contour1)

 if (im .eq. 0) then
   txtstr = models
   txtstr = where(models .eq. "ERA-I-nr", "ERA-I (ml)", txtstr)

   dr_txt2(wks,mstr,xf(0)+wi*0.6,yf(0)-8,22,6.0,0.,1,1)

 end if
 dr_txt2(wks,txtstr(im),xf(ipl)+1,yf(ipl)-2,22,5.5,0.,1,1)


 delete( [/ u, u0, ul /] )
 delete( [/ t, date /] )
 delete( [/ lon0, lat0, lon, lat, lev, levp, zp /] )
 if ( isdefined("lev2") ) then
   delete( lev2 )
 end if
 if ( any( vnheads .eq. "uw" ) ) then
   delete( [/ uw_u0, uw_w0, fc_u, fc_w /] )
 end if
 if ( any( vnheads .eq. "D" ) ) then
   delete( [/ div_u0, div_v0, div, idlon2, idlat2 /] )
 end if

end do  ; im

 frame(wks)

end

