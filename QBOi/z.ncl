load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.0

load "read_dims.ncl"

begin

;--------------------------------------- parameters
; models = (/"ERA-I-nr","MRI","CMAM","MIROC-ESM","ERA-I"/)
 models = (/"ERA-I-nr","MRI","CMAM","ERA-I"/)
; yrs = (/1998,2003/)
 yrs = 2003  ;  96 / 98 / 01 / 03 / 05 / 07
; mmi = "05"
 mmi = "11"
 dd = "02"
 hh = "00"
 latb = 15.
 zp1 = 22. ;24.5  ; 30.2 hPa
 vnheads = (/"u","t"/)
 nv = dimsizes(vnheads)
 vnames = conform_dims((/dimsizes(models),nv/),vnheads+"a",1)
 do iv=0, nv-1
   vnames(:,iv) = where(models .eq. "ERA-I-nr",  \
                        str_upper(vnheads(iv))+"_GDS0_HYBL", vnames(:,iv))
   vnames(:,iv) = where(models .eq. "ERA-I"   ,  \
                        str_upper(vnheads(iv))+"_GDS0_ISBL", vnames(:,iv))
 end do


 wks = cr_wks("z-"+yrs+mmi,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 nm = dimsizes(models)
 nt = dimsizes(yrs)

do im=0, nm-1

;--------------------------------------- read data
 tag = 0
 do iv=0, nv-1
 do iyr=0, nt-1
   yr = yrs(iyr)

   if (models(im) .eq. "ERA-I-nr") then
     file_r = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
          "era-int_f."+vnheads(iv)+".anal."+hh+".ml."+yr+""+mmi+".nc"
     f0 = addfile(file_r,"r")
   else if (models(im) .eq. "ERA-I") then
     file_r = "$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mmi+"/"+  \
          "era-int."+vnheads(iv)+".anal."+hh+".pl."+yr+""+mmi+".nc"
     f0 = addfile(file_r,"r")
   else
     sysout = systemfunc("env_qboi.sh 5 "+vnames(im,iv)+" 6hr XYPT r1i1p1 "+yr+mmi)

     i0 = ind( sysout .eq. "start_"+models(im) ) + 2
     nfile = stringtoint(sysout(i0-1))
     files = sysout(i0:i0+nfile-1)

     f0 = addfile(files(0),"r")
   end if
   end if

   ; time: n1
   tname = f0->$vnames(im,iv)$!0
   t = f0->$tname$
   date = cd_calendar(t,-3)
   n1 = ind( date+"" .eq. yr+mmi+dd+hh )
   if ( ismissing(n1) ) then
     print(" MISSING: "+models(im)+" for "+(yr+mmi+dd+hh))
     exit
   end if

   if (iyr .eq. 0) then
     tag = 1

     lon0 = getlon_q(models(im),f0,vnames(im,iv))
     nx = dimsizes(lon0)

     lat0 = getlat_q(models(im),f0,vnames(im,iv))
     ny = dimsizes(lat0)

     if (models(im) .eq. "ERA-I-nr") then
       lon = lon0(::2)
       lat = lat0(::2)
       nx = dimsizes(lon)
       ny = dimsizes(lat)
     else
       lon = lon0
       lat = lat0
     end if

     if ( iv .gt. 0 ) then
       lev2 = lev
     end if
     lev = getlev_q(models(im),f0,vnames(im,iv))
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

   end if

   if (models(im) .eq. "ERA-I-nr") then
     u0 = f0->$vnames(im,iv)$(n1,k1,::2,::2)
   else
     u0 = f0->$vnames(im,iv)$(n1,k1,:,:)
   end if

   u(iv,iyr,:,:) = u0

 end do  ; iyr
 end do  ; iv

;--------------------------------------- calculation
 u = u - conform_dims(dimsizes(u),dim_avg(u),(/0,1,2/))

 ul = new((/nv,nt,ny,nx/),float)
 ul = 0.
 twopi = 6.2831853
 fc = ezfftf(u)
 do ik=0, 9
   ul = ul + ( conform_dims(dimsizes(ul),fc(0,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(ul),  \
                  cos(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) +  \
               conform_dims(dimsizes(ul),fc(1,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(ul),  \
                  sin(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) )
 end do

;--------------------------------------- cn level
 cnlev = new((/nv,21/),float)
 ltck  = new((/nv,21/),float)
 do iv=0, nv-1
   if (vnheads(iv) .eq. "u") then
     cnlev(iv,:) = tofloat(ispan(-10,10,1))
     ltck (iv,:) = where(mod(cnlev(iv,:),5) .eq. 0, 1.8, 0.8)
   else if (vnheads(iv) .eq. "t") then
     cnlev(iv,:) = tofloat(ispan(-10,10,1))
     ltck (iv,:) = where(mod(cnlev(iv,:),5) .eq. 0, 1.8, 0.8)
     cnlev(iv,:) = cnlev(iv,:)*0.5
   end if
   end if
 end do
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; lcol = (/2,3,4,5,6/)*0+1

;--------------------------------------- contour1
 wi = 120
 he = 30
 xf = 10 + ispan(0,7,1)*(wi+10)*0
 yf = 10 + ispan(0,7,1)*(he+8)

 if (mmi .eq. "11") then
   mmis = "Nov"
 else
   mmis = "May"
 end if
 mstr = where(vnheads(0) .eq. "t", "T", vnheads(0))+"'"
 do iv=1, nv-1
   mstr = mstr+", "+where(vnheads(iv) .eq. "t", "T", vnheads(iv))+"'"
 end do
 mstr = mstr+" at "+zp1+" km, initial time (1 "+mmis+" "+yrs(0)+")"
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 if (im .ne. 0) then
   mstr = ""
 end if
 if (im .ne. nm-1) then
   xstr = ""
 end if
 lbopt = new(nv,string)
 lbopt = ""
 if (im .eq. nm-1) then
   lbopt = where(vnheads .eq. "u", "[m s:S:-1:N:]", lbopt)
   lbopt = where(vnheads .eq. "t", "[K]", lbopt)
 end if
 

;do ipl=0, nt-1
ipl = im

 iv = 1
 contour1 = cnshadeplot(wks,lon,lat,ul(iv,0,:,:),cnlev(iv,:),cnfil,lbopt(iv))
 iv = 0
 cnaddline(contour1,wks,lon,lat,ul(iv,0,:,:),cnlev(iv,:),0,  \
           "D",(/0,0,5/),"T",ltck(iv,:),"C",(/1,-1,1/))

 set_size(contour1,xf(ipl),yf(ipl),wi,he)
 axis_range(contour1,0,360,-latb,latb,False,False)
 axis_str(contour1,mstr,xstr,ystr)
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
  "pmLabelBarWidthF"         : 0.4
  "pmLabelBarHeightF"        : 0.05
;  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.008
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.010
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

;end do  ; ipl

 delete(contour1)

 if (im .eq. 0) then
   txtstr = models
   txtstr = where(models .eq. "ERA-I-nr", "ERA-I (ml)", txtstr)
 end if
 dr_txt2(wks,txtstr(im),xf(ipl)+2,yf(ipl)+he-2.5,22,8,0.,1,1)


 if ( isdefined("sysout") ) then
   delete( [/ sysout, files /] )
 end if
 delete( [/ u, u0, ul /] )
 delete( [/ fc /] )
 delete( [/ t, date /] )
 delete( [/ lon0, lat0, lon, lat, lev, levp, zp /] )
 if ( isdefined("lev2") ) then
   delete( lev2 )
 end if

end do  ; im

 frame(wks)

end

