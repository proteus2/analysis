load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.0

load "read_dims.ncl"

begin

;--------------------------------------- parameters
 models = (/"ERA-I-nr","MRI","CMAM","MIROC-ESM","ERA-I"/)
; models = (/"CMAM","MRI"/)
; yrs = (/1998,2003/)
 yrs = 1998
; mmi = "05"
 mmi = "11"
 dd = "01"
 hh = "06"
 latb = 15.
 zp1 = 24.5  ; 30.2 hPa
; vnhead = "u"
 vnhead = "t"
 vnames = conform_dims(dimsizes(models),vnhead+"a",-1)
 vnames = where(models .eq. "ERA-I-nr", str_upper(vnhead)+"_GDS0_HYBL", vnames)
 vnames = where(models .eq. "ERA-I"   , str_upper(vnhead)+"_GDS0_ISBL", vnames)


 wks = cr_wks("fig/xy_init-"+yrs+mmi,600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 nm = dimsizes(models)
 nt = dimsizes(yrs)

do im=0, nm-1

;--------------------------------------- read data
 tag = 0
 do iyr=0, nt-1
   yr = yrs(iyr)

   if (models(im) .eq. "ERA-I-nr") then
     file_r = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
          "era-int_f."+vnhead+".anal."+hh+".ml."+yr+""+mmi+".nc"
     f0 = addfile(file_r,"r")
   else if (models(im) .eq. "ERA-I") then
     file_r = "$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mmi+"/"+  \
          "era-int."+vnhead+".anal."+hh+".pl."+yr+""+mmi+".nc"
     f0 = addfile(file_r,"r")
   else
     sysout = systemfunc("env_qboi.sh 5 "+vnames(im)+" 6hr XYPT r1i1p1 "+yr+mmi)

     i0 = ind( sysout .eq. "start_"+models(im) ) + 2
     nfile = stringtoint(sysout(i0-1))
     files = sysout(i0:i0+nfile-1)

     f0 = addfile(files(0),"r")
   end if
   end if

   ; time: n1
   tname = f0->$vnames(im)$!0
   t = f0->$tname$(0:4)
   date = cd_calendar(t,-3)
   n1 = ind( date+"" .eq. yr+mmi+dd+hh )
   if ( ismissing(n1) ) then
     print(" MISSING: "+models(im)+" for "+(yr+mmi+dd+hh))
     exit
   end if

   if (iyr .eq. 0) then
     tag = 1

     lon = getlon_q(models(im),f0,vnames(im))
     nx = dimsizes(lon)

     lat = getlat_q(models(im),f0,vnames(im))
     ny = dimsizes(lat)

     if (models(im) .eq. "ERA-I-nr") then
       tmp = lon
       delete(lon)
       lon = tmp(::2)
       delete(tmp)
       tmp = lat
       delete(lat)
       lat = tmp(::2)
       delete(tmp)
       nx = dimsizes(lon)
       ny = dimsizes(lat)
     end if

     lev = getlev_q(models(im),f0,vnames(im))
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

     u = new((/nt,ny,nx/),float)

   end if

   if (models(im) .eq. "ERA-I-nr") then
     u0 = f0->$vnames(im)$(n1,k1,::2,::2)
   else
     u0 = f0->$vnames(im)$(n1,k1,:,:)
   end if

   u(iyr,:,:) = u0

 end do  ; iyr

;--------------------------------------- calculation
 u = u - conform_dims(dimsizes(u),dim_avg(u),(/0,1/))

 ul = new((/nt,ny,nx/),float)
 ul = 0.
 twopi = 6.2831853
 fc = ezfftf(u)
 do ik=0, 9
   ul = ul + ( conform_dims(dimsizes(ul),fc(0,:,:,ik),(/0,1/))*  \
               conform_dims(dimsizes(ul),  \
                  cos(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),2) +  \
               conform_dims(dimsizes(ul),fc(1,:,:,ik),(/0,1/))*  \
               conform_dims(dimsizes(ul),  \
                  sin(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),2) )
 end do

;--------------------------------------- cn level
 if (vnhead .eq. "u") then
   cnlev = ispan(-10,10,1)
   ltck = where(mod(cnlev,5) .eq. 0, 1.6, 0.8)
 else if (vnhead .eq. "t") then
   cnlev = ispan(-10,10,1)*1.
   ltck = where(mod(cnlev,5) .eq. 0, 1.6, 0.8)
   cnlev = cnlev*0.5
 end if
 end if

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
 mstr = " at "+zp1+" km, initial time (1 "+mmis+" "+yrs(0)+") "
 if (vnhead .eq. "u") then
   mstr = "u' "+mstr+" [m s:S:-1:N:]"
 else if (vnhead .eq. "t") then
   mstr = "T' "+mstr+" [K]"
 end if
 end if
 xstr = "lon [deg]"
 ystr = "lat [deg]"
 if (im .ne. 0) then
   mstr = ""
 end if
 if (im .ne. nm-1) then
   xstr = ""
 end if
 

;do ipl=0, nt-1
ipl = im

 contour1 = cnlineplot(wks,lon,lat,ul(0,:,:),cnlev,0)
 cnline_att(contour1,"D",(/0,0,5/),"T",ltck,"C",(/1,0,1/))
; setvalues contour1
;  "cnLineLabelFontHeightF" : 0.022
; end setvalues

 set_size(contour1,xf(ipl),yf(ipl),wi,he)
 axis_range(contour1,0,360,-latb,latb,False,False)
 axis_str(contour1,mstr,xstr,ystr)
; if (im .ne. 0) then
;   setvalues contour1
;    "tiMainOn"  : False
;    "tiXAxisOn" : False
;    "tiYAxisOn" : False
;   end setvalues
; end if

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
 delete( t )
 delete( [/ lon, lat, lev, levp, zp /] )

end do  ; im

 frame(wks)

end

