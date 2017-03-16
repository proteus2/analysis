load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 0.9

load "read_dims.ncl"

begin

;--------------------------------------- parameters
 models = (/"ERA-I-nr","CMAM","MIROC-ESM","MRI"/)
; models = (/"CMAM","MRI"/)
 y1 = 1993
 y2 = 2007
; mmi = "05"
 mmi = "11"
 dd = "01"
 hh = "06"
 lata = 5.
 vnames = conform_dims(dimsizes(models),"ua",-1)
 vnames = where(models .eq. "ERA-I-nr", "U_GDS0_HYBL", vnames)
 vnames = where(models .eq. "ERA-I"   , "U_GDS0_ISBL", vnames)
 

 wks = cr_wks("fig/u_init1-"+mmi,600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 nm = dimsizes(models)
 nt = y2 - y1 + 1

do im=0, nm-1

;--------------------------------------- read data
 tag = 0
 do iyr=0, nt-1
   yr = y1 + iyr

   if (models(im) .eq. "ERA-I-nr") then
     file_r = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
          "era-int_f.u.anal."+hh+".ml."+yr+""+mmi+".nc"
     f0 = addfile(file_r,"r")
   else
     sysout = systemfunc("env_qboi.sh 5 "+vnames(im)+" 6hr XYPT r1i1p1 "+yr+mmi)

     i0 = ind( sysout .eq. "start_"+models(im) ) + 2
     nfile = stringtoint(sysout(i0-1))
     files = sysout(i0:i0+nfile-1)

     f0 = addfile(files(0),"r")
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

     lat = getlat_q(models(im),f0,vnames(im))
     lat_bnd1 = lat@bnd1
     lat_bnd2 = lat@bnd2
     dlat = abs(lat_bnd2 - lat_bnd1)

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

     if (lat(0) .lt. lat(1)) then
       j1 = min( ind( lat .ge. -lata ) )
       j2 = max( ind( lat .le. lata  ) )
     else
       j1 = min( ind( lat .le. lata  ) )
       j2 = max( ind( lat .ge. -lata ) )
     end if
     dlat(j1) = dlat(j1)*0.5
     dlat(j2) = dlat(j2)*0.5
     sumdlat = sum(dlat(j1:j2))

     dlat1b = lata - abs(lat(j1))
     dlat2b = lata - abs(lat(j2))

     if ( dlat1b .ne. 0. .or. dlat2b .ne. 0. ) then
       tag = 2
       dlat1o = lat(j1) - lat(j1-1)
       dlat2o = lat(j2+1) - lat(j2)
     end if

     u = new((/nt,nz/),float)

   end if

   u0 = dim_avg_Wrap( f0->$vnames(im)$(n1,:,:,:) )

   u(iyr,:) = tofloat( dim_sum(  \
       u0(:,j1:j2)*conform_dims(dimsizes(u0(:,j1:j2)),dlat(j1:j2),1) ) )

   if (tag .eq. 2) then
     if (iyr .eq. 0) then
;       sumdlat = sumdlat + (dlat1b + dlat2b)
       sumdlat = lata*2.
     end if

     ub1 = ( (0.5*dlat1b)*u0(:,j1-1) + (dlat1o - 0.5*dlat1b)*u0(:,j1) )/dlat1o
     ub2 = ( (0.5*dlat2b)*u0(:,j2+1) + (dlat2o - 0.5*dlat2b)*u0(:,j2) )/dlat2o

     u(iyr,:) = u(iyr,:) + tofloat(ub1*dlat1b + ub2*dlat2b)
   end if

   u(iyr,:) = u(iyr,:)/tofloat(sumdlat)

 end do  ; iyr

;print(lat_bnd1+"   "+lat_bnd2)

;--------------------------------------- cn level
 lcol = (/1,2,3,4,5,6/)
 ltck = 0.8

;--------------------------------------- plot
 wi = 24
 he = 24
 xf = 10 + (/0,1,2,3,4,0,1,2,3,4,0,1,2,3,4/)*(wi+5)
 yf = 10 + (/0,0,0,0,0,1,1,1,1,1,2,2,2,2,2/)*(he+5)

do iyr=0, nt-1
 yr = y1 + iyr

 if (mmi .eq. "11") then
   mmis = "Nov"
 else
   mmis = "May"
 end if
 mstr = "1 "+mmis+" "+yr
 xstr = "U [m s:S:-1:N:]"
 ystr = "z [km]"
 if (xf(iyr) .ne. xf(0)) then
   ystr = ""
 end if
 if (yf(iyr) .ne. yf(nt-1)) then
   xstr = ""
 end if
 

 plot1 = xylineplot(wks,u(iyr,:),zp,"D",0,"T",ltck,"C",lcol(im))
 xyaddmarker(plot1,1,"S",0.8,"T",0.6,"C",1)

 set_size(plot1,xf(iyr),yf(iyr),wi,he)
 axis_range(plot1,-45,25,16,42,False,False)
 axis_str(plot1,"",xstr,ystr)
 dr_txt(wks,mstr,xf(iyr)+7,yf(iyr)+he-2,21,5,0.) 
 if (im .ne. 0) then
   setvalues plot1
    "tiMainOn"  : False
    "tiXAxisOn" : False
    "tiYAxisOn" : False
   end setvalues
 end if

 draw(plot1)

end do

 delete(plot1)

 if (im .eq. 0) then
   lg = lgplot(wks,models,"D",0,"T",ltck,"C",lcol(:nm-1),2)
   set_size_lg(lg,xf(0),yf(nt-1)+he+10,wi*0.5,3.5*nm)
   draw(lg)
 end if


 delete( [/ sysout, files /] )
 delete( [/ u, u0, ub1, ub2 /] )
 delete( t )
 delete( [/ lat, lev, levp, zp /] )
 delete( [/ lat_bnd1, lat_bnd2, dlat /] )

end do  ; im

 frame(wks)

end

