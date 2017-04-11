load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.0

load "read_dims.ncl"
load "lat_intg.ncl"

begin

;--------------------------------------- parameters
 models = (/"ERA-I-nr"/)
; models = (/"ERA-I-nr","MRI","CMAM"/)  ; ,"MIROC-ESM"/)
; models = (/"ERA-I-nr","MRI","CMAM","ERA-I"/)
 yrs = ispan(1993,2007,1)
 mmsi = (/"05","11"/)
 dd = "01"
 hh = "00"
 lata = 10.
 vnheads = (/"u","uw"/)
 nv = dimsizes(vnheads)
 vnames = conform_dims((/dimsizes(models),nv/),vnheads+"a",1)
 do iv=0, nv-1
   vnames(:,iv) = where(models .eq. "ERA-I-nr",  \
                        str_upper(vnheads(iv))+"_GDS0_HYBL", vnames(:,iv))
   vnames(:,iv) = where(models .eq. "ERA-I"   ,  \
                        str_upper(vnheads(iv))+"_GDS0_ISBL", vnames(:,iv))
 end do
 vfname_ra = where(vnheads .eq. "w", "ome", vnheads)
 

 wks = cr_wks("fig/crit1",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

 nm = dimsizes(models)
 ne = dimsizes(yrs)*dimsizes(mmsi)
 nt = 2

 cri2d = new((/dimsizes(yrs),dimsizes(mmsi)/),logical)
 cri2d = False

do im=0, nm-1

ie = 0
do iyr=0, dimsizes(yrs)-1
do imn=0, dimsizes(mmsi)-1
 yr = yrs(iyr)
 mmi = mmsi(imn)

if (ie .eq. 16) then
 frame(wks)
 delete(wks)

 wks = cr_wks("fig/crit2",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)
end if

;--------------------------------------- read data
 tag = 0
 do iv=0, nv-1

   if (vnheads(iv) .ne. "uw") then

     if (models(im) .eq. "ERA-I-nr") then
       file_r = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
            "era-int_f."+vfname_ra(iv)+".anal."+hh+".ml."+yr+""+mmi+".nc"
       f0 = addfile(file_r,"r")
     else if (models(im) .eq. "ERA-I") then
       file_r = "$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mmi+"/"+  \
            "era-int."+vfname_ra(iv)+".anal."+hh+".pl."+yr+""+mmi+".nc"
       f0 = addfile(file_r,"r")
     else
       sysout = systemfunc("env_qboi.sh 5 "+vnames(im,iv)+" 6hr XYPT r1i1p1 "+yr+mmi)

       i0 = ind( sysout .eq. "start_"+models(im) ) + 2
       nfile = stringtoint(sysout(i0-1))
       files = sysout(i0:i0+nfile-1)

       f0 = addfile(files(0),"r")
     end if
     end if

   end if

   if (vnheads(iv) .eq. "uw") then
     if (models(im) .eq. "ERA-I-nr") then
       file_ur = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
            "era-int_f.u.anal."+hh+".ml."+yr+""+mmi+".nc"
       file_wr = "$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mmi+"/"+  \
            "era-int_f.ome.anal."+hh+".ml."+yr+""+mmi+".nc"
       f_u = addfile(file_ur,"r")
       f_w = addfile(file_wr,"r")
     else if (models(im) .eq. "ERA-I") then
       file_ur = "$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mmi+"/"+  \
            "era-int.u.anal."+hh+".pl."+yr+""+mmi+".nc"
       file_wr = "$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mmi+"/"+  \
            "era-int.ome.anal."+hh+".pl."+yr+""+mmi+".nc"
       f_u = addfile(file_ur,"r")
       f_w = addfile(file_wr,"r")
     else
       sysout_u = systemfunc("env_qboi.sh 5 ua 6hr XYPT r1i1p1 "+yr+mmi)
       sysout_w = systemfunc("env_qboi.sh 5 wa 6hr XYPT r1i1p1 "+yr+mmi)

       i0_u = ind( sysout_u .eq. "start_"+models(im) ) + 2
       nfile_u = stringtoint(sysout_u(i0_u-1))
       files_u = sysout_u(i0_u:i0_u+nfile_u-1)
       i0_w = ind( sysout_w .eq. "start_"+models(im) ) + 2
       nfile_w = stringtoint(sysout_w(i0_w-1))
       files_w = sysout_w(i0_w:i0_w+nfile_w-1)

       f_u = addfile(files_u(0),"r")
       f_w = addfile(files_w(0),"r")
     end if
     end if
     f0 = f_u
     vnames(im,iv) = "ua"
     vname_w = "wa"
     if (models(im) .eq. "ERA-I-nr") then
       vnames(im,iv) = "U_GDS0_HYBL"
       vname_w = "W_GDS0_HYBL"
     else if (models(im) .eq. "ERA-I") then
       vnames(im,iv) = "U_GDS0_ISBL"
       vname_w = "W_GDS0_ISBL"
     end if
     end if
   end if

   ; time: n1
   tname = f0->$vnames(im,iv)$!0
   if ( isdefined("t") ) then
     delete( [/t, date/] )
   end if
   t = f0->$tname$
   date = cd_calendar(t,-3)
   n1 = ind( date+"" .eq. yr+mmi+dd+hh )
   if ( ismissing(n1) ) then
     print(" MISSING: "+models(im)+" for "+(yr+mmi+dd+hh))
     exit
   end if

   if (ie .eq. 0) then
     tag = 1

     lon0 = getlon_q(models(im),f0,vnames(im,iv))
     nx = dimsizes(lon0)

     lat0 = getlat_q(models(im),f0,vnames(im,iv))
     lat_bnd1 = lat0@bnd1
     lat_bnd2 = lat0@bnd2
     ny = dimsizes(lat0)

     lon = lon0
     lat = lat0
     deg2rad = 0.01745329
     coslat = cos(lat*deg2rad)

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
print(models(im)+"  "+tofloat(lat(1)-lat(0)))

     if (iv .eq. 0) then
       u = new((/nv,nt,nz,nx/),float)
     end if

     if (vnheads(iv) .eq. "uw") then
       uw_u0 = new((/nt,nz,ny,nx/),float)
       uw_w0 = new((/nt,nz,ny,nx/),float)
     end if

   end if

   if (vnheads(iv) .ne. "uw") then

     u0 = f0->$vnames(im,iv)$(n1:n1+nt-1,:,:,:)

     if (vnheads(iv) .eq. "w") then
       if ( models(im) .eq. "ERA-I-nr" .or. models(im) .eq. "ERA-I" ) then
         u0 = u0*conform_dims(dimsizes(u0),(-7.e3/(levp*100.)),1)
       end if
       u0 = u0*1.e3  ; [mm/s]
     end if

     u(iv,:,:,:) = lat_intg(lat,lat_bnd1,lat_bnd2,(/-lata,lata/),u0,2,  \
                   coslat)

     sumdlat = u@sumdlat

     u(iv,:,:,:) = u(iv,:,:,:)/sumdlat

   end if

   if (vnheads(iv) .eq. "uw") then

     uw_u0 = f0->$vnames(im,iv)$(n1:n1+nt-1,:,:,:)
     uw_w0 = f_w->$vname_w$(n1:n1+nt-1,:,:,:)

     if ( models(im) .eq. "ERA-I-nr" .or. models(im) .eq. "ERA-I" ) then
       uw_w0 = uw_w0*conform_dims(dimsizes(uw_w0),(-7.e3/(levp*100.)),1)
     end if

     u(iv,:,:,:) = 0.

   end if

 end do  ; iv

;--------------------------------------- calculation
 u = u - conform_dims(dimsizes(u),dim_avg(u),(/0,1,2/))

 nw = 10

 ul = new((/nv,nt,nz,nx/),float)
 ul = 0.
 twopi = 6.2831853
 fc = ezfftf(u)
 do ik=0, nw-1
   ul = ul + ( conform_dims(dimsizes(ul),fc(0,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(ul),  \
                  cos(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) +  \
               conform_dims(dimsizes(ul),fc(1,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(ul),  \
                  sin(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) )
 end do

 do iv=0, nv-1
   if (vnheads(iv) .eq. "uw") then
     uw_u0 = uw_u0 - conform_dims(dimsizes(uw_u0),dim_avg(uw_u0),(/0,1,2/))
     uw_w0 = uw_w0 - conform_dims(dimsizes(uw_w0),dim_avg(uw_w0),(/0,1,2/))

     uw_u0 = 0.5*(uw_u0 + uw_u0(:,:,ny-1:0,:))
     uw_w0 = 0.5*(uw_w0 + uw_w0(:,:,ny-1:0,:))

     uwl0 = new((/nt,nz,ny,nx/),float)
     uwl0 = 0.
     fc_u = ezfftf(uw_u0)
     fc_w = ezfftf(uw_w0)
     do ik=0, nw-1
       uwl0 = uwl0 + (  \
             ( conform_dims(dimsizes(uwl0),fc_u(0,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(uwl0),  \
                  cos(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) +  \
               conform_dims(dimsizes(uwl0),fc_u(1,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(uwl0),  \
                  sin(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) )*   \
             ( conform_dims(dimsizes(uwl0),fc_w(0,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(uwl0),  \
                  cos(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) +  \
               conform_dims(dimsizes(uwl0),fc_w(1,:,:,:,ik),(/0,1,2/))*  \
               conform_dims(dimsizes(uwl0),  \
                  sin(ispan(0,nx-1,1)*(ik+1)/tofloat(nx)*twopi),3) ) )
     end do

     uwl = lat_intg(lat,lat_bnd1,lat_bnd2,(/-lata,lata/),uwl0,2,  \
                    coslat)

     sumdlat_uw = uwl@sumdlat

     uwl = uwl/sumdlat_uw

     uwl = uwl*conform_dims(dimsizes(uwl),exp(-zp/7.),1)

     ul(iv,:,:,:) = uwl*1.e3  ; [mPa]
   end if
 end do


 if ( True ) then
   cri2d(iyr,imn) = True
 end if

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
     cnlev(iv,:) = cnlev(iv,:)*0.1
   end if
   end if
   end if
 end do
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

; lcol = (/2,3,4,5,6/)*0+1

;--------------------------------------- contour1
 wi = 60
 he = 20
 xf = 10 + (/0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1/)*(wi+10)
 yf = 10 + (/0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7/)*(he+8)

ipl = ie
if (ie .ge. 16) then
  ipl = ie - 16
end if

 if (mmi .eq. "11") then
   mmis = "Nov"
 else
   mmis = "May"
 end if
 mstr = where(vnheads(0) .eq. "t", "T", vnheads(0))+"'"
 do iv=1, nv-1
   mstr = mstr+", "+where(vnheads(iv) .eq. "t", "T", vnheads(iv))+"'"
 end do
 mstr = mstr+" at "+hh+" UTC "+dd+" "+mmis+" "+yr
 xstr = "lon [deg]"
 ystr = "z [km]"
; if (im .ne. 0) then
;   mstr = ""
; end if
 if (mod(ipl+1,8) .ne. 0) then
   xstr = ""
 end if
 lbopt = new(nv,string)
 lbopt = ""
 if (ie .eq. ne-1) then
   lbopt = where(vnheads .eq. "u", "[m s:S:-1:N:]", lbopt)
   lbopt = where(vnheads .eq. "t", "[K]", lbopt)
   lbopt = where(vnheads .eq. "w", "[mm s:S:-1:N:]", lbopt)
   lbopt = where(vnheads .eq. "uw", "[mPa]", lbopt)
 end if
 

 iv = 1
 contour1 = cnshadeplot(wks,lon,zp,ul(iv,0,:,:),cnlev(iv,:),cnfil,lbopt(iv))
 iv = 0
 cnaddline(contour1,wks,lon,zp,ul(iv,0,:,:),cnlev(iv,:),0,  \
           "D",(/0,0,5/),"T",ltck(iv,:),"C",(/1,-1,1/))

 set_size(contour1,xf(ipl),yf(ipl),wi,he)
 axis_range(contour1,0,360,16,34,False,False)
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
  "pmLabelBarWidthF"         : 0.2
  "pmLabelBarHeightF"        : 0.03
;  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.008
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.010
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

; if (im .eq. 0) then
;   txtstr = models
;   txtstr = where(models .eq. "ERA-I-nr", "ERA-I (ml)", txtstr)
; end if
; dr_txt2(wks,txtstr(im),xf(ipl)+2,yf(ipl)+he-2.5,22,8,0.,1,1)


 ie = ie + 1
end do  ; imn
end do  ; iyr
 

 if ( isdefined("sysout") ) then
   delete( [/ sysout, files /] )
 end if
 if ( isdefined("sysout_u") ) then
   delete( [/ sysout_u, files_u /] )
 end if
 delete( [/ u, u0, ul /] )
 delete( [/ fc /] )
 delete( [/ t, date /] )
 delete( [/ lon0, lat0, lon, lat, lev, levp, zp /] )
 delete( [/ lat_bnd1, lat_bnd2 /] )
 delete( coslat )
 if ( isdefined("lev2") ) then
   delete( lev2 )
 end if
 if ( any( vnheads .eq. "uw" ) ) then
   delete( [/ uw_u0, uw_w0, uwl0, uwl, fc_u, fc_w /] )
 end if

end do  ; im

 frame(wks)

end

