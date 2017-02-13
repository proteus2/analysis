load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin
;--------------------------------------- parameters
 cexp    = "CGW"
 expname = "tc"
 ppid    = "j"
 varname = (/"u","theta"/)

 y1 = 1955
 y2 = y1 - 1 + 40
 ydif = 100
 mons = 13   ; 1--12 : each month  13 : annual mean  ; 20--23 : djf--son
 lat0 = 15.

 wks = cr_wks("fig/bg_xz-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_merge_colormaps(wks,"StepSeq25","cyclic")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

;--------------------------------------- read data
 nv = dimsizes(varname)

 fdir  = "/hippo0/HG2CMIP/L60"+cexp+"/pj"
 fhead = "uan"+expname+"a.p"+ppid

 nmon = 1
 if (mons .eq. 13) then
   nmon = 12
 end if
 if ( mons .ge. 20 .and. mons .le. 23 ) then
   nmon = 3
 end if
 nt = (y2-y1+1)*nmon

 f0 = addfile(fdir+"/"+fhead+"_"+y1+"12.nc","r")
 tmp = f0->$varname(0)$(0,:,0,:)
 z = f0->$tmp!0$({:25.e3}) / 1.e3
 x = f0->$tmp!1$
 nz = dimsizes(z)
 nx = dimsizes(x)
 delete(tmp)

 ft0 = addfile("/hippo0/HG2CMIP/L60CGW/invariant/L60_z.nc","r")
 ; z at theta lev
 dzr3 = 1./(ft0->ht(0,1:nz,{-lat0:lat0},:) - ft0->ht(0,0:nz-1,{-lat0:lat0},:))

 var = new((/nv,3,nz,nx/),float)
 var = 0.
 do iv=0, nv-1
 do ip=0, 1
 do yr=y1+(ip*ydif), y2+(ip*ydif)
   mon = mons
   do imon=0, nmon-1
     if (mons .eq. 13) then
       mon = imon+1
     else if ( mons .ge. 20 .and. mons .le. 23 ) then
       mon = imon+(mons-20)*3
       if (mon .lt. 1) then
         mon = mon+12  ; not continuous
       end if
     end if
     end if

     f0 = addfile(fdir+"/"+fhead+"_"+yr+sprinti("%2.2i",mon)+".nc","r")
     var(iv,ip,:,:) = var(iv,ip,:,:) + dim_avg_n(  \
        f0->$varname(iv)$(0,:nz-1,{-lat0:lat0},:), 1 )/tofloat(nt)

   end do 
 end do 
 end do 
 end do 

;--------------------------------------- process 0
 tmp = var(1,:,:,:)
 ; if theta on w lev.
 do k=1, nz-1
   var(1,0:1,k,:) = 9.8*log(tmp(0:1,k,:)/tmp(0:1,k-1,:))*  \
      conform_dims((/2,nx/),dim_avg_n(dzr3(k,:,:),0),1)
 end do
 ; if T on ...
 ;  (g/T) * (dT/dz + g/cp)

 var(1,0:1,0,:) = var(1,0:1,1,:)

 var(1,0:1,:,:) = var(1,0:1,:,:)*1.e4

 var(:,2,:,:) = var(:,1,:,:) - var(:,0,:,:)

 delete(tmp)

;--------------------------------------- cn level
; cnlev = (/ fspan(-21,21,15),  \
 cnlev = (/ fspan(-14,14,15),  \
            fspan(0,16,15) /)
 cnlevd = (/ fspan(-10,10,21),  \
             fspan(-4,4,21) /)
; cnlev = fspan(-40,28,18)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = (/20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 20 + ispan(0,2,1)*(wi+25)
 yf = 10 + ispan(0,5,1)*(he+5)

 nfig = nv

 xstr  = new(nfig,string)
 lbopt = new(nfig,string)

 mstr = ":F25:"+(/"(a)","(b)"/)
; xstr(:) = ""
; xstr(nfig-1) = "longitude [deg]"
 xstr(:) = "longitude [deg]"
 ystr = "height [km]"
 lbopt(:) = ""
 lbopt(:) = (/"[m s:S:-1:N:]","[10:S:-4:N: s:S:-2:N:]"/)


do ip=0, 0  ; 1
do ifig=0, 0  ; nfig-1

 contour1 = cnshadeplot(wks,x,z,var(ifig,ip,:,:),cnlev(ifig,:),cnfil,lbopt(ifig))
 if (ip .eq. 0) then
   cnaddline(contour1,wks,x,z,var(ifig,ip,:,:),cnlev(ifig,:),0,"D",0,"T",0.1,"C",(/colind,-1,colind/))
   cnaddline(contour1,wks,x,z,var(ifig,1,:,:),cnlev(ifig,:),0,"D",0,"T",0.8,"C",(/27,-1,29/))
 end if
 if (ip .eq. 1) then
   cnaddline(contour1,wks,x,z,var(ifig,2,:,:),cnlevd(ifig,:),0,"D",(/0,0,1/),"T",0.5,"C",(/1,-1,1/))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.05 ; 0.18 ;0.12
  "pmLabelBarParallelPosF"   : 0.50
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.018
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.015
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(ip),yf(ifig),wi,he)
 axis_range(contour1,0,360,3,18,False,False)
 axis_str(contour1,"",xstr(ifig),ystr)

 ps = (/500,300,200,150,100,70,50,40,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf(ip)+wi+2,yf(ifig),1,he)
 axis_range(plotx,0,360,3,18,False,False)
 axis_str(plotx,"","","p [hPa]")
 setvalues plotx
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : True
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tiYAxisSide" : "Right"
  "tiYAxisAngleF" : 270
  "tmYRLabelsOn" : True
  "tmYLMode" : "Explicit"
  "tmYLValues" : zs
  "tmYLLabels" : ps
 end setvalues

 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 90
  "tmXBMinorPerMajor" : 8
 end setvalues


 cnfinalize(contour1,wks)

 draw(contour1)
; draw(plotx)

; if (nfig .gt. 1) then
;   dr_txt2(wks,mstr(ifig),xf(ip)-wi-11,yf(0)-1,22,6,0.,"CenterLeft",1)
; end if

end do  ; ifig
end do  ; ip

 frame(wks)

end

