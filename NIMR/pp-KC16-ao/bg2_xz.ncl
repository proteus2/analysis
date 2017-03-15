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
 varname = (/"u","theta","p"/)

 y1 = 1955
 y2 = y1 - 1 + 40
 ydif = 100
 mons = 13   ; 1--12 : each month  13 : annual mean  ; 20--23 : djf--son
 lat0 = 15.

 wks = cr_wks("fig/bg2_xz-"+expname(0),600,0)
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
 tmp = f0->$varname(1)$(0,:,0,0)
 zt = f0->$tmp!0$(:nz-1) / 1.e3
 delete(tmp)

 ft0 = addfile("/hippo0/HG2CMIP/L60CGW/invariant/L60_z.nc","r")
 ; z at theta lev
 dzr3 = 1./(ft0->ht(0,1:nz,{-lat0:lat0},:) - ft0->ht(0,0:nz-1,{-lat0:lat0},:))

 var = new((/nv,3,nz,nx/),float)
 var = 0.
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
     do iv=0, nv-2
       var(iv,ip,:,:) = var(iv,ip,:,:) + dim_avg_n(  \
          f0->$varname(iv)$(0,:nz-1,{-lat0:lat0},:), 1 )/tofloat(nt)
     end do
     iv = nv-1
     ptmp = f0->$varname(iv)$(0,:nz-1,{-lat0:lat0},:)
     ptmp = (ptmp/1.e5)^0.286
     ttmp = f0->$varname(1)$(0,:nz-2,{-lat0:lat0},:)*  \
            (0.5*( ptmp(:nz-2,:,:) + ptmp(1:nz-1,:,:) ))
     var(iv,ip,:nz-2,:) = var(iv,ip,:nz-2,:) + dim_avg_n(  \
          ttmp, 1 )/tofloat(nt)
     var(iv,ip,nz-1,:) = var(iv,ip,nz-2,:)
   end do 
 end do 
 end do 
 delete(ptmp)
 delete(ttmp)

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
            fspan(0,16,15),  \
            fspan(190,350,15) /)
 cnlevd = (/ fspan(-10,10,21),  \
             fspan(-4,4,21),  \
             fspan(-4,4,21) /)
; cnlev = fspan(-40,28,18)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = (/20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 20 + ispan(0,2,1)*(wi+8)
 yf = 10 + ispan(0,5,1)*(he+5)*0
 wi2 = wi*1.2

 nfig = nv

 xstr  = new(nfig,string)
 ystr  = new(nfig,string)
 lbopt = new(nfig,string)

 mstr = ":F25:"+(/"(a)","(b)","(c)"/)
; xstr(:) = ""
; xstr(nfig-1) = "longitude [deg]"
 xstr(:) = "longitude [deg]"
 ystr(:) = ""
 ystr(0) = "height [km]"
 lbopt(:) = ""
 lbopt(:) = (/"[m s:S:-1:N:]","[10:S:-4:N: s:S:-2:N:]","[K]"/)

; txtstr = (/y1+"-"+y2,(y1+ydif/2)+"-"+(y2+ydif/2),(y1+ydif)+"-"+(y2+ydif)/)
 txtstr = (/y1+"-"+y2,(y1+ydif)+"-"+(y2+ydif)/)

do ip=0, 0  ; 1
do ifig=0, 0  ; nfig-1

 if (ifig .ne. 2) then
   zd = z
 else
   zd = zt
 end if

 contour1 = cnshadeplot(wks,x,zd,var(ifig,ip,:,:),cnlev(ifig,:),cnfil,lbopt(ifig))
 if (ip .eq. 0) then
   cnaddline(contour1,wks,x,zd,var(ifig,ip,:,:),cnlev(ifig,:),0,"D",0,"T",0.1,"C",(/colind,-1,colind/))
   cnaddline(contour1,wks,x,zd,var(ifig,1,:,:),cnlev(ifig,:),0,"D",0,"T",0.8,"C",(/27,-1,29/))
 end if
 if (ip .eq. 1) then
   cnaddline(contour1,wks,x,zd,var(ifig,2,:,:),cnlevd(ifig,:),0,"D",(/0,0,1/),"T",0.5,"C",(/1,-1,1/))
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

 set_size(contour1,xf(ifig),yf(ifig),wi,he)
 axis_range(contour1,0,360,0,20,False,False)
 axis_str(contour1,"",xstr(ifig),ystr)

 ps = (/500,300,200,150,100,70,50,40,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf(ifig)+wi+2,yf(ifig),1,he)
 axis_range(plotx,0,360,0,20,False,False)
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

 plot = xylineplot(wks,dim_avg(var(2,ip,:,:)),zt,"D",0,"T",1.2,"C",1)
; xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(il))
 xyaddline(plot,dim_avg(var(2,1,:,:)),zt,"D",16,"T",1.2,"C",1)

 set_size(plot,xf(1),yf(0),wi2,he)
 axis_range(plot,190,305,0,20,False,False)
 axis_str(plot,"","T [K]","")

 draw(plot)

 dr_txt2(wks,mstr(0),xf(0)-5,yf(0)-1,22,6,0.,"CenterLeft",1)
 dr_txt2(wks,mstr(1),xf(1)-5,yf(0)-1,22,6,0.,"CenterLeft",1)

 lg = lgplot(wks,txtstr,"D",(/0,16/),"T",1.2,"C",1,2)
 set_size_lg(lg,xf(1)+wi2*0.5,yf(0),wi2*0.3,9)
 draw(lg)

 delete(zd)

end do  ; ifig
end do  ; ip

 frame(wks)

end

