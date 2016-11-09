load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 5.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*12

 wks = cr_wks("fig/u_tp_web",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,:,{-lat0:lat0},0) )
   p1 = f1->p
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )

   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   nz3 = dimsizes(p)
   nt3 = nmon0

   u = new((/2,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:)/)
   end do

   ur = new((/nz3,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t = 1979 + (fspan(1,nt,nt)-0.5)/12.
   t3 = 1979 + (fspan(1,nt3,nt3)-0.5)/12.

   uss = new((/2,nz,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) ;- uss(:,:,mod(n,12))
   end do

   uss0 = new((/nz3,12/),float)
   do n=0, 11
     uss0(:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,n) = ur(:,n) ;- uss0(:,mod(n,12))
   end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72
 he = 28
 xf = 10 + (/0,0,0/)*78
 yf = 10 + (/0,1,2/)*32

 mstr = (/"CTL","CC","MERRA"/)
 xstr = (/"","","year"/)
 ystr = (/"height (km)","height (km)","height (km)"/)
 lbopt = (/"","","(m s:S:-1:N:)"/)


do iv=2, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,t,zp1,udes(iv,:,:),cnlev,cnfil,lbopt(iv))
 else
   contour1 = cnshadeplot(wks,t3,zp,udes0,cnlev,cnfil,lbopt(iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.12;20
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,1979,1991,15,65,False,False)
 axis_str(contour1,"Tropical QBO and SAO in the zonal wind",xstr(iv),ystr(iv))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 1979
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

end do


 frame(wks)

end



