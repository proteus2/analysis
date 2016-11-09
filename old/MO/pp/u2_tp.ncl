load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 5.
 p0 = 50.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*12

 wks = cr_wks("fig/u2_tp-"+expname(0),600,0)
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
     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

   uss0 = new((/nz3,12/),float)
   do n=0, 11
     uss0(:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,n) = ur(:,n) - uss0(:,mod(n,12))
   end do

   tss = fspan(0,13,14)

   uss2 = new((/2,nz,14/),float)
   do n=1, 12
     uss2(:,:,n) = uss(:,:,n-1) - dim_avg( uss )
   end do
   uss2(:,:,0 ) = uss2(:,:,12)
   uss2(:,:,13) = uss2(:,:,1 )

   uss20 = new((/nz3,14/),float)
   do n=1, 12
     uss20(:,n) = uss0(:,n-1) - dim_avg( uss0 )
   end do
   uss20(:,0 ) = uss20(:,12)
   uss20(:,13) = uss20(:,1 )

;--------------------------------------- process 0
 do k=0, nz-1
   if (abs(p1(k)-p0)/p0 .le. 0.1) then
     kk1 = k
   end if
 end do
 do k=0, nz3-1
   if (abs(p(k)-p0)/p0 .le. 0.1) then
     kk = k
   end if
 end do


 ts = new((/3,12/),float,1.e20)
 zz = new((/3,12/),float,1.e20)
 do iv=0, 1
   cnt = 0
   do n=0, nt-2
     if ( udes(iv,kk1,n)*udes(iv,kk1,n+1) .lt. 0. ) then
       x = ( abs(udes(iv,kk1,n))*t(n+1) + abs(udes(iv,kk1,n+1))*t(n) ) / (abs(udes(iv,kk1,n))+abs(udes(iv,kk1,n+1)))
       ts(iv,cnt) = x
       zz(iv,cnt) = 7.0*log(1.e3/p1(kk1))
       cnt = cnt + 1
     end if
   end do
 end do
 cnt = 0
 do n=0, nt3-2
   if ( udes0(kk,n)*udes0(kk,n+1) .lt. 0. ) then
     x = ( abs(udes0(kk,n))*t3(n+1) + abs(udes0(kk,n+1))*t3(n) ) / (abs(udes0(kk,n))+abs(udes0(kk,n+1)))
     ts(2,cnt) = x
     zz(2,cnt) = 7.0*log(1.e3/p(kk))
     cnt = cnt + 1
   end if
 end do

;--------------------------------------- cn level
 cnlev   = fspan(-40,40,21)
 cnlevss = fspan(-35,35,15)
 cnfil   = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfilss = (/20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3/)

; cnlevss = cnlev

;--------------------------------------- contour
 wi = 72
 he = 20
 wi2 = 20
 he2 = 20;16
 xf = 10 + (/0,1,2/)*(wi2+6)
 yf = 50 + (/0,1,2/)*(he+6)
 yf2 = yf(0) - he2 - 18

 mstr = (/"CTL","CC","MERRA"/)
 xstr = (/"","","year"/)
 lbopt = (/"","","(m s:S:-1:N:)"/)


do iv=0, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,t,zp1,udes(iv,:,:),cnlev,cnfil,lbopt(iv))
 else
   contour1 = cnshadeplot(wks,t3,zp,udes0,cnlev,cnfil,lbopt(iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(0),yf(iv),wi,he)
 axis_range(contour1,1979,1991,15,65,False,False)
 axis_str(contour1,"",xstr(iv),"height (km)")
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 1979
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)

 plot = xylineplot(wks,ts(iv,:),zz(iv,:),"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot,16,"S",1.2,"T",1.0,"C",1)

 set_size(plot,xf(0),yf(iv),wi,he)
 axis_range(plot,1979,1991,15,65,False,False)
 axis_str(plot,"","","")
 setvalues plot
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues


 draw(contour1)
 draw(plot)

 dr_txt2(wks,mstr(iv),xf(0)+1.5,yf(iv)+3,22,7,0.,"CenterLeft",1)


 ystr = (/"height (km)","",""/)
 lbopt2 = (/"","(m s:S:-1:N:)",""/)


 if (iv .le. 1) then
   contour2 = cnshadeplot(wks,tss,zp1,uss2(iv,:,:),cnlevss,cnfilss,lbopt2(iv))
 else
   contour2 = cnshadeplot(wks,tss,zp,uss20,cnlevss,cnfilss,lbopt2(iv))
 end if
 setvalues contour2
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6*wi/wi2*dimsizes(cnlevss)/dimsizes(cnlev)
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015*wi/wi2
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018*wi/wi2
 end setvalues
; setvalues contour2
;  "tiMainOffsetXF" : 0.078
;  "pmLabelBarOrthogonalPosF" : 0.12
;  "pmLabelBarParallelPosF"   : 1.1
;  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.5
;  "pmLabelBarHeightF"        : 0.12
;;  "lbLabelAutoStride" : False
;  "lbLabelFontHeightF" : 0.036
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.01
;  "lbTitleFontHeightF" : 0.0432
; end setvalues

 set_size(contour2,xf(iv),yf2,wi2,he2)
 axis_range(contour2,0.5,12.5,15,65,False,False)
 axis_str(contour2,"","month",ystr(iv))
 setvalues contour2
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
;  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 0
;  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour2,wks)

 draw(contour2)

 dr_txt2(wks,mstr(iv),xf(iv)+1,yf2+he2-2,22,6,0.,"CenterLeft",1)

end do

 dr_txt(wks,"(a)",xf(0)-5,yf2  -3,25,6,0.)
 dr_txt(wks,"(b)",xf(0)-5,yf(0)-3,25,6,0.)

 frame(wks)

end



