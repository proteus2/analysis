load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+12*12
 lat0 = 5.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 m2 = 0*12
 m3 = 0*12
 nmon0 = 12*12

 wks = cr_wks("fig/uvar_fp-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,:,{-lat0:lat0},0) )
   lat1 = f1->latitude({-lat0:lat0})
   p1 = f1->p
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )

   ny = dimsizes(lat1)
   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
   u4 = dim_avg( f3->u(m2:m2+nmon0-1,:,{-lat0:lat0}) )
   u5 = dim_avg( f3->u(m3:m3+nmon0-1,:,{-lat0:lat0}) )
   lat = f3->latitude({-lat0:lat0})
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   ny3 = dimsizes(lat)
   nz3 = dimsizes(p)
   nt3 = nmon0

   u = new((/2,nz,nt/),float)
   do k=0, nz-1
     u(:,k,:) = (/u2(:,k),u1(:,k)/)
   end do

   ur = new((/3,nz3,nt3/),float)
   do k=0, nz3-1
     ur(:,k,:) = (/u3(:,k),u4(:,k),u5(:,k)/)
   end do

;--------------------------------------- process 0
 u = dtrend(u ,False)
 ur = dtrend(ur,False)

 win = sin(3.141592*fspan(0,1,nt))        ; cosine window
 winr = sin(3.141592*fspan(0,1,nt3)) 

 mswin = avg((win)^2)
 mswinr = avg((winr)^2)

 do n=0, nt-1
   u(:,:,n) = u(:,:,n)*win(n)
 end do
 do n=0, nt3-1
   ur(:,:,n) = ur(:,:,n)*winr(n)
 end do

   nf  = nt/2
   nf3 = nt3/2
   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
   fr3 = fspan(1,nf3,nf3)/(int2flt(nt3)/12.)
print(fr(0)+"  "+fr3(0))

   fftu  = ezfftf(u )
   fftur = ezfftf(ur)

   pws  = ( fftu (0,:,:,:)^2 + fftu (1,:,:,:)^2 )*0.5*nt /12. / mswin
   pwsr0= ( fftur(0,:,:,:)^2 + fftur(1,:,:,:)^2 )*0.5*nt3/12. / mswinr
   pwsr = dim_avg_n(pwsr0,0)

;--------------------------------------- cn level
 cnlev = (/2,5,10,20,50,100,200/)*0.1
 cnfil = (/0,10,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 30
 he = 16
 xf = 12 + (/0,1,2/)*38 + (/0,0,0/)*6.6
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","CC","MERRA"/)
 xstr = (/"","","frequency (cycle yr:S:-1:N:)"/)
 ystr = (/"height (km)","",""/) ;"log-p height (km)"/)
 lbopt = new((/3,3/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(10:S:2:N: m:S:2:N: s:S:-2:N: / cycle yr:S:-1:N:)","",""/) ;(/"(m s:S:-1:N:)","",""/)


do iv=0, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,fr,zp1,pws(iv,:,:)/1.e2,cnlev,cnfil,lbopt(2,iv))
 end if
 if (iv .eq. 2) then
   contour1 = cnshadeplot(wks,fr3,zp,pwsr/1.e2,cnlev,cnfil,lbopt(2,iv))
 end if
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.2
  "pmLabelBarHeightF"        : 0.2
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.25
  "lbTitleFontHeightF" : 0.038
 end setvalues

 set_size(contour1,xf(iv),yf,wi,he)
 axis_range(contour1,-999,2.5,15,65,False,False)
 axis_str(contour1,"",xstr(2),ystr(iv))

 cnfinalize(contour1,wks)


 draw(contour1)

 dr_txt2(wks,mstr(0,iv),xf(iv)+wi-1,yf+he-2,22,6,0.,"CenterRight",1)

end do


 frame(wks)

end



