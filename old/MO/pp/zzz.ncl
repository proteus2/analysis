load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

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

 wks = cr_wks("fig/zzz-"+expname(0),600,0)
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
win=win*0.+1.
winr=winr*0.+1.

 mswin = avg((win)^2)
 mswinr = avg((winr)^2)

 do n=0, nt-1
   u(:,:,n) = u(:,:,n)*win(n)
   ur(:,:,n) = ur(:,:,n)*winr(n)
 end do

;   nf  = nt/2
;   nf3 = nt3/2
;   fr  = fspan(1,nf ,nf )/(int2flt(nt )/12.)
;   fr3 = fspan(1,nf3,nf3)/(int2flt(nt3)/12.)

;   fftu  = ezfftf(u )
;   fftur = ezfftf(ur)

   k  = 9 ;12
   k0 = 7 ;9
print(zp1(k)+"  "+zp(k0))
   jtot = floattoint(1 + ((log(nt*0.5))/0.25)/log(2.))
   w1 = wavelet(u (0,k ,:),0,1,6.0,2.,0.25,jtot,nt+36, \
                     1,0,0.05,0)
   w2 = wavelet(u (1,k ,:),0,1,6.0,2.,0.25,jtot,nt+36, \
                     1,0,0.05,0)
   wr = wavelet(ur(0,k0,:),0,1,6.0,2.,0.25,jtot,nt+36, \
                     1,0,0.05,0)

   pws1 = ( w1(0,:,:)^2 + w1(1,:,:)^2 )*0.5*nt /12. / mswin
   pws2 = ( w2(0,:,:)^2 + w2(1,:,:)^2 )*0.5*nt /12. / mswin
   pws = (/pws1,pws2/)
   pwsr = ( wr(0,:,:)^2 + wr(1,:,:)^2 )*0.5*nt3/12. / mswinr

pr = wr@period;/12.
fr = 1./pr

; fftu (:,:,:,7:) = 0.
; fftur(:,:,:,7:) = 0.

; urec  = ezfftb(fftu ,fftu@xbar)
; urrec = ezfftb(fftur,fftur@xbar)

; do n=1, nt-1
;   urec(:,:,n) = urec(:,:,n)/win(n)
;   urrec(:,:,n) = urrec(:,:,n)/winr(n)
; end do

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(urec(:,:,ispan(n,nt-1,12)))
;   end do

;   udes = new(dimsizes(urec),float)
;   do n=0, nt-1
;     udes(:,:,n) = urec(:,:,n) - uss(:,:,mod(n,12))
;   end do

;   uss0 = new((/nz3,12/),float)
;   do n=0, 11
;     uss0(:,n) = dim_avg(urrec(0,:,ispan(n,nt3-1,12)))
;   end do

;   udes0 = new(dimsizes(urrec(0,:,:)),float)
;   do n=0, nt3-1
;     udes0(:,n) = urrec(0,:,n) - uss0(:,mod(n,12))
;   end do

   t = fspan(0,nt/12-1,nt)+1979
   t3 = fspan(0,nt3/12-1,nt3)+1979

;--------------------------------------- cn level
 cnlev = (/2,5,10,20,50,100,200/)*0.1
 cnfil = (/0,10,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 42
 he = 24
 xf = 12 + (/0,1,2/)*50 + (/0,0,1/)*8
 yf = 10 + (/0,1,2/)*30

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = (/"","","frequency (cycle yr:S:-1:N:)"/)
 ystr = (/"log-p height (km)","","log-p height (km)"/)
 lbopt = new((/3,3/),string)
 lbopt(:,:) = ""
 lbopt(2,:) = (/"(10:S:2:N: m:S:2:N: s:S:-2:N: / cycle yr:S:-1:N:)","",""/) ;(/"(m s:S:-1:N:)","",""/)


;--------------------------------------- cn level
 cnlev2 = (/2,5,10,20,50,100,200,500/)*0.1*10.
 cnfil2 = (/0,10,7,6,4,2,26,24,22,21,19,17/)

;--------------------------------------- contour
 wi = 72
 he = 20
 xf = 10 + (/0,0,0/)*78
 yf = 50 + (/0,1,2/)*32

 xstr = (/"","","time"/)
 ystr = (/"frequency (cycle yr:S:-1:N:)","frequency (cycle yr:S:-1:N:)","frequency (cycle yr:S:-1:N:)"/)
 lbopt(0,:) = (/"","","(m s:S:-1:N:)"/)


do iv=0, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,t,pr,pws(iv,:,:)/1.e2,cnlev2,cnfil2,lbopt(0,iv))
 else
   contour1 = cnshadeplot(wks,t3,pr,pwsr/1.e2,cnlev2,cnfil2,lbopt(0,iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.5
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.5
  "pmLabelBarHeightF"        : 0.08
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.008
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.01
 end setvalues

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,-999,999,999,4,False,True)
 axis_str(contour1,mstr(0,iv),xstr(iv),ystr(iv))
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



