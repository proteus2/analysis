load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 expname = (/"tk"      ,"tc"    /)
 nick    = (/"AOL60CTL","AOL60CGW"/)

 lat0 = 5.

 y1 = 1960 ;2058 ;1960
 y2 = 1999 ;2097 ;1999
 y0 = 1951

 wks = cr_wks("fig/uvar_fz-"+y1+"-"+y2,600,0)
; wks = cr_wks("z",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/"+nick+"/xm.uan"+expname+"a.pj.nc"

   nt = (y2 - y1 + 1)*12
   it1 = (y1 - y0)*12
   it2 = it1 + nt - 1

   f1 = addfile(fnam1(0),"r")
   u1 = dim_avg( f1->u(it1:it2,:,{-lat0:lat0}) )
   lat1 = f1->latitude({-lat0:lat0})
   z1 = f1->$f1->u!1$ / 1.e3

   ny = dimsizes(lat1)
   nz = dimsizes(z1)

   f2 = addfile(fnam1(1),"r")
   u2 = dim_avg( f2->u(it1:it2,:,{-lat0:lat0}) )

   u = new((/dimsizes(expname),nz,nt/),float)
   do k=0, nz-1
     u(0,k,:) = u1(:,k)
     u(1,k,:) = u2(:,k)
   end do

;--------------------------------------- process 0
 u = dtrend(u ,False)

 win = sin(3.141592*fspan(0,1,nt))        ; cosine window

 mswin = avg((win)^2)

 do n=0, nt-1
   u(:,:,n) = u(:,:,n)*win(n)
 end do

   nf = nt/2
   fr = fspan(1,nf ,nf )/(int2flt(nt)/12.)

   fftu  = ezfftf(u )

   pws  = ( fftu (0,:,:,:)^2 + fftu (1,:,:,:)^2 )*0.5*nt/12. / mswin

;--------------------------------------- cn level
 cnlev = (/2,5,10,20,50,100,200,500/)*0.1
 cnfil = (/0,10,7,6,4,2,26,24,22,21,19/)

;--------------------------------------- contour
 wi = 30
 he = 20
 xf = 12 + (/0,1,2/)*36
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new(2,string)
 mstr = ""
 mstr = (/"CTL","CC"/)
 xstr = "frequency (cycle yr:S:-1:N:)"
 ystr = (/"",""/)
 ystr(0) = (/"height (km)"/) ;"log-p height (km)"/)
 lbopt = (/"(10:S:2:N: m:S:2:N: s:S:-2:N: / cycle yr:S:-1:N:)",""/) ;(/"(m s:S:-1:N:)","",""/)


do iv=0, dimsizes(expname)-1

 contour1 = cnshadeplot(wks,fr,z1,pws(iv,:,:)/1.e2,cnlev,cnfil,lbopt(iv))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.11
  "pmLabelBarParallelPosF"   : 1.0
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.9
  "pmLabelBarHeightF"        : 0.2
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(iv),yf,wi,he)
 axis_range(contour1,0,2.1,15,50,False,False)
 axis_str(contour1,mstr(iv),xstr,ystr(iv))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF" : 0.5
  "tmXBMinorPerMajor" : 4
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

; dr_txt2(wks,mstr(iv),xf(iv)+wi-1,yf+he-2,22,6,0.,"CenterRight",1)

end do


 frame(wks)

end



