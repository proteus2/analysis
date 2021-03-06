load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.1

begin

;--------------------------------------- parameters
 expname = (/"tc"/)

 imon0 = 2*12
 nmon = 18*12
 lat0 = 30.
 lat0avg = 5
 z = 17.41e3

 wks = cr_wks("fig/gwmfss_ty-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/AOL60CGW/mmxm.uan"+expname+"a.pd.nc"

   f1 = addfile(fnam1,"r")
   mfe1 = f1->F_e_gwconv(imon0:,{z},{-lat0:lat0})
   mfw1 = f1->F_w_gwconv(imon0:,{z},{-lat0:lat0})*(-1.)
   lat1 = f1->latitude({-lat0:lat0})

   fnam2 = "../dat/AOL60CGW/mmxm.uan"+expname+"a.pg.nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->Fp_e_gwussp(imon0:,{z},{-lat0:lat0})
   mfw2 = f2->Fp_w_gwussp(imon0:,{z},{-lat0:lat0})
;   mfe2 = mfe2 + mfe1
;   mfw2 = mfw2 + mfw1

   ny = dimsizes(lat1)
   nt = nmon

   u = new((/2,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = (/ (/mfe1(n,:),mfw1(n,:)/), \
                     (/mfe2(n,:),mfw2(n,:)/) /)*1.e3
   end do

   t = fspan(0,13,14)

   uss = new((/3,2,ny,14/),float)
   do n=1, 12
     uss(0:1,:,:,n) = dim_avg(u(0:1,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(2,:,:,:) = uss(0,:,:,:) + uss(1,:,:,:)
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )
   uss!2 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0

;--------------------------------------- cn level
; cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
; cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7/)
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18, 17, 16, 15, 14/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38 + (/0,0,1/)*8
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Eastward-MF","Westward-MF",lat0avg+"N:S:o:NF29:-:F22:"+lat0avg+":S:o:N:S"/)
 xstr = (/"","month",""/)
 ystr = (/"lat (deg)","","MF (mPa)"/)
 txtstr = (/"CGW", "BGW", "BGW (CTL)"/)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(1,:) = (/"(mPa)",""/)


do iv=0, 1

do idr=0, 1

 contour1 = cnshadeplot(wks,t,lat1,uss(iv,idr,:,:),cnlev,cnfil,lbopt(iv,idr))
; setvalues contour1
;  "cnFillMode" : "RasterFill"
; end setvalues
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.18 + 0.2
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6 + 0.4
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02 - 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

  dr_txt2(wks,txtstr(iv),xf(0)+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; idr

 plot = xylineplot(wks,t,dim_avg_n(uss(iv*2,:,{-lat0avg:lat0avg},:),1),"D",(/0,14/),"T",1.2,"C",(1-iv)*2+1)
 set_size(plot,xf(2),yf(0),wi,he)
 axis_range(plot,0.5,12.5,0,1.20,False,False)
 axis_str(plot,mstr(iv,2),xstr(iv),ystr(2))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

end do  ; iv

 lg = lgplot(wks,(/"All E-GW","All W-GW","E-CGW","W-CGW"/),"D",(/0,14,0,14/),"T",1.2,"C",(/1,1,3,3/),2)
 set_size_lg(lg,xf(2),yf(1)+3,wi-6,20)
; setvalues lg
;  "lgLabelOffsetF" : 0.2
;  "lgLineDashSegLenF"  : 0.030
;  "lgLabelFontHeightF" : 0.010
; end setvalues
 draw(lg)

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,8,0.)
;  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end

