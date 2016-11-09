load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"tc"/)
 expdir  = (/"AOL60CGW"/)

 lat0 = 25.
 lat0avg = 5
 z0    = 17.4  ;18.59  ; 18.59  ; 17.41
 ymax0 = 1.0    ; 1.0    ; 1.1
 y1 = 2058 ;2058 ;1960
 y2 = 2097 ;2097 ;1999
 y0 = 1951

 x_wmrgw = 3

 wks = cr_wks("fig/fz-pw_ud_my-z"+numeric2int(z0,1)+"-"+y1+"-"+y2+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.2,0.2,0.2)

;--------------------------------------- read data
   fnam1 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pd.nc"
   fnam2 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pg.nc"
   fnam3 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"
print(fnam1+"")
print(fnam2+"")

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")
   lat1 = f1->latitude({-lat0:lat0})

   nt = (y2 - y1 + 1)*12
   n1 = (y1 - y0)*12
   n2 = n1 + nt - 1

   ny = dimsizes(lat1)
   nw = 3

   rw = new((/nw,2,nt,ny/),float)
   rw(0,0,:,:) = f1->F_e_gwconv (n1:n2,{z0*1.e3},{-lat0:lat0})
   rw(0,1,:,:) = f1->F_w_gwconv (n1:n2,{z0*1.e3},{-lat0:lat0})
   rw(1,0,:,:) = f2->Fp_e_gwussp(n1:n2,{z0*1.e3},{-lat0:lat0})
   rw(1,1,:,:) = f2->Fp_w_gwussp(n1:n2,{z0*1.e3},{-lat0:lat0})*(-1.)
   rw(2,0,:,:) = rw(0,0,:,:) + rw(1,0,:,:)
   rw(2,1,:,:) = rw(0,1,:,:) + rw(1,1,:,:)

   rw(:,:,:,:) = rw(:,:,:,:)*(-1.)*1.e3

   w = new((/nw,2,ny,nt/),float)
   do n=0, nt-1
     w(:,:,:,n) = rw(:,:,n,:)
   end do

   um = f3->u(n1:n2,{4.2e3:z0*1.e3},{-lat0:lat0})
   ud = new((/nw,2,ny,nt/),float)
   do n=0, nt-1
   do j=0, ny-1
     k1 = maxind( um(n,:,j) )
     k0 = minind( um(n,:,j) )
;     ud(0,0,j,n) = um(n,k1,j) - min( um(n,:k1,j) )
;     ud(0,1,j,n) = um(n,k0,j) - max( um(n,:k0,j) )
     ud(1,0,j,n) = max( um(n,1:,j) ) ;- um(n,0,j)
     ud(1,1,j,n) = min( um(n,1:,j) ) ;- um(n,0,j)
   end do
   end do
   ud(0,:,:,:) = ud(1,:,:,:)

   t = fspan(0,13,14)

   wss = new((/nw,2,ny,14/),float)
   do n=1, 12
     wss(:,:,:,n) = dim_avg(w(:,:,:,ispan(n-1,nt-1,12)))
   end do
   wss(:,:,:,0 ) = wss(:,:,:,12)
   wss(:,:,:,13) = wss(:,:,:,1 )
   wss!2 = "lat"
   wss&lat = lat1

   uss = new((/nw,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,n) = dim_avg(ud(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )
   uss!2 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0
 do j=0, ny-1
   wss(:,:,j,:) = wss(:,:,j,:)*cos(lat1(j)*3.141592/180.)
 end do

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10,15/)*0.1
 cnfil = (/5,6,7,8,10,12,14,15,16,17,18,19,20,21,0/)

 ldsh = (/0,16,0/)
 ltck = (/1.2,1.2,2.0/)

;--------------------------------------- contour
 wname = (/ (/"E-CGW","E-BGW","E-GW"/), (/"W-CGW","W-BGW","W-GW"/) /)
 wi = 30
 he = 22
 xf = 12 + (/0,1,0,1/)*38 + (/0,0,0,0/)*8
 yf = 10 + (/0,1,2,3/)*28 + (/0,0,1,1/)*12

 mstr = new((/2,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = ":F10:"+(/"-F:B:z:N:","F:B:z:N:"/)
 xstr = (/"","","month"/)
 ystr = (/"lat (deg)","",":F10:F:B:z:NF22: (mPa)",""/)
 txtstr = wname
 lbopt = new((/2,2/),string)
 lbopt(:,:) = ""
 lbopt(1,:) = (/"(mPa)",""/)
 ymax = ymax0*(/-1,1/)

do idr=0, 1

do iv=0, 1

 contour1 = cnshadeplot(wks,t,lat1,abs(wss(iv,idr,:,:)),cnlev,cnfil,lbopt(iv,idr))
 cnaddline(contour1,wks,t,lat1,uss(iv,idr,:,:),ispan(-60,60,6)/tofloat(idr+1),1,"D",2,"T",1.2,"C",colind)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.01
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

 dr_txt2(wks,txtstr(idr,iv),xf(idr)+0.6,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; iv

 plot = xylineplot(wks,t,dim_avg_n(wss(:,idr,{-lat0avg:lat0avg},:),1),"D",ldsh,"T",ltck,"C",1)
 set_size(plot,xf(2+idr),yf(2),wi,he)
 axis_range(plot,0.5,12.5,0,ymax(idr),False,False)
 axis_str(plot,"",xstr(2),ystr(2+idr))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
; draw(plot)

 lg = lgplot(wks,wname(idr,:),"D",ldsh,"T",ltck,"C",1,2)
 set_size_lg(lg,xf(2+idr),yf(2)+he+8,wi*0.6,12)
; draw(lg)

end do  ; idr

;  txtstrb = ":F22:"+lat0avg+":S:o:N:N:NF29:-:F22:"+lat0avg+":S:o:N:S"
;  dr_txt(wks,txtstrb,xf(3)+wi-8,yf(2)+2,22,6,0.)

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,6,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(2)-2,25,6,0.)


 frame(wks)

end



