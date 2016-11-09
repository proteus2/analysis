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
 varname = "unspecified"

 lat0 = 20.
 z0   = 17.41  ; 18.59  ; 17.41
 y1 = 1955 ;2058 ;1960
 y2 = 1994 ;2097 ;1999
 y0 = 1951
 ydif = 100
 nc = 61
 dc = 2.
 ang = 45.

 wks = cr_wks("fig/fz-pw_cy-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
   fnam1 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.ph.nc"
;   fnam2 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pg.nc"
   fnam3 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"
print(fnam1+"")
;print(fnam2+"")

   f1 = addfile(fnam1,"r")
;   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")
   lat1 = f1->latitude({-lat0:lat0})
print(f1->$f1->unspecified!1$({z0*1.e3-200})+"  "+(z0*1.e3)+"  "+  \
      f1->$f1->unspecified!1$({z0*1.e3+200}))

   nt = (y2 - y1 + 1)*12
   n1 = (y1 - y0)*12
   n2 = n1 + nt - 1
   n1f = (y1+ydif - y0)*12
   n2f = n1f + nt - 1

   ny = dimsizes(lat1)
   nw = 3

   rw = new((/nw,2,ny,nc/),float)
do ic=0, nc-1
   if (ic .eq. 0) then
     varnamei0 = varname
   else
     varnamei0 = varname+"_"+ic
   end if
   varnamei1 = varname+"_"+(ic+nc)
   rw(0,0,:,ic) = dim_avg_n( 0.5*  \
                 ( f1->$varnamei0$(n1:n2,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei0$(n1:n2,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(0,1,:,ic) = dim_avg_n( (-0.5)*  \
                 ( f1->$varnamei1$(n1:n2,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei1$(n1:n2,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(1,0,:,ic) = dim_avg_n( 0.5*  \
                 ( f1->$varnamei0$(n1f:n2f,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei0$(n1f:n2f,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(1,1,:,ic) = dim_avg_n( (-0.5)*  \
                 ( f1->$varnamei1$(n1f:n2f,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei1$(n1f:n2f,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
end do
   rw = abs(rw)

;   rw(2,:,:,:) = 0.5*(rw(0,:,:,:)+rw(1,:,:,:))
   rw(2,:,:,:) = rw(0,:,:,:)

   rw@_FillValue = -999.
   rw(2,:,:,:) = where(rw(2,:,:,:) .eq. 0., rw@_FillValue, rw(2,:,:,:))

;   rw(2,0,:,:) = rw(0,0,:,:) + rw(1,0,:,:)
;   rw(2,1,:,:) = rw(0,1,:,:) + rw(1,1,:,:)

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

   t = ispan(y1,y2,1)

   wss = new((/nw,2,ny,nc/),float)
   wss = rw*1.e3
   wss!2 = "lat"
   wss&lat = lat1

   uss = dim_avg( ud )
   uss!2 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0
 do j=0, ny-1
   wss(:,:,j,:) = wss(:,:,j,:)*cos(lat1(j)*3.141592/180.)
 end do
 rwrdif = (wss(1,:,:,:) - wss(0,:,:,:))/wss(2,:,:,:)
 rwrdif = where(wss(0,:,:,:) .lt. 0.5e-4, rwrdif@_FillValue, rwrdif)

; wss = runave(wss,5,0)
; uss = runave(uss,5,0)

 c = ispan(-nc/2,nc/2,1)*dc*cos(ang*3.141592/180.)

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
; cnlev = (/0.05,0.1,0.25,0.5,1,2.5,3.5,5,6,7,8,9,10,15,20/)*0.2
 cnlev = (/0.1,0.2,0.4,0.7,1,2,4,7,10,20,30,40/)*0.1
 cnfil = (/0,5,6,7,8,12,14,15,16,17,19,20,21/)

 ldsh = (/0,16,0/)
 ltck = (/1.2,1.2,2.0/)

;--------------------------------------- contour
 wname = (/ (/"E-CGW","E-BGW","E-GW"/), (/"W-CGW","W-BGW","W-GW"/) /)
 wi = 30
 he = 22
 xf = 12 + (/0,1,2,3/)*36
 yf = 10 + (/0,1,2,3/)*28 + (/0,0,0,0/)*12

 mstr = new((/nw,2/),string)
 mstr(:,:) = ""
; mstr(0,:) = ":F10:"+(/"-F:B:z:N:","F:B:z:N:"/)
 mstr(0,:) = (/"Eastward flux","Westward flux"/)
 xstr = new(nw,string)
 xstr = ""
; xstr(nw-1) = "phase speed [m s:S:-1:N:]"
 xstr = "phase speed [m s:S:-1:N:]"
 ystr = (/"latitude [deg]","",":F10:F:B:z:NF22: [mPa]",""/)
; txtstr = wname
 txtstr = (/"CGW","BGW"/)
; txtcol = (/ (/0,0,0/), (/0,1,1/) /)
 lbopt = new((/nw,2/),string)
 lbopt(:,:) = ""
; lbopt(nw-1,0) = "[mPa]"
 lbopt(0,0) = "[10:S:-2:N: mPa / m s:S:-1:N:]"

do idr=0, 1

do iv=0, 0 ; 1

 contour1 = cnshadeplot(wks,c,lat1,abs(wss(iv,idr,:,:))*1.e2,cnlev,cnfil,lbopt(iv,idr))
 cnaddline(contour1,wks,c,lat1,rwrdif(idr,:,:)*100.,ispan(-10,10,1)*10.,1,"D",(/0,0,14/),"T",(/0.75,1.0,0.75/),"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.10
  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.4
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.025
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.01
  "lbTitleFontHeightF" : 0.03
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 if (idr .eq. 0) then
   axis_range(contour1,-15,40,-lat0,lat0,False,False)
 else
   axis_range(contour1,-40,15,-lat0,lat0,False,False)
 end if
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))

 cnfinalize(contour1,wks)

 draw(contour1)

; dr_txt2(wks,txtstr(idr,iv),xf(idr)+0.6,yf(iv)+he-2,22,5,0.,"CenterLeft",0)
;if (idr .eq. 0) then
; dr_txt2(wks,txtstr(iv),xf(idr)+0.6,yf(iv)+he-1.5,22,5,0.,"CenterLeft",0)
;end if

end do  ; iv

end do  ; idr


 frame(wks)

end

