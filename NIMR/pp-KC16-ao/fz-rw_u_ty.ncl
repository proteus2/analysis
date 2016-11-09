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
 z0    = 17.4   ; 18.59  ; 17.41  ; 16.3
 ymax0 = 0.5    ; 0.5    ; 0.8
 month = 1
 y1 = 1952 ;2058 ;1960
 y2 = 2097 ;2097 ;1999
 y0 = 1951

 x_wmrgw = 3

 wks = cr_wks("fig/fz-rw_u_ty-m"+month+"-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/"+expdir(0)+"/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/"+expdir(0)+"/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"
   fnam3 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"
print(fnam1+"")
print(fnam2+"")

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nw = 8

   nt = (y2 - y1 + 1)*12 ;647 + 1
   n1 = (y1 - y0)*12 - 1  ; Jan. of y0 is absent.
   n2 = n1 + nt - 1

;   rw_i = new((/nt-1,nw,ny/),float)
   rw_i = new((/nt,nw,ny/),float)
   rw_i(:,:,:) = f1->f_z(n1:n2,:nw-1,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     f2->f_z(n1:n2,8:9,{z0*1.e3},{-lat0:lat0})
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(n1:n2,2,{z0*1.e3},{-lat0:lat0}) + \
                      f2->f_z(n1:n2,3,{z0*1.e3},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(n1:n2,13,{z0*1.e3},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(n1:n2,16:17,{z0*1.e3},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

;   rw0 = new((/nt,nw,ny/),float)
;   rw0(:nt-2,:,:) = rw_i(:,:,:)
;   rw0(nt-1,:,:) = rw0(nt-2,:,:)
   rw0 = rw_i

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   w = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     w(:,:,:,n) = rw(:,n,:,:)
   end do

   u = transpose( f3->u(n1:n2,{z0*1.e3},{-lat0:lat0}) )

   t = ispan(y1,y2,1)

   wss = new((/3,2,12,ny,(nt/12)/),float)
   do n=0, 11
     wss(:,:,n,:,:) = w(:,:,:,ispan(n,nt-1,12))
   end do
   wss!3 = "lat"
   wss&lat = lat1

   uss = new((/12,ny,(nt/12)/),float)
   do n=0, 11
     uss(n,:,:) = u(:,ispan(n,nt-1,12))
   end do
   uss!1 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0
  wss = runave(wss,5,0)
  uss = runave(uss,5,0)

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10,12/)*0.1
 cnfil = (/4,5,6,7,8,10,12,14,15,16,17,18,19,20,21/)

 ldsh = (/0,16,14/)

;--------------------------------------- contour
 wname = (/ (/"K","E-MRG","E-IG"/), (/"R","W-MRG","W-IG"/) /)
 wi = 30
 he = 22
 xf = 12 + (/0,1,0,1/)*36
 yf = 10 + (/0,1,2,3,4/)*28 + (/0,0,0,1,1/)*12

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = ":F10:"+(/"-F:B:z:N:","F:B:z:N:"/)
 xstr = (/"","","month"/)
 ystr = (/"lat (deg)","",":F10:F:B:z:NF22: (mPa)",""/)
 txtstr = (/ (/"Kelvin","E-MRG","E-IG"/), \
             (/"Rossby","W-MRG (:F34:4:F22: "+x_wmrgw+")","W-IG"/) /)
 lbopt = new((/3,2/),string)
 lbopt(:,:) = ""
 lbopt(2,0) = "(mPa)"
 ymax = ymax0*(/-1,1/)

do idr=0, 1

do iv=0, 2

 cnvar = wss(iv,idr,month-1,:,:)
 if ( idr .eq. 0 ) then
   cnvar = cnvar*(-1.)
 end if
 if ( wname(idr,iv) .eq. "W-MRG" ) then
   cnvar = cnvar*x_wmrgw
 end if
; cnvar = abs(wss(iv,idr,:,:))

 contour1 = cnshadeplot(wks,t,lat1,cnvar,cnlev,cnfil,lbopt(iv,idr))
 cnaddline(contour1,wks,t,lat1,uss(month-1,:,:),ispan(-30,30,3),1,"D",0,"T",0.6,"C",1)
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
 axis_range(contour1,-999,999,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,idr),xstr(iv),ystr(idr))

 cnfinalize(contour1,wks)

 draw(contour1)

 dr_txt2(wks,txtstr(idr,iv),xf(idr)+0.6,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; iv

end do  ; idr

 frame(wks)

end



