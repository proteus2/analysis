load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 y1 = 1953
 nmon = 54*12
 lat0 = 25.
 lat0avg = 5
 z0 = 17.5e3

 wks = cr_wks("fig/rwmf100des_t-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/epf/uanuj.epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   lat1 = f1->lat({-lat0:lat0})
   z  = f1->z / 1.e3

   ny = dimsizes(lat1)
   nz = dimsizes(z)
   nt = 648
   nw = dimsizes(f1->wg)

   rw_i = f1->f_z(:,:,{z0},{-lat0:lat0}) / 6370.e3 * 1.e3
   rw0 = new((/nt,nw,ny/),float)
   rw0(:646,:,:) = rw_i(:646,:,:)
   rw0(647,:,:) = rw0(646,:,:)

   rw0a = rw0
   rw0a(1:646,:,:) = 0.5*(rw0(0:645,:,:) + rw0(2:647,:,:))

   rw0 = where(abs(rw0) .ge. 1.e10, rw0a, rw0)

   rw = (/ rw0(:,4:5,:)+rw0(:,6:7,:), rw0(:,8:9,:), rw0(:,0:1,:) /)

   rw(:,:,0,:) = rw(:,:,0,:)*(-1.)

   u = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     u(:,:,:,n) = rw(:,n,:,:)
   end do

   uss = new((/3,2,ny,14/),float)
   do n=1, 12
     uss(:,:,:,n) = dim_avg(u(:,:,:,ispan(n-1,nt-1,12)))
   end do
   uss(:,:,:,0 ) = uss(:,:,:,12)
   uss(:,:,:,13) = uss(:,:,:,1 )
   uss!2 = "lat"
   uss&lat = lat1

   u!2 = "lat"
   u&lat = lat1
   ussa = dim_avg_n(uss(:,:,{-lat0avg:lat0avg},:),2)
   udes = dim_avg_n(u  (:,:,{-lat0avg:lat0avg},:),2)
   do n=0, nt-1
     udes(:,:,n) = udes(:,:,n) - ussa(:,:,mod(n+1,12))
   end do

   t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 udes = runave(udes,3,0)

;--------------------------------------- cn level
; cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
; cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.5/)
 cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)

;--------------------------------------- contour
 wi = 140
 he = 24
 xf = 12
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = (/"Deseasonalized F_z","",""/)
 xstr = (/"","","year"/)
 ystr = (/"F_z (hPa)"/)
 txtstr = (/"Large-scale waves","Inertia-GWs","Total resolved waves"/)


do iv=0, 2


 plot = xylineplot(wks,t,udes(iv,:,:),"D",0,"T",1.2,"C",(/3,23/))
 xyaddline(plot,t,0,"D",0,"T",0.8,"C",1)
 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,1953,2007,-0.3,0.4,False,False)
 axis_str(plot,mstr(iv),xstr(iv),ystr)
; setvalues plot
;  "tmXBMode" : "Explicit" ;"Manual"
;  "tmXBValues" : ispan(1,12,1)
;  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
; end setvalues
 draw(plot)

  dr_txt2(wks,txtstr(iv),xf+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

end do  ; iv

;  dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
;  dr_txt(wks,"(b)",xf(0)-8,yf(1)-2,25,8,0.)
;  dr_txt(wks,"(c)",xf(0)-8,yf(2)-2,25,8,0.)


 frame(wks)

end



