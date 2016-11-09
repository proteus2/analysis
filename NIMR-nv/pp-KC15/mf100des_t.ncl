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
 lat0 = 30.
 lat0avg = 5
 z = 15.e3

 wks = cr_wks("fig/mf100des_t-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/mmxm.uanuja.pd.nc"

   f1 = addfile(fnam1,"r")
   mfe1 = f1->F_e_gwconv(:,{z},{-lat0:lat0})
   mfw1 = f1->F_w_gwconv(:,{z},{-lat0:lat0})*(-1.)
   lat1 = f1->latitude({-lat0:lat0})

   fnam2 = "../dat/mmgwusspz_xm.uanuja.pg.nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->Fp_e_gwussp(:,{z},{-lat0:lat0})
   mfw2 = f2->Fp_w_gwussp(:,{z},{-lat0:lat0})
;   mfe2 = mfe2 + mfe1
;   mfw2 = mfw2 + mfw1

   ny = dimsizes(lat1)
   nt = nmon

   u = new((/3,2,ny,nt/),float)
   do n=0, nt-1
     u(0:1,:,:,n) = (/ (/mfe1(n,:),mfw1(n,:)/), \
                       (/mfe2(n,:),mfw2(n,:)/) /)*1.e3
   end do
   u(2,:,:,:) = u(0,:,:,:) + u(1,:,:,:)

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

 mstr = (/"Deseasonalized MF","",""/)
 xstr = (/"","","year"/)
 ystr = (/"MF (hPa)"/)
 txtstr = (/"CGW", "BGW", "Total GW"/)


do iv=0, 2


 plot = xylineplot(wks,t,udes(iv,:,:),"D",0,"T",1.2,"C",(/3,23/))
 xyaddline(plot,t,0,"D",0,"T",0.8,"C",1)
 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,1953,2007,-0.2,0.25,False,False)
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



