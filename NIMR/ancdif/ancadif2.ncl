load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.5

begin

;--------------------------------------- parameters
 expname0 = (/"l38","l60"/)
 expname = (/"l38","l38l60"/)

 wks = cr_wks("ancadif2",600,0)
 gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "o3_"+expname0+".nc"
   fnam2 = "oxi_"+expname+".nc"

   f10 = addfile(fnam1(0),"r")
   f11 = addfile(fnam1(1),"r")
   f20 = addfile(fnam2(0),"r")
   f21 = addfile(fnam2(1),"r")

   o3_38 = dim_avg_n( f10->O3(:,:,:,0), 0 )
   o3_60 = dim_avg_n( f11->O3(:,:,:,0), 0 )
   oh_38   = dim_avg_n( f20->field580(:,:,:,0), 0 )
   oh_60   = dim_avg_n( f21->field580(:,:,:,0), 0 )
   ho2_38  = dim_avg_n( f20->field600(:,:,:,0), 0 )
   ho2_60  = dim_avg_n( f21->field600(:,:,:,0), 0 )
   h2o2_38 = dim_avg_n( f20->field581(:,:,:,0), 0 )
   h2o2_60 = dim_avg_n( f21->field581(:,:,:,0), 0 )
   uv_38   = dim_avg_n( f20->uv      (:,:,:,0), 0 )
   uv_60   = dim_avg_n( f21->uv      (:,:,:,0), 0 )

   var_38 = (/o3_38,oh_38,ho2_38,h2o2_38,uv_38/)
   var_60 = (/o3_60,oh_60,ho2_60,h2o2_60,uv_60/)

   z_38 = f10->hybrid_ht/1.e3
   z_60 = f11->hybrid_ht/1.e3

   lat = f10->latitude

   ny = dimsizes(lat)
   nz_38 = dimsizes(z_38)
   nz_60 = dimsizes(z_60)

   wa = new(ny,float)
   wa(:) = cos(lat(:)/180.*3.141592)
   wa(:) = wa(:)/sum(wa)

   var_38a = new((/5,nz_38/),float)
   var_60a = new((/5,nz_60/),float)
   do iv=0, 4
     do k=0, nz_38-1
       var_38a(iv,k) = sum(var_38(iv,k,:)*wa(:))
     end do
     do k=0, nz_60-1
       var_60a(iv,k) = sum(var_60(iv,k,:)*wa(:))
     end do
   end do

;--------------------------------------- process 0
 cnlev = fspan(1,20,20)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = ispan(36,201,8)

;--------------------------------------- plot
 scale = (/1.e-6,1.e6,1.e7,1.e-10,1.e-6/)
 xmax = (/15,30,30,20,15/)
 ymax = (/40.,85.,85./)
 wi = (/50,50,30/)
 he = 28.*ymax/85.
 xf = 15 + (/0,1,2/)*60
 yf = 12 + (/0,1,2,3,4/)*42 + (/0,1,1,1,1/)*15

 txtstr = (/"O3","OH","HO2","H2O2","UV"/)
 xstr = (/"lat (deg)","lat (deg)","lat (deg)","lat (deg)","lat (deg)"/)
 ystr = (/"height (km)","",""/)
 xstr2 = (/"O:B:3:N:","OH","HO:B:2:N:","H:B:2:N:O:B:2:N:","UV"/)
 unit2 = " "+(/"(10:S:-6:N:)","(10:S:6:N: mol/cc)","(10:S:7:N: mol/cc)","(10:S:-10:N:)","(10:S:-6:N: m:S:2:N: s:S:-2:N:)"/)
 lbopt = (/"","","","",""/)

 do iv=0, 4

   plot = new(3,graphic)

   plot(0) = cnshadeplot(wks,lat,z_38,var_38(iv,:,:)/scale(iv),cnlev,cnfil,lbopt(iv))
   plot(1) = cnshadeplot(wks,lat,z_60,var_60(iv,:,:)/scale(iv),cnlev,cnfil,lbopt(iv))
   plot(2) = xylineplot(wks,var_60a(iv,:)/scale(iv),z_60,"D",0,"T",1.5,"C",170)
   xyaddline(plot(2),var_38a(iv,:)/scale(iv),z_38,"D",2,"T",1.5,"C",1)

   do i=0, 2
     set_size(plot(i),xf(i),yf(iv)+max(he)-he(i),wi(i),he(i))
     if (i .le. 1) then
       axis_range(plot(i),-90,90,0,ymax(i),False,False)
       axis_str(plot(i),"",xstr(iv),ystr(i))
       cnfinalize(plot(i),wks)
     else
       axis_range(plot(i),0,xmax(iv),0,ymax(i),False,False)
       axis_str(plot(i),"",xstr2(iv)+unit2(iv),ystr(i))
     end if
   end do

   draw(plot)

 end do

 dr_txt(wks,"Ozone"                          ,xf(1)+10,yf(0)-8,22,10,0.)
 dr_txt(wks,"Sulphur-cycle Chemical Oxidants",xf(1)+10,yf(1)-8,22,10,0.)
 do iv=1, 4
   dr_txt(wks,txtstr(iv),xf(1)+7,yf(iv)+4,22,8,0.)
 end do

 frame(wks)

end



