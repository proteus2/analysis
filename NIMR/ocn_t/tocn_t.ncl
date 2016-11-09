load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.5

begin

;--------------------------------------- parameters
 expname = (/"b","d","e"/)
 stime = "i512-"+(/"o012","l712","l712"/)

 ismon = 1  ; 1 for Dec.

 kz  = (/0,12,25,38/)

 wks = cr_wks("tocn_t",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/pt_"+expname+"_"+stime+".nc"

   f1 = addfile(fnam1(0),"r")
   f2 = addfile(fnam1(1),"r")
   f3 = addfile(fnam1(2),"r")

   z = f1->depth(kz)
   print(z)

   ny = dimsizes(f1->latitude)
   nz = dimsizes(kz)
   ne = dimsizes(expname)

   nmon_all = new(ne,integer)
   nmon_all = (/dimsizes(f1->t),dimsizes(f2->t),dimsizes(f3->t)/) - ismon
   nmon = max(nmon_all)
   nt = nmon/12

   toa = new((/ne,nmon,nz,ny/),float,1.e20)
   toa(0,:nmon_all(0)-1,:,:) = f1->temp(ismon:,kz,:,0)
   toa(1,:nmon_all(1)-1,:,:) = f2->temp(ismon:,kz,:,0)
   toa(2,:nmon_all(2)-1,:,:) = f3->temp(ismon:,kz,:,0)

   f0 = addfile("../dat/area_ocn.nc","r")
   wa = f0->wa4xavg({z},:)

   toa_a = new((/nz,ne,nmon/),float,1.e20)
   do i=0, ne-1
   do n=0, nmon-1
     toa_a(:,i,n) = dim_sum( toa(i,n,:,:)*wa(:,:) )
   end do
   end do

   to = new((/nz,ne*3,nt/),float,1.e20)
   do ie=0, 2
   do n=0, nt-1
     n1 = n*12
     n2 = (n+1)*12-1
     to(:,ie  ,n) = dim_avg( toa_a(:,ie,n1:n2) )
     to(:,ie+3,n) = dim_max( toa_a(:,ie,n1:n2) )
     to(:,ie+6,n) = dim_min( toa_a(:,ie,n1:n2) )
   end do
   end do

   t = fspan(1,nt,nt) - 0.5

;--------------------------------------- process 0
 ldsh = (/0,0,0, 2,2,2, 2,2,2/)
 ltck = (/1.8,1,1, 1.8,1,1, 1.8,1,1/)
 lcol = (/1,2,3, 1,2,3, 1,2,3/)

;--------------------------------------- plot
 ymin = (/17.2,13.4,4.6,1.28/)
 ymax = (/18.8,13.9,4.8,1.38/)
 wi = 70
 he = 42
 xf = 15 + (/0,1,0,1/)*80
 yf = 10 + (/0,0,1,1/)*50

 xstr = (/"","","time","time"/)
 ystr = (/"T (:S:o:N:C)","","T (:S:o:N:C)",""/)

 do k=0, nz-1

   if (k .lt. 2) then
     plot = xylineplot(wks,t,to(k,:,:),"D",ldsh,"T",ltck,"C",lcol)
   else
     plot = xylineplot(wks,t,to(k,0:2,:),"D",ldsh,"T",ltck,"C",lcol)
;     xyaddline(plot,t,to(k,:,:),"D",0,"T",0.5,"C",1)
   end if

   set_size(plot,xf(k),yf(k),wi,he)
   axis_range(plot,-999,999,ymin(k),ymax(k),False,False)
   axis_str(plot,"",xstr(k),ystr(k))

   dr_txt(wks,z(k)+" m",xf(k)+wi-10,yf(k)+3,22,8,0.)

   draw(plot)

  end do


 frame(wks)

end



