load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 149;+1
 lat0 = 5.
 kz  = (/8,9/)
 kz0 = (/5,7/)

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*12

 wks = cr_wks("fig/qbo-test2",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(0:nmon-1,kz,{-lat0:lat0},0) )
   p1 = f1->p(kz)
print(p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(0:nmon-1,kz,{-lat0:lat0},0) )

   nz = dimsizes(p1)
   nt = nmon

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,kz0,{-lat0:lat0}) )
   p  = f3->levels(kz0)
print(p)

   nt3 = nmon0

   u = new((/2,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:)/)
   end do

   ur = new((/nz,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t = fspan(-4,nt-4,nt)
   t3 = fspan(0,nt3-1,nt3)

   uss = new((/3,nz,12/),float)
   do n=0, 11
     uss(0:1,:,n) = dim_avg(u(0:1,:,ispan(n+4,nt-2,12)))
     uss(2  ,:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
   end do
   umean = dim_avg(uss)

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(0:1,:,mod(n+8,12))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,n) = ur(:,n) - uss(2,:,mod(n,12))
   end do

   do iv=0, 1
   do k=0, nz-1
     x0 = 0.
     do n=0, nt-2
       if ( udes(iv,k,n)*udes(iv,k,n+1) .lt. 0. ) then
         x = ( abs(udes(iv,k,n))*t(n+1) + abs(udes(iv,k,n+1))*t(n) ) / (abs(udes(iv,k,n))+abs(udes(iv,k,n+1)))
         print(iv+"  "+p1(k)+"  "+(x-x0))
         x0 = x
       end if
     end do
   end do
   end do

   do k=0, nz-1
     x0 = 0.
     do n=0, nt3-2
       if ( udes0(k,n)*udes0(k,n+1) .lt. 0. ) then
         x = ( abs(udes0(k,n))*t3(n+1) + abs(udes0(k,n+1))*t3(n) ) / (abs(udes0(k,n))+abs(udes0(k,n+1)))
         print(p(k)+"  "+(x-x0))
         x0 = x
       end if
     end do
   end do

;--------------------------------------- process 0

;--------------------------------------- plot
 ymin = (/-25,-50,-50,-50,-50/)
 ymax = (/25,40,40,40,40/)
 wi = 72.*(/1.,1.,1./)
 he = 20
 xf = 12 + (/0,1,2/)*78
 yf = 10 + (/0,1,2,3,4/)*31

 mstr = new((/nz,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = new(3,string)
 xstr(:) = ""
 xstr(2) = "month"
 ystr = (/"U (m s:S:-1:N:)","",""/)


 do k=0, nz-1
 do iv=0, 2

   if (iv .le. 1) then
     plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t,udes(iv,k,:),"D",0,"T",0.5,"C",1)
   else
     plot = xylineplot(wks,t3,0.,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t3,udes0(k,:),"D",0,"T",0.5,"C",1)
   end if

   set_size(plot,xf(k),yf(iv),wi(iv),he)
   axis_range(plot,-4,144,ymin(k),ymax(k),False,False)
;   axis_range(plot,82,105,ymin(k),ymax(k),False,False)
   axis_str(plot,mstr(k,iv),xstr(iv),ystr(k))
   setvalues plot
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 12
    "tmXBTickStartF"    : -12
    "tmXBMinorPerMajor" : 11
   end setvalues

   dr_txt(wks,p(k)+" hPa",xf(k)+8,yf(0)+3,22,8,0.)

   draw(plot)

  end do
  end do


 frame(wks)

end



