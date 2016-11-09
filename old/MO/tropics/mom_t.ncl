load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+144
 lat0 = 5.
 lat1 = 20.
 kz  = 47
 zb = 20.e3
 zt = 30.e3
 i1 = 1

 wks = cr_wks("fig/mom_t-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
   fnam1 = "../dat/u_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,{zb-2.e3:zt+3.e3},{-lat0:lat0},0) )
   z1 = f1->hybrid_ht({zb-2.e3:zt+3.e3})/1.e3

   fnam2 = "../dat/u_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,{zb-2.e3:zt+3.e3},{-lat0:lat0},0) )

   fnam1 = "../dat/cg_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   gc1w= dim_avg( f1->unspecified_2(4:nmon-1,kz,{-lat1:lat1},0) )*(-1.)*1.e3
   fnam1 = "../dat/ng_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   gn1w= dim_avg( f1->field422(4:nmon-1,kz,{-lat1:lat1},0) )*1.e3
   g1w= gc1w + gn1w

   fnam2 = "../dat/ng_"+expname(1)+".nc"
   f2 = addfile(fnam2,"r")
   gn2w = dim_avg( f2->field422(4:nmon-1,kz,{-lat1:lat1},0) )*1.e3
   g2w= gn2w

   fnam1 = "../dat/mfct0_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   dc1  = dim_avg_n( f1->unspecified_2(4:nmon-1,0,{-lat1:lat1},:), (/1,2/) )*86400.
   gc1w0= dim_avg_n( f1->unspecified_1(4:nmon-1,0,{-lat1:lat1},:), (/1,2/) )*(-1.)*1.e3
   fnam1 = "../dat/mfct_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   gc1w1= dim_avg_n( f1->unspecified_1(4:nmon-1,0,{-lat1:lat1},:), (/1,2/) )*(-1.)*1.e3

   nz0 = dimsizes(z1)
   nt = nmon-4

   u = new((/2,nz0,nt/),float)
   do k=0, nz0-1
     u(:,k,:) = (/u2(:,k),u1(:,k)/)
   end do

   uss = new((/2,nz0,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
   end do
   do n=0, nt-1
     u(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

   var = new((/6,3,nt/),float,1.e20)
   do n=0, nt-1
     do ie=0, 1
     do k=0, nz0-2
;       if ( ismissing(var(ie,0,n)) ) then
       if ( u(ie,k,n) .gt. 0. .and. u(ie,k+1,n) .le. 0. ) then
         tmp = ( z1(k+1)*u(ie,k,n) - z1(k)*u(ie,k+1,n) ) / (u(ie,k,n) - u(ie,k+1,n) )
;         if ( tmp .ge. zb/1.e3 .and. tmp .le. zt/1.e3 ) then
           var(ie,0,n) = tmp
;         end if
       end if
;       end if
     end do
     end do
   end do
   var(2:3,:,:) = (/ (/g2w,g2w*0.,g2w*0./), (/g1w,gc1w,gn1w/) /)
   var(4:5,:,:) = (/ (/g2w,g2w,g2w/)*0., (/gc1w1*0.,gc1w1,gc1w1*0./) /)

   do n=0, nt-2
   do ie=0, 1
     if ( .not. ismissing(var(ie,0,n)) .and. .not. ismissing(var(ie,0,n+1)) ) then
     if (  ( var(ie,0,n  ) .ge. zb/1.e3 .and. var(ie,0,n  ) .le. zt/1.e3 ) .or. \
           ( var(ie,0,n+1) .ge. zb/1.e3 .and. var(ie,0,n+1) .le. zt/1.e3 ) ) then
       var(ie,1,n) = (var(ie,0,n+1) - var(ie,0,n))/1.
     end if
     end if
   end do
   end do

   t = 1979 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0

;--------------------------------------- plot
 ymin = (/0,0,0,0/)
 ymax = (/-1.5,2.0,2.0,20/)
 if (i1(0) .eq. 0) then
   ymin(0) = zb/1.e3
   ymax(0) = zt/1.e3
 end if
 wi = 72
 he = (/20,20,20/)
 xf = 12 + (/1,0,0/)*78
 yf = 10 + (/0,1,2,3,4/)*26

 mstr = new((/3,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","CC"/)
 xstr = new(3,string)
 xstr(:) = ""
 xstr(2) = "year"
 ystr = (/"","Desc. rate (km mon:S:-1:N:)","","W-MF (mPa)","","W-MF (mPa)","","MDCH:S:2:N: (K:S:2:N: day:S:-2:N:)"/)


 do iv=0, 2
 do ie=1,1;0, 1

   if (iv .eq. 0) then
     plot = xylineplot(wks,t+0.5/12.,var(iv*2+ie,i1(iv),:),"D",0,"T",1.0,"C",1)
   else
     plot = xylineplot(wks,t,0,"D",0,"T",1.0,"C",1)
     xyaddline(plot,t,var(iv*2+ie,:,:),"D",(/0,0,0/),"T",1.0,"C",(/1,2,4/))
   end if
   set_size(plot,xf(ie),yf(iv),wi,he(iv))
   axis_range(plot,1979,1991,ymin(iv),ymax(iv),False,False)
   axis_str(plot,mstr(iv,ie),xstr(iv),ystr(iv*2+ie))
   setvalues plot
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 1
    "tmXBTickStartF"    : 0
    "tmXBMinorPerMajor" : 11
   end setvalues

   if ( iv .eq. 2 .and. ie .eq. 1 ) then
     setvalues plot
      "tmYROn" : False
     end setvalues

     plot2 = xylineplot(wks,t,dc1^2,"D",14,"T",1.0,"C",1)
     set_size(plot2,xf(ie),yf(iv),wi,he(iv))
     axis_range(plot2,0,12,ymin(3),ymax(3),False,False)
     axis_str(plot2,"","",ystr(3*2+ie))
     setvalues plot2
      "tmXBOn" : False
      "tmXTOn" : False
      "tmYLOn" : False
      "tmYRLabelsOn"  : True
      "tmYUseLeft"    : False
      "tiYAxisSide"   : "Right"
      "tiYAxisAngleF" : 270
     end setvalues
   end if

   gsres             = True
   gsres@gsFillColor = colind

   dummy = new(nt/2,graphic)
   tag = 0
   cnt = 0
   do n=0, nt-1
     if ( .not. ismissing(var(ie,0,n)) ) then
       if ( var(ie,0,n) .ge. zb/1.e3 .and. var(ie,0,n) .lt. zt/1.e3 ) then
         if (tag .eq. 0) then
           tag = 1
           xp0 = t(n)
         end if
         xp1 = t(n)
       else
         if (tag .eq. 1) then
           tag = 2
         end if
       end if
     else
       if (tag .eq. 1) then
         tag = 2
       end if
     end if
     if ( tag .eq. 2 .or. ( n .eq. nt-1 .and. tag .eq. 1 ) ) then
       tag = 0
       dummy(cnt) = gsn_add_polygon(wks,plot,(/xp0,xp0,xp1,xp1/),(/ymin(iv),ymax(iv),ymax(iv),ymin(iv)/),gsres)
       cnt = cnt + 1
     end if
   end do

;   dr_txt(wks,p+" hPa",xf(0)+8,yf(iv)+3,22,8,0.)

   if ( .not. ( iv .eq. 2 .and. ie .eq. 0 ) ) then
     draw(plot)
   end if
   if ( iv .eq. 2 .and. ie .eq. 1 ) then
     draw(plot2)
   end if

  end do
  end do


 frame(wks)

end



