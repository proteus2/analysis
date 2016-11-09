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
 zb = 20.
 zt = 30.
 iz = 1

 wks = cr_wks("fig/desc",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
print(zb+"  "+zt)
   pb_2 = 1.e3*exp(-(zb-2.0)/7.0)
   pt_2 = 1.e3*exp(-(zt+3.0)/7.0)

   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,{pb_2:pt_2},{-lat0:lat0},0) )
   p1 = f1->p({pb_2:pt_2})
   zp1 = 7.0*log(1.e3/p1)
print(p1+"  "+zp1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,{pb_2:pt_2},{-lat0:lat0},0) )

   fnam0 = "../dat/merra.u_ypt.1979-2010.nc"

   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(0:nmon-5,{pb_2:pt_2},{-lat0:lat0}) )
   p  = f3->levels({pb_2:pt_2})
   zp = 7.0*log(1.e3/p)
print(p+"  "+zp)

   nz0 = dimsizes(zp1)
   nzr = dimsizes(zp)
   nt = nmon-4

   u = new((/2,nz0,nt/),float)
   do k=0, nz0-1
     u(:,k,:) = (/u2(:,k),u1(:,k)/)
   end do
   ur = new((/nzr,nt/),float)
   do k=0, nzr-1
     ur(k,:) = u3(:,k)
   end do

   uss = new((/2,nz0,12/),float)
   ussr= new((/nzr,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
     ussr(:,n) = dim_avg(ur(:,ispan(n,nt-1,12)))
   end do
   do n=0, nt-1
     u(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
     ur(:,n) = ur(:,n) - ussr(:,mod(n,12))
   end do

   var = new((/3,2,nt/),float,1.e20)
   do n=0, nt-1
     do ie=0, 1
     do k=0, nz0-2
;       if ( ismissing(var(ie,0,n)) ) then
       if ( u(ie,k,n) .gt. 0. .and. u(ie,k+1,n) .le. 0. ) then
         tmp = ( zp1(k+1)*u(ie,k,n) - zp1(k)*u(ie,k+1,n) ) / (u(ie,k,n) - u(ie,k+1,n) )
;         if ( tmp .ge. zb .and. tmp .le. zt ) then
           var(ie,0,n) = tmp
;         end if
       end if
;       end if
     end do
     end do
     do k=0, nzr-2
;       if ( ismissing(var(2,0,n)) ) then
       if ( ur(k,n) .gt. 0. .and. ur(k+1,n) .le. 0. ) then
         tmp = ( zp(k+1)*ur(k,n) - zp(k)*ur(k+1,n) ) / (ur(k,n) - ur(k+1,n) )
;         if ( tmp .ge. zb .and. tmp .le. zt ) then
           var(2,0,n) = tmp
;         end if
       end if
;       end if
     end do
   end do

   do n=0, nt-2
   do ie=0, 2
     if ( .not. ismissing(var(ie,0,n)) .and. .not. ismissing(var(ie,0,n+1)) ) then
     if (  ( var(ie,0,n  ) .ge. zb .and. var(ie,0,n  ) .le. zt ) .or. \
           ( var(ie,0,n+1) .ge. zb .and. var(ie,0,n+1) .le. zt ) ) then
       var(ie,1,n) = (var(ie,0,n+1) - var(ie,0,n))/1.
     end if
     end if
   end do
   end do

   t = 0 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0

;--------------------------------------- plot
 ymin = (/0,0.,-0.5/)
 ymax = (/-1.5,1.8,0.5/)
 if (iz .eq. 0) then
   ymin(0) = zb
   ymax(0) = zt
 end if
 wi = 72
 he = (/20,20,20/)
 xf = 12 + (/0,0,0/)*78
 yf = 10 + (/0,1,2/)*36

 mstr = new((/3,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","CC","MERRA"/)
 xstr = new(3,string)
 xstr(:) = ""
 xstr(2) = "year"
 ystr = (/"descent rate (km mon:S:-1:N:)","","(m s:S:-1:N: mon:S:-1:N:)","","(m s:S:-1:N: mon:S:-1:N:)",""/)


; do iv=0, 2
 iv = 0
 do ie=0, 2

   if (iv .eq. 0) then
     plot = xylineplot(wks,t+0.5/12.,var(iv*2+ie,iz,:),"D",0,"T",0.5,"C",1)
   else
     plot = xylineplot(wks,t,0,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t,var(iv*2+ie,:,:),"D",0,"T",0.5,"C",(/1,4,2,6/))
   end if

   set_size(plot,xf(ie),yf(ie),wi,he(iv))
   axis_range(plot,0,12,ymin(iv),ymax(iv),False,False)
   axis_str(plot,mstr(iv,ie),xstr(iv),ystr(iv*2+ie))
   setvalues plot
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 1
    "tmXBTickStartF"    : 0
    "tmXBMinorPerMajor" : 11
   end setvalues

   gsres             = True
   gsres@gsFillColor = colind

   dummy = new(nt/2,graphic)
   tag = 0
   cnt = 0
   do n=0, nt-1
     if ( .not. ismissing(var(ie,0,n)) ) then
       if ( var(ie,0,n) .ge. zb .and. var(ie,0,n) .lt. zt ) then
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

   draw(plot)

  end do
;  end do


 frame(wks)

end



