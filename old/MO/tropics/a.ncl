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
 kz  = (/9,12/)
 kz0 = (/7,9/)

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 13*12

 wks = cr_wks("fig/a-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,kz,{-lat0:lat0},0) )
   p1 = f1->p(kz)
print(p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,kz,{-lat0:lat0},0) )

   nz = dimsizes(p1)
   nt = nmon-4

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

   t = fspan(1,nt,nt)
   t3 = fspan(1,nt3,nt3)

   uss = new((/3,nz,12/),float)
   do n=0, 11
     uss(0:1,:,n) = dim_avg(u(0:1,:,ispan(n,nt-1,12)))
     uss(2  ,:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
   end do
   umean = dim_avg(uss)

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(0:1,:,mod(n,12))
   end do

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,n) = ur(:,n) - uss(2,:,mod(n,12))
   end do

   ifilt = nt/24-1
   ifilt0 = nt3/24-1

   fftu = ezfftf(udes)
   fftu(:,:,:,ifilt+1:) = 0.
   ulow = ezfftb(fftu,0.)

   fftu0 = ezfftf(udes0)
   fftu0(:,:,ifilt0+1:) = 0.
   ulow0 = ezfftb(fftu0,0.)

   do k=0, nz-1
     udes(0,k,:) = udes(0,k,:) + umean(0,k)
     udes(1,k,:) = udes(1,k,:) + umean(1,k)
     udes0 (k,:) = udes0 (k,:) + umean(2,k)
     ulow(0,k,:) = ulow(0,k,:) + umean(0,k)
     ulow(1,k,:) = ulow(1,k,:) + umean(1,k)
     ulow0 (k,:) = ulow0 (k,:) + umean(2,k)
   end do
print(udes0(0,:))

;--------------------------------------- process 0

;--------------------------------------- plot
 xmax = nt3
 ymin = (/-25,-50,-50,-50,-50/)
 ymax = (/25,40,40,40,40/)
 wi = 72.*((/1.,1.,int2flt(nt3)/xmax/))
 he = 20
 xf = 12 + (/0,0,0/)*78
 yf = 10 + (/0,1,2,3,4/)*30

 mstr = new((/nz,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "month"
 ystr = (/"U (m s:S:-1:N:)","",""/)


 do k=0, 0
 do iv=0, 2

   if (iv .le. 1) then
     plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t,udes(iv,k,:),"D",0,"T",0.5,"C",1)
     if (k .gt. 1) then
       xyaddline(plot,t,ulow(iv,k,:),"D",0,"T",0.5,"C",2)
     end if
   else
     plot = xylineplot(wks,t3,0.,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t3,udes0(k,:),"D",0,"T",0.5,"C",1)
     if (k .gt. 1) then
       xyaddline(plot,t3,ulow0(k,:),"D",0,"T",0.5,"C",2)
     end if
   end if

   set_size(plot,xf(iv),yf(iv),wi(iv),he)
   axis_range(plot,1,xmax,-40,40,False,False)
   axis_str(plot,mstr(k,iv),xstr(k),ystr(iv))
   setvalues plot
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 12
    "tmXBTickStartF"    : 0
    "tmXBMinorPerMajor" : 11
   end setvalues

   dr_txt(wks,p(k)+" hPa",xf(0)+8,yf(k)+3,22,8,0.)

   draw(plot)

  end do
  end do


 frame(wks)

end



