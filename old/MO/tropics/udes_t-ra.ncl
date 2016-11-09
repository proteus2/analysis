load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 lat0 = 5.
 kz0 = (/5,9,13,15,19/)

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 25*12

 wks = cr_wks("udes_t",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,kz0,{-lat0:lat0}) )
   p  = f3->levels(kz0)
print(p)

   nt3 = nmon0
   nz = dimsizes(p)

   ur = new((/nz,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t3 = fspan(1,nt3,nt3)

   uss = new((/nz,12/),float)
   do n=0, 11
     uss(:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
   end do
   umean = dim_avg(uss)

   udes0 = new(dimsizes(ur),float)
   do n=0, nt3-1
     udes0(:,n) = ur(:,n) - uss(:,mod(n,12))
   end do

   ifilt0 = nt3/24-1

   fftu0 = ezfftf(udes0)
   fftu0(:,:,ifilt0+1:) = 0.
   ulow0 = ezfftb(fftu0,0.)

   do k=0, nz-1
     udes0(k,:) = udes0(k,:) + umean(k)
     ulow0(k,:) = ulow0(k,:) + umean(k)
   end do

;--------------------------------------- process 0

;--------------------------------------- plot
 xmax = (/144,144,nt3/)
 ymin = (/-25,-50,-50,-50,-50/)
 ymax = (/25,40,40,40,40/)
 wi = 72.*((/1.,1.,int2flt(nt3)/xmax(0)/))
 he = 20
 xf = 12 + (/0,1,0/)*78
 yf = 10 + (/0,1,2,3,4/)*26

 mstr = new((/nz,3/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC","MERRA"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "month"
 ystr = (/"U (m s:S:-1:N:)","","U (m s:S:-1:N:)"/)


 do k=0, nz-1
 do iv=2, 2

   plot = xylineplot(wks,t3,0.,"D",0,"T",0.5,"C",1)
   xyaddline(plot,t3,udes0(k,:),"D",0,"T",0.5,"C",1)
   if (k .gt. 1) then
     xyaddline(plot,t3,ulow0(k,:),"D",0,"T",0.5,"C",2)
   end if

   set_size(plot,xf(iv),yf(k),wi(iv),he)
   axis_range(plot,1,xmax(2),ymin(k),ymax(k),False,False)
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



