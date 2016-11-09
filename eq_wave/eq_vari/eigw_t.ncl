load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = "era-int_f"
 dirname = "ERA-I-nr"

 lat0 = 25.
 lat0avg = 20
 p0    = 96.  ; 54, 66, 80, 96
 p0u   = 25   ; 22, 23, 24, 25
 ymax0 = 0.5
 month = 1

 xmin = 1978
 xmax = 2014
 xmax2 = 2020

 wks = cr_wks("fig/eigw_t-m"+month+"-p"+numeric2int(p0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "$DATD/clim/"+dirname(0)+"/"+expname(0)+".u_ypt.nc"
   fnam1 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")

   fmjo = "omi.nc"

   f0 = addfile(fnam0,"r")
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nt = 1 + 407
   nw = 8

   u = f0->u(month-1:nt-1:12,{p0u},{-lat0:lat0})

   rw_i = new((/nt-1,nw,ny/),float)
   rw_i(:,:,:) = f1->f_z(:,:nw-1,{p0},{-lat0:lat0})
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     f2->f_z(:,8:9,{p0},{-lat0:lat0})
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(:,2,{p0},{-lat0:lat0}) + \
                      f2->f_z(:,3,{p0},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(:,13,{p0},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(:,16:17,{p0},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,nw,ny/),float,1.e20)
   rw0!2 = "lat"
   rw0&lat = lat1
   rw0(1:,:,:) = rw_i(:,:,:)
;   rw0(0,:,:) = rw0(1,:,:)
   rw0(0,:,:) = 1.e20

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   nt1m = nt/12

   fz = new((/ny,nt1m/),float,1.e20)
   fz!0 = "lat"
   fz&lat = lat1
   do n=0, nt1m-1
     fz(:,n) = (/ rw(2,month-1+n*12,0,:) /)
   end do

   nt1m = nt/12

   t = 1979 + ispan(0,nt1m-1,1)

   fm = addfile(fmjo,"r")
   mjoi1 = fm->coef1
   mjoi2 = fm->coef2
   tmjo = fm->t
   mon_mjo = fm->mon
   yr0 = fm->year(0)

;--------------------------------------- process 0
 um = tofloat( dim_avg(u(:,{-lat0avg:lat0avg})) )
 fzm = dim_avg_n(fz({-lat0avg:lat0avg},:),0)

 um = where(ismissing(fzm), 1.e20, um)
 um@_FillValue = 1.e20

 ums = (um-min(um))/(max(um)-min(um))*abs(max(fzm)-min(fzm)) + abs(min(-fzm))

 print("um:  "+min(um)+" ~ "+max(um))
 print("corr:  "+esccr(-fzm,um,0))

 mjoi = sqrt(mjoi1^2 + mjoi2^2)

 yr_mjo = ispan(yr0,2019,1)
 mjo = new(2020-yr0,float,1.e20)
 mjo = 0.
 tag = 0
 yr = 0
 do n=0, dimsizes(tmjo)-1
   if (mon_mjo(n) .eq. month) then
     tag = tag + 1
     mjo(yr) = mjo(yr) + mjoi(n)
   else
     if (tag .ne. 0) then
       mjo(yr) = mjo(yr)/tofloat(tag)
       yr = yr + 1
       tag = 0
     end if
   end if
 end do
 if (mon_mjo(dimsizes(tmjo)-1) .eq. month) then
   mjo(yr) = mjo(yr)/tofloat(tag)
   yr = yr + 1
   tag = 0
 end if
 mjo = where(mjo .eq. 0., 1.e20, mjo)
 mjo!0 = "t"
 mjo&t = yr_mjo
 mjo({t}) = where(.not. ismissing(fzm), mjo({t}), 1.e20)

; mjo = runave(mjo,91,0)
 print("mjo:  "+min(mjo)+" ~ "+max(mjo))
 print("corr:  "+esccr(-fzm,mjo({t}),0))
; mjos = (mjo-min(mjo))/(max(mjo)-min(mjo))*abs(max(fzm)-min(fzm)) + abs(min(-fzm))

;--------------------------------------- cn level
 imrk = 1
 lcol = (/1,2,3/)

 smon = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 15
 yf = 10
 wi2 = tofloat(wi)/tofloat(xmax-xmin)*tofloat(xmax2-xmin)

 ystr = (/"|F:B:z:N:| (mPa)","OMI","U (m s:S:-1:N:)"/)
 xstr = "year"

 plot = xylineplot(wks,t,-fzm,"D",0,"T",1,"C",lcol(0))
 xyaddmarker(plot,imrk,"S",2.5,"T",5.0,"C",lcol(0))

; xyaddline(plot,t,ums,"D",0,"T",0.5,"C",lcol(2))

; xyaddline(plot,yr_mjo,mjos,"D",0,"T",0.5,"C",2)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,xmin,xmax,-999,999,False,False)
 axis_str(plot,"",xstr,ystr(0))
 setvalues plot
  "tmXBMaxTicks" : 8
  "tmYROn"       : False
  "tmYRBorderOn" : False
 end setvalues

 plot2 = xylineplot(wks,yr_mjo,mjo({t}),"D",0,"T",0.8,"C",lcol(1))
 xyaddmarker(plot2,imrk,"S",2.0,"T",5.0,"C",lcol(1))

 set_size(plot2,xf,yf,wi,he)
 axis_range(plot2,xmin,xmax,-999,999,False,False)
 axis_str(plot2,smon(month),"",ystr(1))
 setvalues plot2
  "tmXBOn"       : False
  "tmXTOn"       : False
  "tmYLOn"       : False
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tmYUseLeft"   : False
  "tmBorderLineColor"  : lcol(1)
  "tmYRMajorLineColor" : lcol(1)
  "tmYRMinorLineColor" : lcol(1)
  "tmYRLabelsOn"     : True
  "tmYRLabelFontColor" : lcol(1)
  "tiYAxisSide"      : "Right"
  "tiYAxisAngleF"    : 270
  "tiYAxisFontColor" : lcol(1)
 end setvalues

 draw(plot2)
 draw(plot)

 plot3 = xylineplot(wks,t,um,"D",0,"T",0.6,"C",lcol(2))
 xyaddmarker(plot3,imrk,"S",1.5,"T",5.0,"C",lcol(2))

 set_size(plot3,xf,yf,wi2,he)
 axis_range(plot3,xmin,xmax2,-999,999,False,False)
 axis_str(plot3,"","",ystr(2))
 setvalues plot3
  "tmXBOn"       : False
  "tmXTOn"       : False
  "tmYLOn"       : False
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tmYUseLeft"   : False
  "tmBorderLineColor"  : lcol(2)
  "tmYRMajorLineColor" : lcol(2)
  "tmYRMinorLineColor" : lcol(2)
  "tmYRLabelsOn"     : True
  "tmYRLabelFontColor" : lcol(2)
  "tiYAxisSide"      : "Right"
  "tiYAxisAngleF"    : 270
  "tiYAxisFontColor" : lcol(2)
 end setvalues

 draw(plot3)

 frame(wks)

end

