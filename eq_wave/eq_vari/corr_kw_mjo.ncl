load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 load "./ftn_extr-eqwave.ncl"

 expname = "era-int_f"
 dirname = "ERA-I-nr"
 uname = "u"
; expname = "jra55c_f"
; dirname = "JRA55C-nr"
; uname = "U_GRD_GDS4_HYBL"

 lat0 = 25.
 lat0avg = 15
; p0    = 100.  ; 54, 66, 80, 96
 p0s   = (/100,70,50/)  ; 54, 66, 80, 96
if (expname .eq. "era-int_f") then
 p0u   = 25   ; 22, 23, 24, 25
else
 p0u   = 36   ; 39, 38, 37, 36
end if
 ymax0 = 0.5

 xmin = 1979 ;1978
 xmax = 2016 ;2013 ;2014
 xmax2 = 2020

 wks = cr_wks("fig/corr_kw_mjo-"+expname(0),600,0)
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "$DATD/clim/"+dirname(0)+"/"+expname(0)+".u_ypt.nc"
   fnam1 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")

   fmjo = "omi.nc"

   f0 = addfile(fnam0,"r")

   nt = 1 + 407 + 36
if (expname .eq. "jra55c_f") then
  nt = 1 + 407 - 12
end if

 do k=0, dimsizes(p0s)-1
   p0 = p0s(k)
 do month=1, 12

   u = f0->$uname(0)$(month-1:nt-1:12,{p0u},{-lat0:lat0})

   rw_i = read_eqwave_rc0("K","f_z",(/fnam1,fnam2/),lat0,p0)
   lat1 = rw_i&lat
   ny = dimsizes(lat1)

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,ny/),float,1.e20)
   rw0(1:,:) = rw_i(:nt-2,:)
;   rw0(0,:) = rw0(1,:)
   rw0(0,:) = 1.e20

   nt1m = nt/12

   fz = new((/ny,nt1m/),float,1.e20)
   do n=0, nt1m-1
     fz(:,n) = rw0(month-1+n*12,:)
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

 if ( .not. isdefined("corr") ) then
   corr = new((/dimsizes(p0s),12/),float)
 end if

 corr(k,month-1) = esccr(-fzm,mjo({t}),0)

 end do
 end do

;--------------------------------------- cn level
 imrk = 1
 lcol = (/1,2,3/)
 ldsh = (/0,1,2/)*0
 ltck = (/1.5,0.8,0.8/)

 smon = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;--------------------------------------- contour
 wi = 30 ;60
 he = 20 ;30
 xf = 15
 yf = 10

 ystr = "Corr."
 xstr = "month"

 plot = xylineplot(wks,ispan(1,12,1),corr,"D",ldsh,"T",ltck,"C",lcol)
 xyaddmarker(plot,imrk,"S",4.0,"T",5.0,"C",lcol)

 xyaddline(plot,ispan(0,13,1),0,"D",0,"T",0.5,"C",colind)

; xyaddline(plot,yr_mjo,mjos,"D",0,"T",0.5,"C",2)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,0.5,12.5,-999,0.65,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : ispan(1,12,1)
 end setvalues

 draw(plot)

 frame(wks)

end

