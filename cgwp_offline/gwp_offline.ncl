load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
; yy = 2005
; mm = "11"
 yy = 1993
 mm = "05"

 wks = cr_wks("fig/gwp_offline_"+yy+mm+"_CC11",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 fname0 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+".nc"
 fname1 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+"_rescale_1.nc"
 fname2 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+"_rescale_10.nc"

 f0 = addfile(fname0,"r")
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 zf = f0->z_f * 1.e-3
 zd = f0->z_d * 1.e-3

 mf = dim_avg( (/ f0->mflx_east, f1->mflx_east, f2->mflx_east, \
                  -f0->mflx_west, -f1->mflx_west, -f2->mflx_west /) )*1.e3

 fd = dim_avg( (/ f0->drag_u, f1->drag_u, f2->drag_u /) )*86400.
 
 do n=0, 2
   print("  E-MF [mPa]    W-MF [mPa]      z [km]")
   print(sprintf(" %10.6f",mf(n,:))+"    "+sprintf("%10.6f",mf(n+3,:))+"    "+sprintf("%10.6f",zf))
   print("Drag [m/s/day]    z [km]")
   print(sprintf(" %10.6f",fd(n,:))+"    "+sprintf("%10.6f",zd))
 end do

 mf((/2,5/),:) = mf((/2,5/),:)*0.1
 fd(2,:) = fd(2,:)*0.1

;--------------------------------------- cn level
 lcol = (/1,2,4,1,2,4/)
 ldsh = (/0,0,0,1,1,1/)

;--------------------------------------- ctr
 wi = 40
 he = 44
 xf = 20 + (/0,1/)*(wi+6)
 yf = 10

 mstr = "Nov "+yy
 if (mm .eq. "05") then
   mstr = "May "+yy
 end if
 xstr = "Zonal momentum flux [mPa]"
 ystr = "z [km]"

 plot1 = xylineplot(wks,mf,zf,"D",ldsh,"T",1,"C",lcol)
 xyaddline(plot1,0,zf,"D",0,"T",0.5,"C",colind)
; xyaddmarker(plot1,1,"S",3.0,"T",0.6,"C",1)

 set_size(plot1,xf(0),yf,wi,he)
 axis_range(plot1,-1.2,1.2,16,60,False,False)
 axis_str(plot1,mstr,xstr,ystr)

 draw(plot1)

 xstr = "Zonal drag [m s:S:-1:N: day:S:-1:N:]"
 ystr = "" ; "z [km]"

 plot1 = xylineplot(wks,fd,zd,"D",0,"T",1,"C",lcol(:2))
 xyaddline(plot1,0,zf,"D",0,"T",0.5,"C",colind)
; xyaddmarker(plot1,1,"S",3.0,"T",0.6,"C",1)

 set_size(plot1,xf(1),yf,wi,he)
 axis_range(plot1,-0.7,0.7,16,60,False,False)
 axis_str(plot1,mstr,xstr,ystr)

 draw(plot1)

 lgtxt = (/"as used","1 mPa","10 mPa (x 0.1)"/)
 lg = lgplot(wks,lgtxt,"D",0,"T",0.5,"C",lcol(:2),2)
 set_size_lg(lg,xf(0)+0.1,yf(0)-0.1,wi*0.20,12)
 draw(lg)


 frame(wks)

end

