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

 load "res/data_from_Andrew"
 load "res/data_from_Andrew_HadGEM2-AC"

 wks = cr_wks("fig/gwp_offline_"+yy+mm+"_sum",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 fname0 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+".nc"
 fname1 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+"_rescale_0.5.nc"
 fname2 = "res/uanuj_ra.cgwp_offline."+yy+"."+mm+"_rescale_5.nc"

 f0 = addfile(fname0,"r")
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 zf = f0->z_f * 1.e-3
 zd = f0->z_d * 1.e-3

 mf = dim_avg( (/ f0->mflx_east, f1->mflx_east, f2->mflx_east, \
                  -f0->mflx_west, -f1->mflx_west, -f2->mflx_west /) )*1.e3

 fd = dim_avg( (/ f0->drag_u, f1->drag_u, f2->drag_u /) )*86400.
 
 if (yy+mm .eq. "199305") then
   mfa = (/ EqZm_Eastward_flux_ss9305, EqZm_Eastward_flux_1mpa9305,  \
            EqZm_Eastward_flux_10mpa9305, -EqZm_Westward_flux_ss9305,  \
            -EqZm_Westward_flux_1mpa9305, -EqZm_Westward_flux_10mpa9305 /)
   fda = (/ 0.2*EqZm_EW_accel_ss9305, 0.2*EqZm_EW_accel_1mpa9305, EqZm_EW_accel_10mpa9305 /)
 else if (yy+mm .eq. "200511") then
   mfa = (/ EqZm_Eastward_flux_ss0511, EqZm_Eastward_flux_1mpa0511,  \
            EqZm_Eastward_flux_10mpa0511, -EqZm_Westward_flux_ss0511,  \
            -EqZm_Westward_flux_1mpa0511, -EqZm_Westward_flux_10mpa0511 /)
   fda = (/ 0.2*EqZm_EW_accel_ss0511, 0.2*EqZm_EW_accel_1mpa0511, EqZm_EW_accel_10mpa0511 /)
 end if
 end if
 mfa((/0,3/),:) = mfa((/0,3/),:)*0.6
 fda(0,:) = fda(0,:)*0.6
 mfa(1:2,:) = mfa(1:2,:)*0.5
 mfa(4:5,:) = mfa(4:5,:)*0.5
 fda(1:2,:) = fda(1:2,:)*0.5
 if (yy+mm .eq. "199305") then
   mfa((/0,3/),:) = (/ EqZm_Eastward_flux_9305, -EqZm_Westward_flux_9305 /)
   fda(0,:) = EqZm_EW_accel_9305
 else if (yy+mm .eq. "200511") then
   mfa((/0,3/),:) = (/ EqZm_Eastward_flux_0511, -EqZm_Westward_flux_0511/)
   fda(0,:) = EqZm_EW_accel_0511
 end if
 end if
 
 dima = dimsizes(fda)
 mf = mf + mfa
 fd = fd + fda(:,:dima(1)-2)

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
 axis_range(plot1,-1.5,1.5,16,60,False,False)
 axis_str(plot1,mstr,xstr,ystr)

 draw(plot1)

 xstr = "Zonal drag [m s:S:-1:N: day:S:-1:N:]"
 ystr = "" ; "z [km]"

 plot1 = xylineplot(wks,fd,zd,"D",0,"T",1,"C",lcol(:2))
 xyaddline(plot1,0,zf,"D",0,"T",0.5,"C",colind)
; xyaddmarker(plot1,1,"S",3.0,"T",0.6,"C",1)

 set_size(plot1,xf(1),yf,wi,he)
 axis_range(plot1,-2.0,2.0,16,60,False,False)
 axis_str(plot1,mstr,xstr,ystr)

 draw(plot1)

 lgtxt = (/"as in model","1 mPa","10 mPa (x 0.1)"/)
 lg = lgplot(wks,lgtxt,"D",0,"T",0.5,"C",lcol(:2),2)
 set_size_lg(lg,xf(0)+0.1,yf(0)-0.1,wi*0.20,12)
 draw(lg)

 lgtxt = (/"as in model","1 mPa","10 mPa (x 0.1)"/)
 lg = lgplot(wks,lgtxt,"D",0,"T",0.5,"C",lcol(:2),2)
 set_size_lg(lg,xf(1)+0.1,yf(0)-0.1,wi*0.20,12)
 draw(lg)

 frame(wks)

end

