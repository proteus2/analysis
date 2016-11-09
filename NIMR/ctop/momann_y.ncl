load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"uj","ui"/)

 kz = (/28,50/)

 wks = cr_wks("fig/momann_y-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read base data
   fnam1 = "$DATD/L60CGW-t/mmxm.uan"+expname(0)+"a.pd.nc"
   fnam2 = "$DATD/L60CGW-t/mmxm.uan"+expname(0)+"a.pg.nc"

   f1  = addfile(fnam1,"r")
   mfe1 = f1->F_e_gwconv(:,kz,:)
   mfw1 = f1->F_w_gwconv(:,kz,:) *(-1.)
   lat = f1->latitude
print(f1->hybrid_ht_1(kz))

   f2  = addfile(fnam2,"r")
   mfe2 = f2->Fp_e_gwussp(:,kz,:)
   mfw2 = f2->Fp_w_gwussp(:,kz,:)

   mfe2 = mfe2+mfe1
   mfw2 = mfw2+mfw1

   fnam3 = "$DATD/L60CTL/mmxm.uan"+expname(1)+"a.pg.nc"

   f3 = addfile(fnam3,"r")
   mfe3 = f3->Fp_e_gwussp(:,kz,:)
   mfw3 = f3->Fp_w_gwussp(:,kz,:)

   mf_avg = dim_avg_n_Wrap( \
             (/ (/mfe1, mfw1/), \
                (/mfe2, mfw2/), \
                (/mfe3, mfw3/) /), 2 )*1.e3

;--------------------------------------- plot

 xf = 20 + (/0,1/)*70
 yf = 10 + (/0,1/)*42
 wi = 60
 he = 30
 ystr = (/ (/"E", "W"/)+"-MF (mPa)", (/"",""/) /)
 llbl = (/"CGW","CGW+BGW","BGW (CTL)"/)
 ymax = (/1.5,0.75/)
 ldsh = (/0,0,1/)
 lcol = (/2,1,1/)

 do k=0, dimsizes(kz)-1
 do idr=0, 1

   plot = xylineplot(wks,lat,mf_avg(:,idr,k,:),"D",ldsh,"T",1.0,"C",lcol)

   set_size(plot,xf(k),yf(idr),wi,he)
   axis_range(plot,-90,90,0,ymax(k),False,False)
   axis_str(plot,"","lat (deg)",ystr(k,idr))

   draw(plot)

 end do
 end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",1.0,"C",lcol,2)
 set_size_lg(lg,xf,yf(0),wi*0.25,he*0.4)
 draw(lg)


 frame(wks)

end



