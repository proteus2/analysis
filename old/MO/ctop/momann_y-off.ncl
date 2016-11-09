load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "0"

 nmon = 4+144
 kz = 47
 kzl = 21

 wks = cr_wks("momann_y-off-"+expname,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read base data
   fnam1 = "../dat/cg_"+expname+".nc"

   f1  = addfile(fnam1,"r")
   mfe1 = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw1 = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   mfe10 = f1->unspecified_3(4:nmon-1,0,:,0)
   mfw10 = f1->unspecified_4(4:nmon-1,0,:,0) *(-1.)

   lat = f1->latitude

   fnam2 = "../dat/ng_"+expname+".nc"

   f2 = addfile(fnam2,"r")
   mfe2 = f2->field420(4:nmon-1,kz,:,0)
   mfw2 = f2->field422(4:nmon-1,kz,:,0)
   mfe20 = f2->field420(4:nmon-1,kzl,:,0)
   mfw20 = f2->field422(4:nmon-1,kzl,:,0)
   zl = int2flt(floattoint((f2->hybrid_ht(kzl)/1.e3)*10))/10.

   mf_avg = dim_avg_n_Wrap( \
             (/ (/ (/mfe20, mfw20/), \
                   (/mfe2 , mfw2 /) /), \
                (/ (/mfe10, mfw10/), \
                   (/mfe1 , mfw1 /) /) /), 3 )*1.e3

;--------------------------------------- plot

 xf = 20 + (/0,1/)*70
 yf = 10 + (/0,1/)*42
 wi = 60
 he = 30
 mstr = (/ (/"NGW","CGW (off-line)"/), (/"",""/) /)
 ystr = (/ (/"E", "W"/)+"-MF (mPa)", (/"",""/) /)
 llbl = (/ (/zl+" km","16 km"/), (/"cloud top","16 km"/) /)
 ldsh = (/1,0/)
 lcol = (/1,1/)

 do iw=0, 1
 do idr=0, 1

   plot = xylineplot(wks,lat,mf_avg(iw,:,idr,:),"D",ldsh,"T",1.0,"C",lcol)

   set_size(plot,xf(iw),yf(idr),wi,he)
   axis_range(plot,-90,90,0,3,False,False)
   axis_str(plot,mstr(idr,iw),"lat (deg)",ystr(iw,idr))

   draw(plot)

 end do
 end do

 do iw=0, 1
   lg = lgplot(wks,llbl(iw,:),"D",ldsh,"T",1.0,"C",lcol,2)
   set_size_lg(lg,xf(iw)+1,yf(0),wi*0.25,he*0.28)
   draw(lg)
 end do


 frame(wks)

end



