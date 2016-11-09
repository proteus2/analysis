load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fa","fd","fb","fc","0"/)

 nmon = 4+96;36;144
 kz = 47

 wks = cr_wks("momann_y2",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read base data
   fnam1 = "../dat/cg_"+expname(0)+".nc"
   f1  = addfile(fnam1,"r")
   mfe  = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw  = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   lat = f1->latitude
   fnam2 = "../dat/ng_"+expname(0)+".nc"
   f2 = addfile(fnam2,"r")
   mfe0 = f2->field420(4:nmon-1,kz,:,0)
   mfw0 = f2->field422(4:nmon-1,kz,:,0)
   mfe1 = mfe+mfe0
   mfw1 = mfw+mfw0

   fnam1 = "../dat/cg_"+expname(1)+".nc"
   f1  = addfile(fnam1,"r")
   mfe  = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw  = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   lat = f1->latitude
   fnam2 = "../dat/ng_"+expname(1)+".nc"
   f2 = addfile(fnam2,"r")
   mfe0 = f2->field420(4:nmon-1,kz,:,0)
   mfw0 = f2->field422(4:nmon-1,kz,:,0)
   mfe2 = mfe+mfe0
   mfw2 = mfw+mfw0

   fnam1 = "../dat/cg_"+expname(2)+".nc"
   f1  = addfile(fnam1,"r")
   mfe  = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw  = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   lat = f1->latitude
   fnam2 = "../dat/ng_"+expname(2)+".nc"
   f2 = addfile(fnam2,"r")
   mfe0 = f2->field420(4:nmon-1,kz,:,0)
   mfw0 = f2->field422(4:nmon-1,kz,:,0)
   mfe3 = mfe+mfe0
   mfw3 = mfw+mfw0

   fnam1 = "../dat/cg_"+expname(3)+".nc"
   f1  = addfile(fnam1,"r")
   mfe  = f1->unspecified_1(4:nmon-1,kz,:,0)
   mfw  = f1->unspecified_2(4:nmon-1,kz,:,0) *(-1.)
   lat = f1->latitude
   fnam2 = "../dat/ng_"+expname(3)+".nc"
   f2 = addfile(fnam2,"r")
   mfe0 = f2->field420(4:nmon-1,kz,:,0)
   mfw0 = f2->field422(4:nmon-1,kz,:,0)
   mfe4 = mfe+mfe0
   mfw4 = mfw+mfw0

   fnam3 = "../dat/ng_"+expname(4)+".nc"
   f3 = addfile(fnam3,"r")
   mfe5 = f3->field420(4:nmon-1,kz,:,0)
   mfw5 = f3->field422(4:nmon-1,kz,:,0)

   mf_avg = dim_avg_n_Wrap( \
             (/ (/mfe1, mfw1/), \
                (/mfe2, mfw2/), \
                (/mfe3, mfw3/), \
                (/mfe4, mfw4/), \
                (/mfe5, mfw5/) /), 2 )*1.e3

;--------------------------------------- plot

 xf = 20
 yf = 10 + (/0,1/)*42
 wi = 60
 he = 30
 ystr = (/"E", "W"/)+"-MF (mPa)"
 llbl = (/"const.","60:S:o:N:,10:S:o:N:","60:S:o:N:,5:S:o:N:","50:S:o:N:,10:S:o:N:","CTL"/)
 ldsh = (/0,0,0,0,1/)
 lcol = (/2,1,3,4,1/)

 do idr=0, 1

   plot = xylineplot(wks,lat,mf_avg(:,idr,:),"D",ldsh,"T",0.8,"C",lcol)

   set_size(plot,xf,yf(idr),wi,he)
   axis_range(plot,-90,90,0,3,False,False)
   axis_str(plot,"","lat (deg)",ystr(idr))

   draw(plot)

 end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",1.0,"C",lcol,2)
 set_size_lg(lg,xf,yf(0),wi*0.25,he*0.6)
 draw(lg)


 frame(wks)

end



