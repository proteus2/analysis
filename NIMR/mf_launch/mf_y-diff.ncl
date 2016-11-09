load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.1
load "$UTIL/ncl/UM/get_um_t.ncl"

begin

;--------------------------------------- parameters
 expname = (/"tc","tk"/)
 expdir  = "AOL60"+(/"CGW","CTL"/)

; z0 = (/8, 20/)
; ymax = (/0.3,0.4/)
 z0 = 17
 ymax = 0.2
 y1 = 1960 ;2058 ;1960
 y2 = 1999 ;2097 ;1999
 ydif = 98

 wks = cr_wks("fig/mf_y-z"+toint(z0(0))+"-diff-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read base data
 fnam1 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pd.nc"
 fnam2 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pg.nc"
 fnam3 = "$DATD/"+expdir(1)+"/mmxm.uan"+expname(1)+"a.pg.nc"

 f1 = addfile(fnam1,"r")
 lat = f1->latitude
 z = f1->F_e_gwconv&$f1->F_e_gwconv!1$ / 1.e3
 z!0 = "hgt"
 z&hgt = z
print(z({z0})+" km")
 y1r = get_um_t(y1  , f1->F_e_gwconv, 0, 360, True )
 y2r = get_um_t(y2+1, f1->F_e_gwconv, 0, 360, False)
 y1rf = y1r + ydif*360.
 y2rf = y2r + ydif*360.

 f1  = addfile(fnam1,"r")
 mfe1 = f1->F_e_gwconv({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
        f1->F_e_gwconv({y1r:y2r},{z0*1.e3:z0*1.e3},:)
 mfw1 = ( f1->F_w_gwconv({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
          f1->F_w_gwconv({y1r:y2r},{z0*1.e3:z0*1.e3},:) )*(-1.)

 f2  = addfile(fnam2,"r")
 mfe2 = f2->Fp_e_gwussp({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
        f2->Fp_e_gwussp({y1r:y2r},{z0*1.e3:z0*1.e3},:)
 mfw2 = f2->Fp_w_gwussp({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
        f2->Fp_w_gwussp({y1r:y2r},{z0*1.e3:z0*1.e3},:)

 mfe2 = mfe2+mfe1
 mfw2 = mfw2+mfw1

 f3 = addfile(fnam3,"r")
 mfe3 = f3->Fp_e_gwussp({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
        f3->Fp_e_gwussp({y1r:y2r},{z0*1.e3:z0*1.e3},:)
 mfw3 = f3->Fp_w_gwussp({y1rf:y2rf},{z0*1.e3:z0*1.e3},:) - \
        f3->Fp_w_gwussp({y1r:y2r},{z0*1.e3:z0*1.e3},:)

 mf_avg = dim_avg_n_Wrap( \
           (/ (/mfe1, mfw1/), \
              (/mfe2, mfw2/), \
              (/mfe3, mfw3/) /), 2 )*1.e3

 mf_avg(:,:,:,0) = mf_avg@_FillValue
 mf_avg(:,:,:,dimsizes(lat)-1) = mf_avg@_FillValue

;--------------------------------------- plot
 xf = 20 + (/0,1,2,3/)*58
 yf = 10 + (/0,1/)*36
 wi = 50
 he = 30
 ystr = new((/dimsizes(z0),2/),string,"")
 ystr(0,:) = (/"Eastward", "Westward"/)+"-MF (mPa)"
 xstr = (/"","lat (deg)"/)
 llbl = (/"CGW","CGW + BGW","BGW (CTL)"/)
 ldsh = (/0,0,1/)
 lcol = (/2,1,1/)

 do k=0, dimsizes(z0)-1
 do idr=0, 1

   plot = xylineplot(wks,lat,mf_avg(:,idr,k,:),"D",ldsh,"T",1.0,"C",lcol)
   xyaddline(plot,lat,0.,"D",0,"T",0.8,"C",colind)

   set_size(plot,xf(k),yf(idr),wi,he)
   axis_range(plot,-90,90,-ymax(k),ymax(k),False,False)
   axis_str(plot,"",xstr(idr),ystr(k,idr))

   draw(plot)

 end do
 end do

 lg = lgplot(wks,llbl,"D",ldsh,"T",1.0,"C",lcol,2)
 set_size_lg(lg,xf,yf(0),wi*0.20,he*0.36)
 draw(lg)


 frame(wks)

end



