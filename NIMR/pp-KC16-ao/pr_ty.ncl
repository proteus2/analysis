load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = (/"tc"/)
 expdir  = (/"AOL60CGW"/)

 lat0 = 20.
 ymax0 = 1.0    ; 1.0    ; 1.1
 y1 = 1951 ;2058 ;1960
 y2 = 2097 ;2097 ;1999
 y0 = 1951

 y1a = (/1955,2005,2055/)
 yap = 40

 wks = cr_wks("fig/pr_ty-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- read data
   fnam3 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"

   f3 = addfile(fnam3,"r")
   lat1 = f3->latitude({-lat0:lat0})

   nt = (y2 - y1 + 1)*12
   n1 = (y1 - y0)*12
   n2 = n1 + nt - 1

   ny = dimsizes(lat1)

   rw = new((/nt,ny/),float)
   rw = f3->precip(n1:n2,0,{-lat0:lat0})
   rw = rw*86400.  ; [mm/d]

   w = new((/ny,nt/),float)
   do n=0, nt-1
     w(:,n) = rw(n,:)
   end do

;   um = f3->u(n1:n2,{4.2e3:z0*1.e3},{-lat0:lat0})

   t = ispan(y1,y2,1)

   wan = new((/ny,(nt/12)/),float)
   do n=0, nt/12-1
     wan(:,n) = dim_avg( w(:,n*12:n*12+11) )
   end do
   wan!0 = "lat"
   wan&lat = lat1
   wan!1 = "time"
   wan&time = t

   wdec = dim_avg( (/ wan(:,{y1a(0):y1a(0)+yap}),  \
                      wan(:,{y1a(1):y1a(1)+yap}),  \
                      wan(:,{y1a(2):y1a(2)+yap}) /) )

;--------------------------------------- process 0
; do j=0, ny-1
;   wan(j,:) = wan(j,:)*cos(lat1(j)*3.141592/180.)
; end do

; wan = runave(wan,5,0)

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
 cnlev = (/4,5,6,7,8,9,10,11,12,13,14,15,16/)*0.5
 cnfil = (/5,6,7,8,10,12,14,15,16,17,18,19,20,21/)

 ldsh = (/0,14,16/)
 ltck = (/0.8,0.8,0.8/)
 lcol = (/colind,1,1/)

;--------------------------------------- contour
 wname = (/ (/"E-CGW","E-BGW","E-GW"/), (/"W-CGW","W-BGW","W-GW"/) /)
 wname = ""
 wi = 30
 wi2 = 15
 he = 22
 xf = 12 + (/0,1,2,3/)*36
 yf = 10 + (/0,1,2,3/)*28 + (/0,0,0,0/)*12

 mstr = new((/2,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = "" ;(/"Precip.","F:B:z:N:"/)
 xstr = (/"year"/)
 ystr = (/"latitude [deg]",""/)
 txtstr = wname
 lbopt = new((/2,2/),string)
 lbopt(:,:) = ""
 lbopt(0,0) = "[mm d:S:-1:N:]"
 ymax = ymax0*(/-1,1/)
 xstr2 = (/"Prec. "+lbopt(0,0)/)

do ifig=0, 0

do iv=0, 0

 contour1 = cnshadeplot(wks,t,lat1,wan,cnlev,cnfil,lbopt(iv,ifig))
; cnaddline(contour1,wks,t,lat1,uss,ispan(-60,60,6)/tofloat(ifig+1),1,"D",0,"T",0.75,"C",colind)
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.10
;  "pmLabelBarParallelPosF"   : 1.10
  "pmLabelBarKeepAspect"     : False
;  "pmLabelBarWidthF"         : 1.4
  "pmLabelBarHeightF"        : 0.12
;  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.025
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.01
  "lbTitleFontHeightF" : 0.03
 end setvalues

 plot1 = xylineplot(wks,wdec,lat1,"D",ldsh,"T",ltck,"C",lcol)

 set_size(contour1,xf(ifig),yf(iv),wi,he)
 axis_range(contour1,y1/5*5,y2+1,-lat0,lat0,False,False)
 axis_str(contour1,mstr(iv,ifig),xstr(iv),ystr(ifig))

 cnfinalize(contour1,wks)

 set_size(plot1,xf(ifig)+wi+6,yf(iv),wi2,he)
 axis_range(plot1,-999,999,-lat0,lat0,False,False)
 axis_str(plot1,"",xstr2(iv),"")
 getvalues plot1
  "xyCoordDataSpec" : spec
 end getvalues
 do i=0, dimsizes(spec)-1
   setvalues spec(i)
    "xyLineDashSegLenF" : 0.03
   end setvalues
 end do

 draw(contour1)
 draw(plot1)

 dr_txt2(wks,txtstr(ifig,iv),xf(ifig)+0.6,yf(iv)+he-2,22,5,0.,"CenterLeft",1)

end do  ; iv

end do  ; ifig

 lg = lgplot(wks,y1a+"-"+(y1a+yap-1),"D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(0)+wi+6,yf(0)+he+7,wi2*0.5,9)
 setvalues lg
  "lgLineDashSegLenF" : 0.025
 end setvalues
 draw(lg)

 frame(wks)

end

