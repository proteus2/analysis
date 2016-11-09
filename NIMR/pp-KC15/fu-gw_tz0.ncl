load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"

 y0 = 1953
 y1 = 1980;+7
 y2 = 2005;+7
 lat0 = 5.

 wks = cr_wks("fig/fu-gw_tz0-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

   fnam1 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   z  = f1->z / 1.e3

   nz = dimsizes(z)
   nt = 648

   z1 = f1->z(0)
   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")
   zg = fg->hybrid_ht_1 / 1.e3

   gw  = new((/nt,nz/),float)
   gwc = new((/nt,nz/),float)
   gwb = new((/nt,nz/),float)
   gw (:,:) = dim_avg( fg->du_gw      (:,{z1:},{-lat0:lat0}) )*72. + \
              dim_avg( fg->dudt_gwconv(:,{z1:},{-lat0:lat0}) )*86400.
   gwc(:,:) = dim_avg( fg->dudt_gwconv(:,{z1:},{-lat0:lat0}) )*86400.
   gwb(:,:) = dim_avg( fg->du_gw      (:,{z1:},{-lat0:lat0}) )*72.


   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"

   f1 = addfile(fnam1,"r")
   zu = f1->hybrid_ht / 1.e3

   nzu = dimsizes(zu)

   u0 = transpose( dim_avg( f1->u(:nt-1,:,{-lat0:lat0}) ) )

   wf = (/ gw, gwb, gwc, gwc /)

   fu = new((/4,nz,nt/),float)
   do n=0, nt-1
     fu(:,:,n) = wf(:,n,:)
   end do
print(dimsizes(fu))

   fu!1 = "z"
   fu&z = z*1.e3

   t = y0 + (fspan(1,nt,nt)-0.5)/12.

;   tt = y0 + (fspan(1,nt-1,nt-1))/12.

;--------------------------------------- process 0
 fu = fu*30.

;--------------------------------------- cn level
 cscale = (/1,1,1,3/)

 cnlev = fspan(-10,10,21)*2.
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 4

 wi = ( 64 + (64./13.*2.) )*(26./12.)
 he = 15
 xf = 22 - (64./13.*2.) + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,nl-1,1)*(he+5)

 xstr  = new(nl,string)
 lbopt = new(nl,string)

 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"height (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/"CGW+BGW","BGW","CGW","CGW"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)
 do il=0, nl-1
   if (cscale(il) .ne. 1) then
     txtstr(il) = txtstr(il)+" (:F34:4:F22: "+cscale(il)+")"
   end if
 end do

ie = 0
do il=0, nl-1

 xstr(3) = "year"
 lbopt(3) = "(m s:S:-1:N: month:S:-1:N:)"

 contour1 = cnshadeplot(wks,t,z,fu(il,:,:)*cscale(il),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12 + 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.55
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05 + 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues

 set_size(contour1,xf(ie),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,"",xstr(il),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(ie),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,19,33,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(il),xf(ie)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)
; dr_txt(wks,subl(il),xf(0)-8,yf(il),25,8,0.)
 dr_txt(wks,"",xf(0)-8,yf(il),25,8,0.)

 delete(contour1)
 delete(contouru)

end do


 frame(wks)

end

