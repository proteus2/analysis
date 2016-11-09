load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 var1 = "W"
 lev = 8000
 nave = 12*12
 lat0 = (/58,42/)
 n1s = (/ ispan(576,1151,nave), ispan(0,575,nave) /)
 ismth = 1
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_z_kl-300",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir = "$DATD/gwjet_wrf/x6_5min/"+var1+"_z/prt_d300km/fcoef_kl"
 fname1 = fdir+"/fft_"+var1name+"_kl__z"+sprinti("%5.5i",lev)+"-1a.nc"
 fname2 = fdir+"/fft_"+var1name+"_kl__z"+sprinti("%5.5i",lev)+"-2a.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 nx = 432
 ny = 270
 ndim = dimsizes(n1s)
 nt = ndim(1)
 dx = 1./10.8 * 111.*cos(lat0/180.*3.141592)
 dy = 1./10.8 * 111.
; ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 do nn=0, nt-1
   v1 = dim_avg_n( \
         ( f1->$("fc_"+var1name+"_r")$(ispan(n1s(0,nn),n1s(0,nn)+nave-1,1),:,:) )^2 + \
         ( f1->$("fc_"+var1name+"_i")$(ispan(n1s(0,nn),n1s(0,nn)+nave-1,1),:,:) )^2, 0 )
   v2 = dim_avg_n( \
         ( f2->$("fc_"+var1name+"_r")$(ispan(n1s(1,nn),n1s(1,nn)+nave-1,1),:,:) )^2 + \
         ( f2->$("fc_"+var1name+"_i")$(ispan(n1s(1,nn),n1s(1,nn)+nave-1,1),:,:) )^2, 0 )
   if (nn .eq. 0) then
     nd1 = dimsizes(v1)
     vs = new((/2,nt,nd1(0),nd1(1)/),float,1.e32)
   end if
   vs(0,nn,:,:) = v1(:,:)
   vs(1,nn,:,:) = v2(:,:)
 end do
 vs(:,:,:,1:(nx-1)/2) = vs(:,:,:,1:(nx-1)/2)*2.
 vs(0,:,:,:) = vs(0,:,:,:)*(nx*dx(0))*(ny*dy)
 vs(1,:,:,:) = vs(1,:,:,:)*(nx*dx(1))*(ny*dy)
 vs(:,:,:,0) = 1.e32

;--------------------------------------- process 0
 nv = dimsizes(vs)
 ndim = dimsizes(nv)

 nk = nx/2+1
 nl = ny
 kwn = new((/2,nk/),float)
 kwn(0,:) = ispan(0,nk-1,1)/tofloat(nx*dx(0))
 kwn(1,:) = ispan(0,nk-1,1)/tofloat(nx*dx(1))
 lwn = ispan(-nl/2,(nl-1)/2,1)/tofloat(ny*dy)

 tmp = vs
 vs(:,:,nl/2:,:) = tmp(:,:,:(nl-1)/2,:)
 vs(:,:,:nl/2-1,:) = tmp(:,:,(nl+1)/2:,:)

 if (ismth .ne. 0) then
   vs = smth9(vs,0.5,0.25,False)
 end if

; xc = ispan(0,nx,1)*dx
; nvc = nv
; nvc(ndim-1) = nvc(ndim-1) + 1
; vc = new(nvc,float)
; vc(:,:,:nx-1) = v0(:,:,:nx-1)
; vc(:,:,nx   ) = v0(:,:,0    )

;--------------------------------------- cn level
 if (var1 .eq. "W") then
;   cnlev = (/ fspan(-2.5,0.0,11), fspan(-1.5,1.0,11) /)*1.
   cnlev = (/ fspan(-2.5,0.0,11), fspan(-2.5,0.0,11) /)*1. + 0.5
   k9 = 3
   l9 = (/3,3/)
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-1.0,4.0,11)*1.
   k9 = 1.5
   l9 = 2.5
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 25
 he = tofloat(wi)/k9*(l9*2)
 xf = 12 + (/0,1,2,3,4,5/)*(wi+6)
 yf = 10 + (/0,1/)*(he(0)+12)

 mstr  = new((/2,nt/),string)
 ystr  = new(nt,string)
 lbopt = new((/2,nt/),string)

 mstr = ""
 mstr(0,:) = (/"Day 6, 00::00-12::00", \
               "Day 6, 12::00-24::00", \
               "Day 7, 00::00-12::00", \
               "Day 7, 12::00-24::00"/)
 xstr = "k [10:S:-2:N: cyc km:S:-1:N:]"
 ystr = ""
 ystr(0) = "l [10:S:-2:N: cyc km:S:-1:N:]"
 lbopt = ""
 lbopt(1,nt-1) = "log:B:10:N:(PSD [m:S:2:N: s:S:-2:N: / cyc km:S:-1:N:])"
 txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

do id=0, 1
 if (id .eq. 1) then
   ystr(1) = ystr(0)
 end if
do il=0, nt-1
 contour1 = cnshadeplot(wks,kwn(id,:)*1.e2,lwn*1.e2,log10(vs(id,il,:,:)),cnlev(id,:),cnfil,lbopt(id,il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.04
  "pmLabelBarParallelPosF"   : -0.8;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.10
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.04
 end setvalues

 if ( id .eq. 0 .and. ( il .ge. 1 .and. il .le. 2 ) ) then
   dumm = new(3,graphic)
   plres = True
   plres@gsLineThicknessF = 0.8
   plres@gsLineDashPattern = 0
   xl1 = (/0.6 ,1.0 ,2.9,2.9,0.6 ,0.6 /)
   yl1 = (/-0.2,-0.2,-1.5 ,-2.4,-2.4,-0.2/)
   dumm(0) = gsn_add_polyline(wks,contour1,xl1,yl1,plres)
   xl2 = (/0.6 ,0.9 ,2.9,2.9,1.687,0.6  ,0.6 /)
   yl2 = (/0.02,0.02,-1 ,0.8,0.8  ,0.173,0.02/)
   dumm(1) = gsn_add_polyline(wks,contour1,xl2,yl2,plres)
   xl3 = (/0.  ,0.4 ,0.4,2.9,2.9,0. ,0. /)
   yl3 = (/-0.1,-0.1,0.2,1.8,2.4,2.4,2.4/)
   dumm(2) = gsn_add_polyline(wks,contour1,xl3,yl3,plres)
 end if

 if ( id .eq. 1 .and. ( il .ge. 2 .and. il .le. 3 ) ) then
   dumm = new(2,graphic)
   plres = True
   plres@gsLineThicknessF = 0.8
   plres@gsLineDashPattern = 0
;   xl4 = (/0. ,0.  ,0.3 ,1.8,0. /)
;   yl4 = (/2.8,0.02,0.02,2.8,2.8/)
   xl4 = (/0. ,0.  ,0.3 ,2.0,0. /)
   yl4 = (/2.9,0.02,0.02,2.9,2.9/)
   dumm(0) = gsn_add_polyline(wks,contour1,xl4,yl4,plres)
;   xl5 = (/0.2 ,0.4 ,1.7,2.9,2.9 ,0.2 ,0.2 /)
;   yl5 = (/-0.1,-0.1,1.2,1.2,-1.5,-1.5,-0.1/)
   xl5 = (/0.25 ,2.0  ,2.0 ,0.25,0.25 /)
   yl5 = (/-0.05,-0.05,-1.5,-1.5,-0.05/)
   dumm(1) = gsn_add_polyline(wks,contour1,xl5,yl5,plres)
 end if

 set_size(contour1,xf(il),yf(id),wi,he(id))
 axis_range(contour1,0,k9,-l9(id),l9(id),False,False)
 axis_str(contour1,mstr(id,il),xstr,ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 4
  "tmYLMode" : "Manual"
  "tmYLTickSpacingF"  : 1
  "tmYLMinorPerMajor" : 4
 end setvalues

; contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
; cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)
;
; set_size(contouru,xf(0),yf(il),wi,he(id))
; axis_range(contouru,y1,y2+1,18,z99,False,False)
; axis_str(contouru,"","","")
; setvalues contouru
;  "tmXBOn" : False
;  "tmXTOn" : False
;  "tmYLOn" : False
;  "tmYROn" : False
; end setvalues

 cnfinalize(contour1,wks)

 if ( .not. (id .eq. 1 .and. il .eq. 0) ) then
   draw(contour1)
 end if

; cnfinalize(contouru,wks)
;
; draw(contouru)
;
; dr_txt2(wks,txtstr(il),xf(0)+1,yf(il)+he(id)-2,22,6,0.,"CenterLeft",1)
;
 delete(contour1)
; delete(contouru)
;
end do
 delete(dumm)
; dr_txt(wks,subl(id),xf(0)-8,yf(id)-3,25,8,0.)
 dr_txt(wks,subl(id),xf(id)-10,yf(id)+1.5,25,8,0.)
end do

 dr_txt(wks,"1",xf(1)+7,yf(0)+he(0)-8,22,6,0.)
 dr_txt(wks,"2",xf(1)+wi-3,yf(0)+he(0)*0.5+5,22,6,0.)
 dr_txt(wks,"3",xf(1)+2,yf(0)+8,22,6,0.)
 dr_txt(wks,"4",xf(2)+2,yf(1)+3,22,6,0.)
 dr_txt(wks,"5",xf(2)+wi*2./3.-2,yf(1)+he(1)*0.5+2.5,22,6,0.)

 frame(wks)

end

