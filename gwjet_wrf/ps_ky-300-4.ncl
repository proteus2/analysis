load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 lev = 300
 nave = 108
 n1s = ispan(720,1151,nave)
 var1name = "prt_"+var1+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_ky-300"+var1name+"-4",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/fcoef_k/fft_"+var1name+"_k.nc"
 f1 = addfile(fname1,"r")
 nx = 432
 ny = 648+2
 nt = dimsizes(n1s)
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 do nn=0, nt-1
   v1 = dim_avg_n( \
         ( f1->$("fc_"+var1name+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2 + \
         ( f1->$("fc_"+var1name+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2, 0 )
   if (nn .eq. 0) then
     nd1 = dimsizes(v1)
     vs = new((/nt,nd1(0),nd1(1)/),float,1.e32)
   end if
   vs(nn,:,:) = v1(:,:)
 end do
 vs(:,:,1:(nx-1)/2) = vs(:,:,1:(nx-1)/2)*2.
 vs = vs*(nx*dx)
 vs(:,:,0) = 1.e32

;--------------------------------------- process 0
 nv = dimsizes(vs)
 ndim = dimsizes(nv)

 nk = nx/2+1
 kwn = ispan(0,nk-1,1)/tofloat(nx*dx)
 y = ys + ispan(0,ny-1,1)*dy

; xc = ispan(0,nx,1)*dx
; nvc = nv
; nvc(ndim-1) = nvc(ndim-1) + 1
; vc = new(nvc,float)
; vc(:,:,:nx-1) = v0(:,:,:nx-1)
; vc(:,:,nx   ) = v0(:,:,0    )

;--------------------------------------- cn level
 if (var1 .eq. "W") then
   cnlev = fspan(-7.5,-2.5,11)*1.
   k9 = 3
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-4,1,11)*1.
   k9 = 1.5
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

 mstr  = new(nt,string)
 xstr  = new(nt,string)
 ystr  = new(nt,string)
 lbopt = new(nt,string)

 mstr(:) = ""
 xstr(:) = ""
 xstr(:) = "k [cyc deg:S:-1:N:]"
 ystr(:) = ""
 ystr(ispan(0,nt-1,nt)) = "lat [deg]"
 lbopt(:) = ""
 lbopt(nt-1) = "log:B:10:N:(PSD[m:S:2:N: s:S:-2:N: / cyc deg:S:-1:N:])"
 txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

do il=0, nt-1

 contour1 = cnshadeplot(wks,kwn,y,log10(vs(il,:,:)),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,k9,20,80,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 0.5
  "tmXBMinorPerMajor" : 4
 end setvalues

; contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
; cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)
;
; set_size(contouru,xf(0),yf(il),wi,he)
; axis_range(contouru,y1,y2+1,18,z99,False,False)
; axis_str(contouru,"","","")
; setvalues contouru
;  "tmXBOn" : False
;  "tmXTOn" : False
;  "tmYLOn" : False
;  "tmYROn" : False
; end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

; cnfinalize(contouru,wks)
;
; draw(contouru)
;
; dr_txt2(wks,txtstr(il),xf(0)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)
;
 delete(contour1)
; delete(contouru)
;
; dr_txt(wks,subl(il),xf(0)-8,yf(il),25,6,0.)

end do

 frame(wks)

end



