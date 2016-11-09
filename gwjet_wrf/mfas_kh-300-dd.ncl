load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 dd = 1
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 var2 = "U"
 var3 = "V"
 lev = 300
 nave = 12*3
 n1s = ispan(0,1151,nave)
 nkh = 160
 khi = fspan(0,8,nkh+1)
 var1name = "prt_"+var1+"_"+lev
 var2name = "prt_"+var2+"_"+lev
 var3name = "prt_"+var3+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/mfas_kh-300"+var1name+"-"+dd,600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/fcoef_kl/fft_"+var1name+"_kl-"+dd+".nc"
 fname2 = fdir+"/x1/prt_d300km/fcoef_kl/fft_"+var2name+"_kl-"+dd+".nc"
 fname3 = fdir+"/x1/prt_d300km/fcoef_kl/fft_"+var3name+"_kl-"+dd+".nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f3 = addfile(fname3,"r")
 nx = 432
 ny = 302
 nt = dimsizes(n1s)
 dx = 1./10.8
 dy = 1./10.8
 dt = 1./12.
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 do nn=0, nt-1
   v1 = dim_avg_n( sqrt( \
          ( f1->$("fc_"+var1name+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:)* \
            f2->$("fc_"+var2name+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) + \
            f1->$("fc_"+var1name+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:)* \
            f2->$("fc_"+var2name+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2 + \
          ( f1->$("fc_"+var1name+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:)* \
            f3->$("fc_"+var3name+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) + \
            f1->$("fc_"+var1name+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:)* \
            f3->$("fc_"+var3name+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2 \
        ), 0 )
   if (nn .eq. 0) then
     nd1 = dimsizes(v1)
     vs = new((/nt,nd1(0),nd1(1)/),float,1.e32)
   end if
   vs(nn,:,:) = v1(:,:)
 end do
 vs(:,:,1:(nx-1)/2) = vs(:,:,1:(nx-1)/2)*2.
 vs = vs*(nx*dx)*(ny*dy)
 vs(:,:,0) = 1.e32

;--------------------------------------- process 0
 t0 = 240.  ; at 300 hPa
 rho = lev*1.e2/287./t0
 vs = vs*rho
 vs(:,:,0) = 1.e32

 nv = dimsizes(vs)
 ndim = dimsizes(nv)

 nk = nx/2+1
 nl = ny
 kwn = ispan(0,nk-1,1)/tofloat(nx*dx)
 lwn = ispan(0,nl-1,1)/tofloat(ny*dy)
 lwn((nl+1)/2:nl-1) = -lwn(nl/2:1)

 khs = khi(0:nkh-1) + 0.5*(khi(1)-khi(0))
 t = n1s*dt

 kh = sqrt( conform_dims((/nl,nk/),kwn^2,1) + conform_dims((/nl,nk/),lwn^2,0) )

 vsh = new((/nt,nkh/),float)
 vsh = 0.
 do j=0, nl-1
 do i=1, nk-1
   do ii=0, nkh-1
     if ( kh(j,i) .gt. khi(ii) .and. kh(j,i) .le. khi(ii+1) ) then
       vsh(:,ii) = vsh(:,ii) + vs(:,j,i)
     end if
   end do
 end do
 end do

 vsh = vsh/(tofloat(nx*dx)*tofloat(ny*dy))/(khi(1)-khi(0))

 vsh = where(vsh .eq. 0., 1.e-30, vsh)
print(max(vsh))

;--------------------------------------- cn level
 cnlev = fspan(-7.5,-2.5,11)*1.
 k9 = 7
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 32
 he = 45
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

; mstr  = new(nt,string)
; xstr  = new(nt,string)
; ystr  = new(nt,string)
; lbopt = new(nt,string)
;
; mstr(:) = ""
; xstr(:) = ""
; xstr(nt/2:) = "k [cyc deg:S:-1:N:]"
; ystr(:) = ""
; ystr(ispan(0,nt-1,nt/2)) = "l [cyc deg:S:-1:N:]"
; lbopt(:) = ""
; lbopt(nt-1) = "log:B:10:N:(MSD[Pa / cyc:S:2:N: deg:S:-2:N:])"
; txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz"/)
; subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)
mstr = ""
xstr = "k:B:h:N: [cyc deg:S:-1:N:]"
ystr = "t [day]"
lbopt = "log:B:10:N:(MSD[Pa / cyc deg:S:-1:N:])"

 il = 0

 contour1 = cnshadeplot(wks,khs,t/24.+4.,log10(vsh),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.4 ;0.6 ;1.5
  "pmLabelBarHeightF"        : 0.08 ;0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.01 ;0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.015 ;0.04
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,k9,-999,999,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmYLMode" : "Manual"
  "tmYLTickSpacingF"  : 1
  "tmYLMinorPerMajor" : 8
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

 frame(wks)

end



