load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 dd = 2
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
; nave = 12*3
 nave = 12
 n1s = ispan(0,1151,nave)
 nkh = 32
 khi = fspan(0,4,nkh+1)
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_khz-"+var1+"-"+dd+".ps",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/fcoef_kl/fft_"+var1+"_kl__mls001-"+dd+".nc"
 f1 = addfile(fname1,"r")
 nx = 432
 ny = 302
 nz = 44
 nt = dimsizes(n1s)
 dx = 1./10.8
 dy = 1./10.8
 dt = 1./(12.*24.)
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 do k=0, nz-1
 fname1 = fdir+"/x6_5min/fcoef_kl/fft_"+var1+"_kl__mls"+sprinti("%3.3i",k)+"-"+dd+".nc"
 f1 = addfile(fname1,"r")
 do nn=0, nt-1
   v1 = dim_avg_n( \
         ( f1->$("fc_"+var1+"_r")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2 + \
         ( f1->$("fc_"+var1+"_i")$(ispan(n1s(nn),n1s(nn)+nave-1,1),:,:) )^2, 0 )
   if (nn .eq. 0 .and. k .eq. 0) then
     nd1 = dimsizes(v1)
     vs = new((/nt,nz,nd1(0),nd1(1)/),float,1.e32)
   end if
   vs(nn,k,:,:) = v1(:,:)
 end do
 end do
 vs(:,:,:,1:(nx-1)/2) = vs(:,:,:,1:(nx-1)/2)*2.
 vs = vs*(nx*dx)*(ny*dy)
 vs(:,:,:,0) = 1.e32

;--------------------------------------- process 0
 nv = dimsizes(vs)
 ndim = dimsizes(nv)

 nk = nx/2+1
 nl = ny
 kwn = ispan(0,nk-1,1)/tofloat(nx*dx)
 lwn = ispan(0,nl-1,1)/tofloat(ny*dy)
 lwn((nl+1)/2:nl-1) = -lwn(nl/2:1)

 khs = khi(0:nkh-1) + 0.5*(khi(1)-khi(0))
 t = n1s*dt + 4.

 z = ispan(0,nz-1,1)

 kh = sqrt( conform_dims((/nl,nk/),kwn^2,1) + conform_dims((/nl,nk/),lwn^2,0) )

 vsh = new((/nt,nz,nkh/),float)
 vsh = 0.
 do j=0, nl-1
 do i=1, nk-1
   do ii=0, nkh-1
     if ( kh(j,i) .gt. khi(ii) .and. kh(j,i) .le. khi(ii+1) ) then
       vsh(:,:,ii) = vsh(:,:,ii) + vs(:,:,j,i)
     end if
   end do
 end do
 end do

 vsh = vsh/(tofloat(nx*dx)*tofloat(ny*dy))/(khi(1)-khi(0))

rho = (/1.309292, \
1.298076, 1.28686, 1.272685, 1.254641, 1.231688, 1.202937, 1.169908, 1.135343, 1.100479, 1.065568, 1.023193, 0.9663156, 0.9037974, 0.8440109, 0.7869166, 0.7325803, 0.6811302, 0.6325406, 0.5867092, 0.5435455, 0.5029063, 0.4646725, 0.4287073, 0.3948925, 0.3631124, 0.3332551, 0.3052172, 0.2789088, 0.254227, 0.2310518, 0.2092237, 0.188495, 0.1687621, 0.1498067, 0.1315297, 0.1142049, 0.09782779, 0.08257273, 0.06877332, 0.05656385, 0.04588582, 0.03675439, 0.02907875/)
;, 0.02265134, 0.01740417, 0.01312034, 0.009641304, 0.006842336, \
;0.004043368/)

fpname = "$DATD/gwjet_wrf/x1/prt_d300km/prt_W__x1__08_00_00__12_00_00__6hrly.nc"
fp = addfile(fpname,"r")
 p = fp->h(:nz-1)*990. + 10.
;p = (/1000., \
;992.0709, 982.1647, 969.286, 953.1859, 933.1161, 908.5852, 880.9703, 852.6151, 824.376, 796.1371, 761.8624, 716.465, 667.5033, 621.3027, 577.7341, 536.6725, 497.9987, 461.5981, 427.3604, 395.1801, 364.9559, 336.5904, 309.991, 285.068, 261.7359, 239.9129, 219.5207, 200.4838, 182.7307, 166.1923, 150.8028, 136.4993, 123.2217, 110.9125, 99.51643, 88.98112, 79.25641, 70.29414, 62.04844, 54.47538, 47.53281, 41.18027, 35.37865, 30.08991, 25.27643, 20.89989, 16.91896, 13.28294, \
;9.64692/)

 do k=0, nz-1
   vsh(:,k,:) = vsh(:,k,:)*rho(k)
 end do

 tmp3 = vsh
 do i=1, nkh-2
   vsh(:,:,i) = (tmp3(:,:,i-1)+2.*tmp3(:,:,i)+tmp3(:,:,i+1))*0.25
 end do

 vsh = where(vsh .eq. 0., 1.e-30, vsh)
print(max(vsh))

;--------------------------------------- cn level
 if (var1 .eq. "W") then
   cnlev = fspan(-9.0,-4.0,11)*1.
   k9 = 3
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-3.5,1.5,11)*1.
   k9 = 3
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "k:B:h:N: [cyc deg:S:-1:N:]"
ystr = "p* [hPa]"
lbopt = "log:B:10:N:(PSD[m:S:2:N: s:S:-2:N: / cyc deg:S:-1:N:])"

 il = 0

 do nn=0, nt-1

   tday = toint(t(nn)+1.e-4)
   tmp = (t(nn)+1.e-4 - tday)*24.
   thour = toint(tmp)
   tmin = toint((tmp - thour)*60.)
   mstr = "Day "+tday+"  "+sprinti("%2.2i",thour)+"::"+sprinti("%2.2i",tmin)

 contour1 = cnshadeplot(wks,khs,p,log10(vsh(nn,:,:)),cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6 ;1.5
  "pmLabelBarHeightF"        : 0.08 ;0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.02 ;0.03
  "lbTitlePosition"    : "Bottom"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.022 ;0.04
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,k9,999,-999,False,True)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
; setvalues contour1
;  "tmYLMode" : "Manual"
;  "tmYLTickSpacingF"  : 1
;  "tmYLMinorPerMajor" : 8
; end setvalues

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
; dr_txt(wks,t(nn)+4,xf(0)-8,yf(il),25,6,0.)

 frame(wks)

 end do

end



