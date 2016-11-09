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
 nave = 12*9
 n1s = ispan(720,1151,nave)
 nkh = 160
 khi = fspan(0,8,nkh+1)
 var1name = "prt_"+var1+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_wl-300"+var1name,600,0)
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/fcoef_kl/fft_"+var1name+"_kl.nc"
 f1 = addfile(fname1,"r")
 nx = 432
 ny = 648
 nt = dimsizes(n1s)
 dx = 1./10.8
 dy = 1./10.8
 dt = 1./12.
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
 vs = vs*(nx*dx)*(ny*dy)
 vs(:,:,0) = 1.e32

;--------------------------------------- process 0
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

 vsh = vsh*log(10.)*conform(vsh,khs,1)

 vsh = where(vsh .eq. 0., 1.e-30, vsh)

 wl = where(khs .ne. 0., 1./khs, 999)

;--------------------------------------- cn level
 lcol = (/1,2,3,4/)

 if (var1 .eq. "W") then
   cnlev = fspan(-9.0,-4.0,11)*1.
   k9 = 7
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-3.5,1.5,11)*1.
   k9 = 3
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 32
 he = 20
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "horizontal wavelength [deg]"
ystr = "PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"

 il = 0

 plot = xylineplot(wks,wl,vsh*1.e5,"D",0,"T",1,"C",lcol)

 set_size(plot,xf(il),yf(il),wi,he)
 axis_range(plot,40,-999,-999,999,True,False)
 axis_str(plot,mstr(il),xstr(il),ystr(il))

 draw(plot)

; dr_txt(wks,subl(il),xf(0)-8,yf(il),25,6,0.)

 frame(wks)

end



