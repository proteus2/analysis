load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 0.8

begin

;--------------------------------------- parameters
 l_cal = 0
 dz = 250 

 wks = cr_wks("fig/ps_wlv-300prt_W-5",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- main body

; read -----------------------------------------------------
if (l_cal .ne. 0) then

do iw=0, 4  ; iw

 fdir = "$DATD/gwjet_wrf/x6_5min"
; fdname = "domain_x6.nc"
 if (iw .lt. 3) then
   ftail = "__x6__10_00_00__11_12_00__1hrly-w"+(iw+1)+".nc"
;   ns = ispan(144,431,1)
   ns = ispan(12,35,1)
 else
   ftail = "__x6__11_00_00__12_00_00__1hrly-w"+(iw+1)+".nc"
;   ns = ispan(0,287,1)
   ns = ispan(0,23,1)
 end if
 var = "W"
 lev = 8000
 lev0 = 0
 lev9 = 25000

 nz = (lev9-lev0)/dz + 1
 nt = dimsizes(ns)

 do k=0, nz-1
   z = lev0 + k*dz
   fname0 = fdir+"/"+var+"_z/prt_d300km/recon_xy/"+ \
                "prt_"+var+"__z"+sprinti("%5.5i",z)+ftail
   ff0 = addfile(fname0,"r")
   wf0 = ff0->$("prt_"+var)$(ns,:,:) * exp(-z/(2.*7.e3))
   if ( iw .eq. 0 .and. k .eq. 0 ) then
     ndim = dimsizes(wf0)
     ny = ndim(1)
     nx = ndim(2)
     wfz = new((/nz,nt,ny,nx/),float)
   end if
   wfz(k,:,:,:) = wf0
 end do

; process --------------------------------------------------
 kk = (lev - lev0)/dz

 wl = wavelet_default(wfz(:,0,0,0),0)
 na = dimsizes(wl@power)
 nl = na/nz

 if (iw .eq. 0) then
   psd = new((/5,nl/),float,1.e20)
 end if

 psd(iw,:) = 0.
 do n=0, nt-1
 do j=0, ny-1
 do i=0, nx-1
   wz = wfz(:,n,j,i)
   wl = wavelet_default(wz,0)
;   power = wl@power(ispan(kk,na-1,nz))
   psd(iw,:) = psd(iw,:) + wl@power(ispan(kk,na-1,nz))
 end do
 end do
 end do
 psd(iw,:) = psd(iw,:)/tofloat(nx*ny*nt)

 psd(iw,:) = psd(iw,:)/wl@scale

 sc = dz*wl@scale
 lz = dz*wl@period

 print(iw+1)
 print(lz+"  "+sc)
 print(sc+"  "+psd(iw,:))

end do  ; iw

else  ; l_cal

 nl = 23

 ws = asciiread("zl-s",(/nl,2/),"float")
 w1 = asciiread("zl-1",(/nl,2/),"float")
 w2 = asciiread("zl-2",(/nl,2/),"float")
 w3 = asciiread("zl-3",(/nl,2/),"float")
 w4 = asciiread("zl-4",(/nl,2/),"float")
 w5 = asciiread("zl-5",(/nl,2/),"float")

 lz = ws(:,0)
 sc = ws(:,1)
 psd = (/ w1(:,1), w2(:,1), w3(:,1), w4(:,1), w5(:,1) /)

end if  ; l_cal

 lz = lz*1.e-3
 do iw=0, 4
   psd(iw,:) = psd(iw,:)*log(10.)*(sc/dz)/lz
   psd(iw,:) = where(sc .le. 16.e3, psd(iw,:), 1.e20)
 end do

 psd@_FillValue = 1.e20

 sumsp = dim_sum(psd(:,1:))
 sumsp1 = psd*0.
 do ik=2, nl-2
   sumsp1(:,ik) = dim_sum(psd(:,ik:nl-1))/sumsp
 end do
 do iw=0, 4
   print(lz(maxind(psd(iw,:))))
   print(lz(minind(abs(sumsp1(iw,:)-0.1)))+"  "+lz(minind(abs(sumsp1(iw,:)-0.9))))
 end do

; plot setting ---------------------------------------------
 lcol = (/18,16,13,10,1/)

; plot -----------------------------------------------------
 wi = 60
 he = 30
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

 mstr = ""
 xstr = "vertical wavelength [km]"
 ystr = "PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"

 il = 0

 psd(3,:) = psd(3,:)*0.5

 plot = xylineplot(wks,lz,psd*1.e5,"D",0,"T",1,"C",lcol)

 set_size(plot,xf(il),yf(il),wi,he)
 axis_range(plot,20.01,1,-999,999,True,False)
 axis_str(plot,mstr(il),xstr(il),ystr(il))

 draw(plot)

 do iw=0, 4
   txt = "W"+(iw+1)
   if (iw .eq. 3) then
     txt = txt+" / 2"
   end if
   dr_txt2(wks,txt,xf(0)+wi-6.5,yf(0)+1+(iw*2.5),22,4,0,0,lcol(iw))
 end do

frame(wks)

end

