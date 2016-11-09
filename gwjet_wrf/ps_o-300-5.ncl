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
 lev = 8000
 dphi = 15.
 var1name = "prt_"+var1
 
 wks = cr_wks("fig/ps_o-300"+var1name+"-5",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-w"+ispan(1,5,1)+".nc"

 f1 = addfile(fname1(0),"r")
 f2 = addfile(fname1(1),"r")
 f3 = addfile(fname1(2),"r")
 f4 = addfile(fname1(3),"r")
 f5 = addfile(fname1(4),"r")

; kh = f1->kh / 111.
; nkh = dimsizes(kh)
 ome = f1->ome / 1.
 no = dimsizes(ome)
 nphi = dimsizes(f1->dir)

 vs1w = dim_sum( f1->psd_phi_o(:,0:nphi/4) ) + dim_sum( f1->psd_phi_o(:,nphi/4*3+1:) )
 vs2w = dim_sum( f2->psd_phi_o(:,0:nphi/4) ) + dim_sum( f2->psd_phi_o(:,nphi/4*3+1:) )
 vs3w = dim_sum( f3->psd_phi_o(:,0:nphi/4) ) + dim_sum( f3->psd_phi_o(:,nphi/4*3+1:) )
 vs4w = dim_sum( f4->psd_phi_o(:,0:nphi/4) ) + dim_sum( f4->psd_phi_o(:,nphi/4*3+1:) )
 vs5w = dim_sum( f5->psd_phi_o(:,0:nphi/4) ) + dim_sum( f5->psd_phi_o(:,nphi/4*3+1:) )

 vs1 = dim_sum( f1->psd_phi_o(:,nphi/4+1:nphi/4*3) )
 vs2 = dim_sum( f2->psd_phi_o(:,nphi/4+1:nphi/4*3) )
 vs3 = dim_sum( f3->psd_phi_o(:,nphi/4+1:nphi/4*3) )
 vs4 = dim_sum( f4->psd_phi_o(:,nphi/4+1:nphi/4*3) )
 vs5 = dim_sum( f5->psd_phi_o(:,nphi/4+1:nphi/4*3) )

 vs1(0) = vs1(0) + vs1w(0)
 vs2(0) = vs2(0) + vs2w(0)
 vs3(0) = vs3(0) + vs3w(0)
 vs4(0) = vs4(0) + vs4w(0)
 vs5(0) = vs5(0) + vs5w(0)

 vs1a = array_append_record( vs1w(no-1:1), vs1, 0 )
 vs2a = array_append_record( vs2w(no-1:1), vs2, 0 )
 vs3a = array_append_record( vs3w(no-1:1), vs3, 0 )
 vs4a = array_append_record( vs4w(no-1:1), vs4, 0 )
 vs5a = array_append_record( vs5w(no-1:1), vs5, 0 )

 omea = array_append_record( (-1.)*ome(no-1:1), ome, 0 )

 vs = (/vs1a,vs2a,vs3a,vs4a,vs5a/) * 1. * dphi

;--------------------------------------- process 0
; vsh = vs*log(10.)*conform(vs,kh,1)
 vsh = vs

; vsh = where(vsh .eq. 0., 1.e-30, vsh)

; wl = where(kh .ne. 0., 1./kh, 999)

 sumsp = dim_sum(vs)
 sumsp1 = vs*0.
 do io=2, no*2-3
   sumsp1(:,io) = dim_sum(vs(:,io:no*2-2))/sumsp
 end do
 do iw=0, 4
   print(omea(maxind(vsh(iw,:))))
   print((1./omea(minind(abs(sumsp1(iw,:)-0.1))))+"  "+ \
         (1./omea(minind(abs(sumsp1(iw,:)-0.9)))))
 end do

;--------------------------------------- cn level
 lcol = (/18,16,13,10,1/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "frequency [cyc hr:S:-1:N:]"
ystr = "PSD [10:S:-8:N: m:S:2:N: s:S:-2:N: / cyc hr:S:-1:N:]"

 il = 0

 vsh(3,:) = vsh(3,:)*0.1

 plot = xylineplot(wks,omea,vsh*1.e8,"D",0,"T",1,"C",lcol)

 set_size(plot,xf(il),yf(il),wi,he)
 axis_range(plot,-0.5,2,-999,999,False,False)
 axis_str(plot,mstr(il),xstr(il),ystr(il))
 setvalues plot
  "tmXUseBottom" : False
  "tmXTLabelsOn" : True
  "tmXTMode"   : "Explicit"
  "tmXTValues" : 1./(/-0.5,-1,-2,-3,-6,-24,24,6,3,2,1,0.5/)
  "tmXTLabels" : (/0.5,1,2,3,6,24,24,6,3,2,1,0.5/)
  "tmXTLabelFontHeightF" : 0.004
  "tmXTMajorLengthF" : -0.0025
 end setvalues
 dr_txt(wks,"period [hr]",xf(0)+wi*0.5,yf(0)-4.2,22,4.5,0)

 draw(plot)

 do iw=0, 4
   txt = "W"+(iw+1)
   if (iw .eq. 3) then
     txt = txt+" / 10"
   end if
   dr_txt2(wks,txt,xf(0)+wi-7.5,yf(0)+1+(iw*2.5),22,4,0,0,lcol(iw))
 end do

 frame(wks)

end

