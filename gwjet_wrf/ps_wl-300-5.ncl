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
 it1 = ispan(144,431,1)
 it2 = ispan(0,287,1)
 
 wks = cr_wks("fig/ps_wl-300"+var1name+"-5",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_wn__z"+sprinti("%5.5i",lev)+"-w"+ispan(1,5,1)+".nc"

 f1 = addfile(fname1(0),"r")
 f2 = addfile(fname1(1),"r")
 f3 = addfile(fname1(2),"r")
 f4 = addfile(fname1(3),"r")
 f5 = addfile(fname1(4),"r")

 kh = f1->kh / 111.
 nkh = dimsizes(kh)

 vs1 = dim_avg_n( dim_sum( f1->psd_phi_kh(it1,:,:) ), 0 )
 vs2 = dim_avg_n( dim_sum( f2->psd_phi_kh(it1,:,:) ), 0 )
 vs3 = dim_avg_n( dim_sum( f3->psd_phi_kh(it1,:,:) ), 0 )
 vs4 = dim_avg_n( dim_sum( f4->psd_phi_kh(it2,:,:) ), 0 )
 vs5 = dim_avg_n( dim_sum( f5->psd_phi_kh(it2,:,:) ), 0 )
 vs = (/vs1,vs2,vs3,vs4,vs5/) * 111. * dphi

 fname0 = fdir+"/x6_5min/"+var1+"_z/prt_d300km/psd/psd_"+var1name+"_wn__z"+sprinti("%5.5i",lev)+"-"+ispan(1,2,1)+"a.nc"

 f01 = addfile(fname0(0),"r")
 f02 = addfile(fname0(1),"r")

 vs01 = dim_avg_n( dim_sum( f01->psd_phi_kh(it1,:,:) ), 0 )
 vs02 = dim_avg_n( dim_sum( f02->psd_phi_kh(it2,:,:) ), 0 )

 print( sum(vs1)/sum(vs01)+" + "+sum(vs2)/sum(vs01)+" + "+sum(vs3)/sum(vs01) )
 print( " = "+((sum(vs1)+sum(vs2)+sum(vs3))/sum(vs01)) )
 print( sum(vs4)/sum(vs02)+" + "+sum(vs5)/sum(vs02))
 print( " = "+((sum(vs4)+sum(vs5))/sum(vs02)) )

;--------------------------------------- process 0
 vsh = vs*log(10.)*conform(vs,kh,1)

 vsh = where(vsh .eq. 0., 1.e-30, vsh)

 wl = where(kh .ne. 0., 1./kh, 999)

 sumsp = dim_sum(vs(:,1:))
 sumsp1 = vs*0.
 do ikh=2, nkh-2
   sumsp1(:,ikh) = dim_sum(vs(:,ikh:nkh-1))/sumsp
 end do
 do iw=0, 4
   print(wl(maxind(vsh(iw,:))))
   print(wl(minind(abs(sumsp1(iw,:)-0.1)))+"  "+wl(minind(abs(sumsp1(iw,:)-0.9))))
 end do

;--------------------------------------- cn level
 lcol = (/18,16,13,10,1/)

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

mstr = ""
xstr = "horizontal wavelength [km]"
ystr = "PSD [10:S:-5:N: m:S:2:N: s:S:-2:N:]"

 il = 0

 vsh(3,:) = vsh(3,:)*0.5

 plot = xylineplot(wks,wl,vsh*1.e5,"D",0,"T",1,"C",lcol)

 set_size(plot,xf(il),yf(il),wi,he)
 axis_range(plot,999,10,-999,999,True,False)
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



