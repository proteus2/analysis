load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj","merra"/)
; expname = (/"uj","era-int"/)

 y0 = 1953
 y1 = 1979
 y2 = 2006
 lat0 = 5.
 lat0r = lat0 + 1.0

;do k0=0, 1

 k0 = 1

 z0 = (/21.0,26.9/)
 p0 = (/50,  20  /)

 z = z0(k0)
 p = p0(k0)

 wks = cr_wks("fig/fu-all_u-1m-z"+numeric2int(z,1)+"-"+expname(0)+"-"+expname(1),600,0)
 colind = NhlNewColor(wks,0.60,0.60,0.60)

 if (k0 .eq. 0) then
   udist = fspan(-16.5,10.5,10)
   ymax0 = (/10,4.5/)
   ymin0 = (/3,2/)
 end if
 if (k0 .eq. 1) then
   udist = fspan(-30.5,14.5,10)
   ymax0 = (/20,5/)
   ymin0 = (/15,4/) ;0.15*30.
 end if
; udist = fspan(-33,27,11)
; udist = fspan(-40,32,10)
 n_int = 10 ; 3-day interpol.
 nwg = 2

;--------------------------------------- read data
 ne = dimsizes(expname)
; nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2-y0)*12 - 1
 if (y2 .eq. 2006) then
   it9 = it9 - 1  ; Dec 2006 is missing
 end if

   nt = it9 - it0 + 1
   ntr= 336 - 1  ; Jan 1979 is missing
   nw = 8

dirn = "MERRA"
varn = "u"
if (expname(1) .eq. "era-int") then
  dirn = "ERA-inter"
  varn = "U_GDS0_ISBL"
end if

   fnam1 = "$DATD/L60CGW-t/epf-wc-1m/uan"+expname(0)+".epfsa_yz_recon0_1m.nc"
   fnam2 = "$DATD/L60CGW-t/epf-wc-1m/uan"+expname(0)+".epf_yz_recon0_1m.nc"
   fnam1r = "$DATD/tem/"+dirn+"/epf/"+expname(1)+".epfsa_yz_recon0.nc"
   fnam2r = "$DATD/tem/"+dirn+"/epf/"+expname(1)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   f1r = addfile(fnam1r,"r")
   f2r = addfile(fnam2r,"r")

   rw0 = new((/nt,nw/),float)
   rw0(:,:) = dim_avg( f1->epd(it0:it9,:nw-1,{z*1.e3},{-lat0:lat0}) )
   rw0(:,nw-2:) = rw0(:,nw-2:) + \
                  dim_avg( f2->epd(it0:it9,8:9,{z*1.e3},{-lat0:lat0}) )
   rw0r = new((/ntr,nw/),float)
   rw0r(:,:) = dim_avg( f1r->epd(:ntr-1,:nw-1,{p},{-lat0r:lat0r}) )
   rw0r(:,nw-2:) = rw0r(:,nw-2:) + \
                   dim_avg( f2r->epd(:ntr-1,8:9,{p},{-lat0r:lat0r}) )

; migrating tide
;   rw0(:,nw-1) = rw0(:,nw-1) + \
;                    ( dim_avg( f2->epd(it0:it9,2,:,{-lat0:lat0}) ) + \
;                      dim_avg( f2->epd(it0:it9,3,:,{-lat0:lat0}) ) )
;   rw0r(:,nw-1) = rw0r(:,nw-1) + \
;                    ( dim_avg( f2r->epd(:ntr-1,2,:,{-lat0r:lat0r}) ) + \
;                      dim_avg( f2r->epd(:ntr-1,3,:,{-lat0r:lat0r}) ) )

   rw0(:,3) = rw0(:,3) + dim_avg( f1->epd(it0:it9,13,{z*1.e3},{-lat0:lat0}) )
   rw0(:,6:7) = rw0(:,6:7) + dim_avg( f1->epd(it0:it9,16:17,{z*1.e3},{-lat0:lat0}) )
   rw0r(:,3) = rw0r(:,3) + dim_avg( f1r->epd(:ntr-1,13,{p},{-lat0r:lat0r}) )
   rw0r(:,6:7) = rw0r(:,6:7) + dim_avg( f1r->epd(:ntr-1,16:17,{p},{-lat0r:lat0r}) )

   fnamg = "$DATD/L60CGW-t/du_xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")

   bgw = dim_avg( fg->du_gw      (it0:it9,{z*1.e3},{-lat0:lat0}) )*72.
   cgw = dim_avg( fg->dudt_gwconv(it0:it9,{z*1.e3},{-lat0:lat0}) )*86400.
   gw = bgw + cgw


   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"
   fnam1r = "$DATD/clim/"+dirn+"/"+expname(1)+".u_ypt.1979-2006.nc"

   f1 = addfile(fnam1,"r")
   f1r = addfile(fnam1r,"r")

   u0 = dim_avg( f1->u(it0:it9,{z*1.e3},{-lat0:lat0}) )
   u0r = dim_avg( f1r->$varn$(:ntr-1,{p},{-lat0r:lat0r}) )

   u0b = dim_avg( f1->u(it0:it9,{21.e3},{-lat0:lat0}) )
   u0br = dim_avg( f1r->$varn$(:ntr-1,{50},{-lat0r:lat0r}) )

   fnam2 = "$DATD/L60CGW-t/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   fnam2r = "$DATD/tem/"+dirn+"/"+expname(1)+".tem_ypt.1979-2006.nc"
   f2 = addfile(fnam2,"r")
   f2r = addfile(fnam2r,"r")

   advz = dim_avg( f2->uadv_z(it0:it9,{z*1.e3},{-lat0:lat0}) )
   advzr = dim_avg( f2r->uadv_z(:ntr-1,{p},{-lat0r:lat0r}) )
;   advz = advz + dim_avg( f2->uadv_y(it0:it9,{z*1.e3},{-lat0:lat0}) )
;   advzr = advzr + dim_avg( f2r->uadv_y(:ntr-1,{p},{-lat0r:lat0r}) )

   ut = dim_avg( f2->u_tend(it0:it9,{z*1.e3},{-lat0:lat0}) )
   utr = dim_avg( f2r->u_tend(:ntr-1,{p},{-lat0r:lat0r}) )

   epd = dim_avg( f2->epd(it0:it9,{z*1.e3},{-lat0:lat0}) )
   epdr = dim_avg( f2r->epd(:ntr-1,{p},{-lat0r:lat0r}) )

   cor = dim_avg( f2->cor(it0:it9,{z*1.e3},{-lat0:lat0}) )
   corr = dim_avg( f2r->cor(:ntr-1,{p},{-lat0r:lat0r}) )

   uforce = dim_avg( f2->u_force(it0:it9,{z*1.e3},{-lat0:lat0}) )
   uforcer = dim_avg( f2r->u_force(:ntr-1,{p},{-lat0r:lat0r}) )

   uforce = uforce + gw

do l0=0, 1

if (l0 .eq. 0) then
   wv = (/ rw0(:,2), rw0(:,3), rw0(:,4), rw0(:,5), rw0(:,6), rw0(:,7), \
           ut, gw, cgw, uforce, epd, advz, u0b /)
   t = y1 + (fspan(1,nt,nt)-0.5)/12.
else
   wv = (/ rw0r(:,2), rw0r(:,3), rw0r(:,4), rw0r(:,5), rw0r(:,6), rw0r(:,7), \
           utr, utr*0., utr*0., uforcer, epdr, advzr, u0br /)
   t = 1979. + (fspan(2,ntr+1,ntr)-0.5)/12.
   nt = ntr
   delete(u0)
   u0 = u0r
end if

   nd2 = dimsizes(wv)
   nvi = nd2(0)

;--------------------------------------- process 0
 wv(:nvi-2,:) = wv(:nvi-2,:)*30.

 ndist = dimsizes(udist)
 itv = udist(1) - udist(0)
 uitv = (/udist-0.5*itv,udist+0.5*itv/)

 nti = (nt-1)*n_int+1
 ti  = fspan(t(0),t(nt-1),nti)

 ui  = linint1(t,u0,False,ti,0)
 wvi = linint1(t,wv,False,ti,0)

; determine criterion vari.
; cr  = sh
 cri = new((/2,nti/),float)

 zeropt = new(nti,integer)
 zeropt = -999

 tag = 0
 nlast = -999
 nlast2 = -999
 cnt = 0
 do n=1, nti-1
   if ( ui(n)*ui(n-1) .le. 0. ) then 
     zeropt(cnt) = n
     if (n-nlast2 .lt. 6*n_int) then
       zeropt(cnt-1:cnt) = -999
       cnt = cnt - 1
     else
       cnt = cnt + 1
     end if
     nlast2 = nlast
     nlast = n
   end if
 end do
 nzero = cnt

 maxmin = new((/2,nzero-1/),integer)
 do cnt=0, nzero-2
   if ( sum(ui(zeropt(cnt):zeropt(cnt+1))) .gt. 0. ) then
     maxmin(0,cnt) = maxind(ui(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = 1 
   else
     maxmin(0,cnt) = minind(ui(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = -1
   end if
 end do

 cri = 0.
 do cnt=0, nzero-3
   cri(0,maxmin(0,cnt):maxmin(0,cnt+1)) = tofloat(maxmin(1,cnt))*(-1.0)
   ; test
   if (maxmin(1,cnt)*maxmin(1,cnt+1) .ge. 0) then
     print("Check here !!")
   end if
 end do
 ncyc = new(2,float)
 ncyc(:) = tofloat((nzero-2)/2)
 if (mod(nzero-2,2) .ne. 0) then
   if (maxmin(1,0) .lt. 0) then
     ncyc(0) = ncyc(0) + 1.0
   else
     ncyc(1) = ncyc(1) + 1.0
   end if
 end if
print(ncyc)

 dura  = new((/2,ndist/),float)
 dura  = 0.
 wvd  = new((/2,nvi,ndist/),float)
 wvd  = 0.
 ncom  = new((/2,ndist/),float)
 ncom  = 0.
 do n=0, nti-1
   if (cri(0,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom(0,ii) = ncom(0,ii) + 1.
         dura(0,ii) = dura(0,ii) + 1.
         wvd(0,:,ii) = wvd(0,:,ii) + wvi(:,n)
       end if
     end do
   end if
   if (cri(0,n) .lt. 0.0) then
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom(1,ii) = ncom(1,ii) + 1.
         dura(1,ii) = dura(1,ii) + 1.
         wvd(1,:,ii) = wvd(1,:,ii) + wvi(:,n)
       end if
     end do
   end if
 end do
 dura@_FillValue = 0.
 ncom@_FillValue = 0.

 do is=0, 1
   dura(is,:) = dura(is,:)*(tofloat(nt)/tofloat(nti))/ncyc(is)
 end do

 do iv=0, nvi-1
   wvd(:,iv,:) = wvd(:,iv,:)/ncom(:,:)
 end do

print(ncom)

 nv1 = 8
 vard = new((/2,nv1,ndist/),float)
 vard(:,0,:) = wvd(:,6,:)  ; UT
 vard(:,1,:) = wvd(:,11,:) ; ADVZ
 vard(:,2,:) = wvd(:,6,:) - (wvd(:,9,:) - wvd(:,10,:) - wvd(:,7,:)) ; REQ
 vard(:,3,:) = wvd(:,0,:)  ; K
 vard(:,4,:) = wvd(:,7,:)  ; GW
 vard(:,5,:) = wvd(:,8,:)  ; CGW
 vard(:,6,:) = wvd(:,6,:) - wvd(:,9,:)  ; RESIDUAL
 vard(:,7,:) = vard(:,2,:) - wvd(:,10,:) ; REQ - EPD
; vard(:,7,:) = wvd(:,12,:) ; U50

 nv2 = 5
 vard2 = new((/2,nv2,ndist/),float)
 vard2(:,0,:) = wvd(:,10,:)  ; EPD
 vard2(:,1,:) = wvd(:,0,:)  ; K
 vard2(:,2,:) = wvd(:,1,:)  ; R
 vard2(:,3,:) = wvd(:,2,:) + wvd(:,3,:)  ; MRG
 vard2(:,4,:) = wvd(:,4,:) + wvd(:,5,:)  ; IG

; vard2(:,5,:) = wvd(:,3,:)  ; test ; will be shown in sky-blue
; nv2 = nv2 + 1

;--------------------------------------- cn level
 lcol = (/1,8,1,1,2,2,6,colind/)
 ldsh = (/16,0,0,0,0,16,0,0/)
 ltck = (/2.5,2.5,2.5,1,1,1,1,2.5/)*1.0
 ltck(0) = 2.5

 lcol2 = (/colind,1,4,7,3,6,8/)
 ldsh2 = (/16,0,0,0,0,0,0/)
 ltck2 = (/3,1,1,1,1,1,1/)*0.8
; ltck2(0) = ltck(1)

 if (l0 .ne. 0) then
   lcol(4:6) = -1
 end if

;--------------------------------------- contour
; l0 = k0

 nfig = 2
 xmin = (/min(uitv(0,:)),max(uitv(1,:))/)
 xmax = (/max(uitv(1,:)),min(uitv(0,:))/)

 wi = 32
 he = (/30,20/)
 xf = 10 + (/0,1/)*40 ;+ l0*88
 yf = 10 + (/0,1,2/)*36 + l0*65

 mstr = new((/nwg,nfig/),string,"")
 xstr = new((/nwg,nfig/),string,"")
 ystr = new(nfig,string,"")
 mstr(0,:) = (/"Easterly to westerly","Westerly to easterly"/)
 if (l0 .ne. 0) then
   mstr(:,:) = ""
 end if
 xstr(nwg-1,:) = (/"U:B:EW:N: (m s:S:-1:N:)","U:B:WE:N: (m s:S:-1:N:)"/)
 ystr(0) = "force (m s:S:-1:N: month:S:-1:N:)"
 txtab = (/"(a)","(b)"/)

do iwg=0, nwg-1
do ifig=0, nfig-1

 ymin = (/-ymin0(iwg),ymin0(iwg)/)
 ymax = (/ymax0(iwg),-ymax0(iwg)/)

 plot = xylineplot(wks,udist*2.,0.,"D",0,"T",0.6,"C",1)
 if (iwg .eq. 0) then
   xyaddline(plot,udist,vard(ifig,:nv1-1,:),"D",ldsh,"T",ltck,"C",lcol)
 else
   xyaddline(plot,udist,vard2(ifig,:,:),"D",ldsh2,"T",ltck2,"C",lcol2)
 end if

 set_size(plot,xf(ifig),yf(iwg),wi,he(iwg))
 axis_range(plot,xmin(ifig),xmax(ifig),ymin(ifig),ymax(ifig),False,False)
 axis_str(plot,mstr(iwg,ifig),xstr(iwg,ifig),ystr(ifig))
 setvalues plot
  "tmXBMode" : "Explicit"
  "tmXBValues" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
  "tmXBLabels" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

end do
end do

if (l0 .eq. 0) then
 lg = lgplot(wks,(/":F34:6:F21:U/:F34:6:F21:t","ADVz","REQ","Kelvin","CGW :F18:+:F21: BGW","CGW", \
                   "D:B:u:H-10V15:_:N:","REQ :F18:-:F21: Resolved"/), \
        "D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(1)+wi+2.0,yf(0)+0.4,wi*0.4,nv1*3)
 setvalues lg
  "lgLineDashSegLenF" : 0.025
 end setvalues
 draw(lg)

 lg2 = lgplot(wks,(/"Resolved","Kelvin","Rossby","MRG","IG"/), \
        "D",ldsh2,"T",ltck2,"C",lcol2,2)
 set_size_lg(lg2,xf(1)+wi+2.0,yf(1)+0.4,wi*0.4,nv2*3)
 setvalues lg2
  "lgLineDashSegLenF" : 0.025
 end setvalues
 draw(lg2)
end if

 delete(maxmin)
 delete(ncom)

 dr_txt(wks,txtab(l0),xf(0)-6,yf(0)-2,25,8,0.)
; dr_txt(wks,numeric2int(z,1)+" km",xf(0)+5,yf(0)+2,22,6,0.)

;end do  ; k0

 delete([/wv,t,ti,ui,wvi,cri,zeropt/])

end do  ; l0

 frame(wks)

end

