load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"tc"/)
 ppid    = (/"j" /)
 varname = "u"

 y1 = 1951
 y2 = 1970
 lat0 = 5.

 k0 = 3

 z0 = (/21.0,23.3,26.9,30.6/)
 p0 = (/50,  30,  20  ,10  /)

 nwg = 2

 z = z0(k0)
 p = p0(k0)

 if (k0 .eq. 0) then
   udist = fspan(-16.5,10.5,10)
   ymax0 = (/5,3/)
   ymin0 = (/3,2/)
 end if
 if (k0 .eq. 1) then
   udist = fspan(-22,14,10)
   ymax0 = (/12,4/)
   ymin0 = (/9,2/) ;0.15*30.
 end if
 if (k0 .eq. 2) then
   udist = fspan(-30.5,19.5,11)
   ymax0 = (/25,5/)
   ymin0 = (/15,3/) ;0.15*30.
 end if
 if (k0 .eq. 3) then
   udist = fspan(-30.5,19.5,11)
   ymax0 = (/35,5/)
   ymin0 = (/25,3/) ;0.15*30.
 end if

; udist = fspan(-33,27,11)
; udist = fspan(-40,32,10)
 n_int = 10 ; 3-day interpol.

 wks = cr_wks("fig/fu-rw_avg"+floattoint(z)+"-"+expname(0),600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new(ne,string)
 do ie=0, ne-1
   fname = "$DATD/AOL60CGW/xm.uan"+expname(ie)+"a.pj.nc"
print(expname(ie)+"")
   f1 = addfile(fname(ie),"r")
   if (ie .eq. 0) then
;     nt = dimsizes(f1->t)
     u  = new((/ne,nt/),float,1.e20)
     ut = new((/ne,nt/),float,1.e20)
     sh = new((/ne,nt/),float,1.e20)
     ht = f1->hybrid_ht
     kk = closest_val(z*1.e3,ht)
     if (ht(kk) .ne. ht({z*1.e3})) then
       print("wrong index")
     end if
   end if
   u (ie,:) = dim_avg( f1->$varname$(:nt-1,{z*1.e3},{-lat0:lat0}) )
   sh(ie,:) = ( dim_avg( f1->$varname$(:nt-1,kk+1,{-lat0:lat0}) ) -  \
                dim_avg( f1->$varname$(:nt-1,kk-1,{-lat0:lat0}) ) )/ \
              ( (ht(kk+1)-ht(kk-1))/1.e3 )  ; m/s/km
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("$DATD/AOL60CGW/tem/uan"+expname(0)+".tem_yzt.nc","r")
 utend  = dim_avg( f2->u_tend (:,{z*1.e3},{-lat0:lat0}) )
 uforce = dim_avg( f2->u_force(:,{z*1.e3},{-lat0:lat0}) )
 epd    = dim_avg( f2->epd    (:,{z*1.e3},{-lat0:lat0}) )
 cor    = dim_avg( f2->cor    (:,{z*1.e3},{-lat0:lat0}) )
 advy   = dim_avg( f2->uadv_y (:,{z*1.e3},{-lat0:lat0}) )
 advz   = dim_avg( f2->uadv_z (:,{z*1.e3},{-lat0:lat0}) )

 f2a = addfile("$DATD/AOL60CGW/epf/uan"+expname(0)+".epfsa0_yz_recon0.nc","r")
 f2b = addfile("$DATD/AOL60CGW/epf/uan"+expname(0)+".epf_yz_recon0.nc","r")
 kw = dim_avg( f2a->epd(:,2,{z*1.e3},{-lat0:lat0}) )
 rw = dim_avg( f2a->epd(:,3,{z*1.e3},{-lat0:lat0}) )
 ergw = dim_avg( f2a->epd(:,4,{z*1.e3},{-lat0:lat0}) )
 wrgw = dim_avg( f2a->epd(:,5,{z*1.e3},{-lat0:lat0}) )
 eigw = dim_avg( f2b->epd(:,8,{z*1.e3},{-lat0:lat0}) )
 wigw = dim_avg( f2b->epd(:,9,{z*1.e3},{-lat0:lat0}) )
print(dimsizes(ut))
print(dimsizes(utend))

 ut(0,:) = utend(:)

 f3 = addfile("$DATD/AOL60CGW/xm.uan"+expname(0)+"a.pe.nc","r")
 du_gwb = dim_avg( f3->du_gw(:,{z*1.e3},{-lat0:lat0}) )*72.
 du_gw = du_gwb
 if (cexp(0) .ne. "CTL") then
   du_gwc = dim_avg( f3->dudt_gwconv(:,{z*1.e3},{-lat0:lat0}) )*86400.
   du_gw = du_gw + du_gwc
 end if

 uforce = uforce + du_gw

print(max(u)+"  "+min(u))

;--------------------------------------- process 0
 ndist = dimsizes(udist)
 itv = udist(1) - udist(0)
 uitv = (/udist-0.5*itv,udist+0.5*itv/)

 nti = (nt-1)*n_int+1
 ti  = fspan(t(0),t(nt-1),nti)

 ui  = linint1(t,u (0,:),False,ti,0)
 uti = linint1(t,ut(0,:),False,ti,0)
 shi = linint1(t,sh(0,:),False,ti,0)
 epdi = linint1(t,epd,False,ti,0)
 kwi = linint1(t(:nt-2),kw,False,ti,0)
 rwi = linint1(t(:nt-2),rw,False,ti,0)
 ergwi = linint1(t(:nt-2),ergw,False,ti,0)
 wrgwi = linint1(t(:nt-2),wrgw,False,ti,0)
 eigwi = linint1(t(:nt-2),eigw,False,ti,0)
 wigwi = linint1(t(:nt-2),wigw,False,ti,0)
 cori = linint1(t,cor,False,ti,0)
 advyi= linint1(t,advy,False,ti,0)
 advzi= linint1(t,advz,False,ti,0)
 gwdi = linint1(t,du_gw,False,ti,0)
 if (cexp(0) .ne. "CTL") then
   gwdci = linint1(t,du_gwc,False,ti,0)
 end if
 uforcei = linint1(t,uforce,False,ti,0)

 vari = (/uti,kwi,gwdi,gwdi,advzi,uti-uforcei/)
 vari2 = (/kwi,rwi,ergwi+wrgwi,eigwi+wigwi/)
 nvi = 6
 nvi2 = 4
 if (cexp(0) .ne. "CTL") then
   vari(3,:) = gwdci(:)
 end if

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
 vard  = new((/2,nvi,ndist/),float)
 vard  = 0.
 vard2 = new((/2,nvi2,ndist/),float)
 vard2 = 0.
 shear = new((/2,ndist/),float,1.e20)
 shear = 0.
 ncom  = new((/2,ndist/),float)
 ncom  = 0.
 do n=0, nti-1
   if (cri(0,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom(0,ii) = ncom(0,ii) + 1.
         dura(0,ii) = dura(0,ii) + 1.
         shear(0,ii) = shear(0,ii) + shi(n)
         vard(0,:,ii) = vard(0,:,ii) + vari(:,n)
         vard2(0,:,ii) = vard2(0,:,ii) + vari2(:,n)
       end if
     end do
   end if
   if (cri(0,n) .lt. 0.0) then
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom(1,ii) = ncom(1,ii) + 1.
         dura(1,ii) = dura(1,ii) + 1.
         shear(1,ii) = shear(1,ii) + shi(n)
         vard(1,:,ii) = vard(1,:,ii) + vari(:,n)
         vard2(1,:,ii) = vard2(1,:,ii) + vari2(:,n)
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
   vard(:,iv,:) = vard(:,iv,:)/ncom(:,:)
 end do
 do iv=0, nvi2-1
   vard2(:,iv,:) = vard2(:,iv,:)/ncom(:,:)
 end do

 shear(:,:) = shear(:,:)/ncom(:,:)
print(ncom)

 vard = vard*30.
 vard2 = vard2*30.

;--------------------------------------- cn level
 lcol = (/1,4,2,2,7,3,6,8/)
 ldsh = (/0,0,0,16,0,0,0,0/)
 ltck = (/1,1,1,1,1,1,1,1/)*1.0
 ltck(0) = 1.8

 lcol2 = (/4,2,3,7,6,8/)
 ldsh2 = (/0,0,0,0,0,0/)
 ltck2 = (/1,1,1,1,1,1/)*0.6
 ltck2(0) = ltck(1)

;--------------------------------------- contour
 nfig = 2
 xmin = (/min(uitv(0,:)),max(uitv(1,:))/)
 xmax = (/max(uitv(1,:)),min(uitv(0,:))/)

 wi = 26
 he = (/30,20/)
 xf = 8 + (/0,1/)*34
 yf = 10 + (/0,1,2/)*36

 mstr = new((/nwg,nfig/),string,"")
 xstr = new(nwg,string,"")
 ystr = new(nfig,string,"")
 mstr(0,:) = (/"Easterly to westerly","Westerly to easterly"/)
 xstr(nwg-1) = "U (m s:S:-1:N:)"
 ystr(0) = "force (m s:S:-1:N: month:S:-1:N:)"

do iwg=0, nwg-1
do ifig=0, nfig-1

 ymin = (/-ymin0(iwg),ymin0(iwg)/)
 ymax = (/ymax0(iwg),-ymax0(iwg)/)

 plot = xylineplot(wks,udist*2.,0.,"D",0,"T",0.6,"C",1)
 if (iwg .eq. 0) then
   xyaddline(plot,udist,vard(ifig,:,:),"D",ldsh,"T",ltck,"C",lcol)
 else
   xyaddline(plot,udist,vard2(ifig,:,:),"D",ldsh2,"T",ltck2,"C",lcol2)
 end if

 set_size(plot,xf(ifig),yf(iwg),wi,he(iwg))
 axis_range(plot,xmin(ifig),xmax(ifig),ymin(ifig),ymax(ifig),False,False)
 axis_str(plot,mstr(iwg,ifig),xstr(iwg),ystr(ifig))
 setvalues plot
  "tmXBMode" : "Explicit"
  "tmXBValues" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
  "tmXBLabels" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

end do
end do

 lg = lgplot(wks,(/"dU/dt","KW","GW","CGW","ADVz","DIFF"/), \
        "D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(1)+wi+2.0,yf(0)+0.4,wi*0.4,nvi*3)
 setvalues lg
  "lgLineDashSegLenF" : 0.025
 end setvalues
 draw(lg)

 lg2 = lgplot(wks,(/"KW","RW","MRGW","I-GW"/), \
        "D",ldsh2,"T",ltck2,"C",lcol2,2)
 set_size_lg(lg2,xf(1)+wi+2.0,yf(1)+0.4,wi*0.4,nvi2*3)
 setvalues lg2
  "lgLineDashSegLenF" : 0.025
 end setvalues
 draw(lg2)


 frame(wks)

end

