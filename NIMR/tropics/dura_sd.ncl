load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.5

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"

 y1 = 1953
 y2 = 2006
 lat0 = 5.

 k0 = 0

 z0 = (/21.0,23.3,26.9/)
 p0 = (/50,  30,  20  /)

 z = z0(k0)
 p = p0(k0)

 if (k0 .eq. 0) then
   udist = fspan(-16.5,10.5,10)
   ymax0 = 4
 end if
 if (k0 .eq. 1) then
   udist = fspan(-26,14,11)
   ymax0 = 4
 end if
 if (k0 .eq. 2) then
   udist = fspan(-30.5,14.5,10)
   ymax0 = 4
 end if
; udist = fspan(-33,27,11)
; udist = fspan(-40,32,10)
 n_int = 10 ; 3-day interpol.

 wks = cr_wks("fig/dura_sd"+floattoint(z)+"-"+expname(0),600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new(ne,string)
 do ie=0, ne-1
   fname = "$DATD/NIMR/u_xm.uan"+expname(ie)+"a.pj.nc"
print(expname(ie)+"")
   f1 = addfile(fname(ie),"r")
   if (ie .eq. 0) then
     nt = dimsizes(f1->t)
     u  = new((/ne,nt/),float,1.e20)
     ut = new((/ne,nt/),float,1.e20)
     sh = new((/ne,nt/),float,1.e20)
     ht = f1->hybrid_ht
     kk = closest_val(z*1.e3,ht)
     if (ht(kk) .ne. ht({z*1.e3})) then
       print("wrong index")
     end if
   end if
   u (ie,:) = dim_avg( f1->$varname$(:,{z*1.e3},{-lat0:lat0}) )
   sh(ie,:) = ( dim_avg( f1->$varname$(:,kk+1,{-lat0:lat0}) ) -  \
                dim_avg( f1->$varname$(:,kk-1,{-lat0:lat0}) ) )/ \
              ( (ht(kk+1)-ht(kk-1))/1.e3 )  ; m/s/km
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("../dat/tem/uan"+expname(0)+".tem_yzt.nc","r")
 utend  = dim_avg( f2->u_tend (:,{z*1.e3},{-lat0:lat0}) )

 ut(0,:) = utend(:)

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

 nvi = 1
 vari = new((/nvi,nti/),float)
 vari = (/uti/)

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
 shear = new((/2,ndist/),float,1.e20)
 shear = 0.
 ncom  = new((/2,ndist/),float)
 ncom  = 0.
 vsign = new(nti,integer)
 vsign = 0
 dura1  = new((/sum(floattoint(ncyc)),2,ndist/),float)
 dura1  = 0.
 vard1  = new((/sum(floattoint(ncyc)),2,nvi,ndist/),float)
 vard1  = 0.
 shear1 = new((/sum(floattoint(ncyc)),2,ndist/),float,1.e20)
 shear1 = 0.
 ncom1  = new((/sum(floattoint(ncyc)),2,ndist/),float)
 ncom1  = 0.

 j = -1
 tag = 0
 do n=0, nti-1
   if (cri(0,n) .gt. 0.0) then
     if (tag .ne. 1) then
       j = j + 1
     end if
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom1(j,0,ii) = ncom1(j,0,ii) + 1.
         dura1(j,0,ii) = dura1(j,0,ii) + 1.
         shear1(j,0,ii) = shear1(j,0,ii) + uti(n)
         vard1(j,0,:,ii) = vard1(j,0,:,ii) + vari(:,n)
         vsign(n) = ii
         ncom(0,ii) = ncom(0,ii) + ncom1(j,0,ii)
         dura(0,ii) = dura(0,ii) + dura1(j,0,ii)
         shear(0,ii) = shear(0,ii) + shear1(j,0,ii)
         vard(0,:,ii) = vard(0,:,ii) + vard1(j,0,:,ii)
       end if
     end do
     tag = 1
   end if
   if (cri(0,n) .lt. 0.0) then
     if (tag .ne. -1) then
       j = j + 1
     end if
     do ii=0, ndist-1
       if ( ui(n) .ge. uitv(0,ii) .and. ui(n) .lt. uitv(1,ii) ) then
         ncom1(j,1,ii) = ncom1(j,1,ii) + 1.
         dura1(j,1,ii) = dura1(j,1,ii) + 1.
         shear1(j,1,ii) = shear1(j,1,ii) + uti(n)
         vard1(j,1,:,ii) = vard1(j,1,:,ii) + vari(:,n)
         vsign(n) = ii
         ncom(1,ii) = ncom(1,ii) + ncom1(j,1,ii) 
         dura(1,ii) = dura(1,ii) + dura1(j,1,ii) 
         shear(1,ii) = shear(1,ii) + shear1(j,1,ii) 
         vard(1,:,ii) = vard(1,:,ii) + vard1(j,1,:,ii) 
       end if
     end do
     tag = -1
   end if
 end do
 dura@_FillValue = 0.
 ncom@_FillValue = 0.
 dura1@_FillValue = 0.
 ncom1@_FillValue = 0.

 do is=0, 1
   dura(is,:) = dura(is,:)*(tofloat(nt)/tofloat(nti))/ncyc(is)
 end do
 dura1(:,:,:) = dura1(:,:,:)*(tofloat(nt)/tofloat(nti))

 do iv=0, nvi-1
   vard(:,iv,:) = vard(:,iv,:)/ncom(:,:)
   vard1(:,:,iv,:) = vard1(:,:,iv,:)/ncom1(:,:,:)
 end do

 shear(:,:) = shear(:,:)/ncom(:,:)
 shear1(:,:,:) = shear1(:,:,:)/ncom1(:,:,:)
print(ncom)

 varsd = new(dimsizes(vard),float)
 varsd = 0.
 if ( cri(0,maxmin(0,0)) .gt. 0. ) then
   tag1 = 0
   tag2 = 1
 else
   tag1 = 1
   tag2 = 0
 end if
 varsd(0,0,:) = dim_stddev_n(dura1(ispan(tag1,floattoint(ncyc(0))-1,2),0,:),0)
 varsd(1,0,:) = dim_stddev_n(dura1(ispan(tag2,floattoint(ncyc(1))-1,2),1,:),0)

 varm = new(dimsizes(varsd),float)
 varm(0,0,:) = dim_avg_n(dura1(ispan(tag1,floattoint(ncyc(0))-1,2),0,:),0)
 varm(1,0,:) = dim_avg_n(dura1(ispan(tag2,floattoint(ncyc(1))-1,2),1,:),0)
print(sum(varm(0,0,:)))
print(sum(varm(1,0,:)))

 varo = (/varm,varsd/)

;--------------------------------------- cn level
 lcol = (/1,1,4,2,2 ,3,6,7,8,9/)
 ldsh = (/0,1,0,0,14,0,0,0,0,0/)
 ltck = (/2.0,1.5,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2/)

;--------------------------------------- contour
 nfig = 2
 xmin = (/min(uitv(0,:)),max(uitv(1,:))/)
 xmax = (/max(uitv(1,:)),min(uitv(0,:))/)
 ymin = (/0,0/)
 ymax = (/ymax0,ymax0/)

 wi = 50
 he = 60
 xf = 15 + (/0,1/)*72
 yf = 10 + (/0,1,2/)*28

 mstr  = (/"Easterly to westerly","Westerly to easterly"/)
 xstr  = "U (m s:S:-1:N:)"
 ystr  = "duration (month)"+(/"",""/)

do ifig=0, nfig-1

 plot = xylineplot(wks,udist*2.,0.,"D",0,"T",0.6,"C",1)
 do k=0, 1
   xyaddline(plot,udist,varo(k,ifig,0,:),"D",ldsh(k),"T",ltck(k),"C",lcol(k))
 end do

 set_size(plot,xf(ifig),yf,wi,he)
 axis_range(plot,xmin(ifig),xmax(ifig),ymin(ifig),ymax(ifig),False,False)
 axis_str(plot,mstr(ifig),xstr,ystr(ifig))
 setvalues plot
  "tmXBMode" : "Explicit"
  "tmXBValues" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
  "tmXBLabels" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

end do

 lg = lgplot(wks,(/"mean","Std. Dev."/), \
        "D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(1)+wi+5.0,yf+0.7,wi*0.40,15)
 setvalues lg
  "lgLabelOffsetF" : 0.2
  "lgLineDashSegLenF"  : 0.030
  "lgLabelFontHeightF" : 0.010
 end setvalues
 draw(lg)


 frame(wks)

end



