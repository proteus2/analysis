load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+144
 lat0 = 5.
 kz  = (/55,57,59/)
 kz0 = (/7,8,9/)
 ncyc = (/ (/ (/5.5,5.8,6./),(/5.5,5.2,5. /) /), \
           (/ (/5. ,5. ,5./),(/5. ,5. ,5. /) /), \
           (/ (/5. ,5. ,5./),(/4.8,4.8,4.8/) /) /)

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*12

 udist = fspan(-36,24,11)
; udist = fspan(-33,27,11)
; udist = fspan(-40,32,10)
 n_int = 15 ; 2-day interpol.

 wks = cr_wks("fig/hist2_u-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/u_"+expname(0)+".nc"
   fnam1a = "../dat/ng_"+expname(0)+".nc"
   fnam1b = "../dat/cg_"+expname(0)+".nc"

   f1 = addfile(fnam1 ,"r")
   f1a= addfile(fnam1a,"r")
   f1b= addfile(fnam1b,"r")
   u1 = dim_avg( f1->u(4:nmon-1,kz,{-lat0:lat0},0) )
   z1 = f1->hybrid_ht(kz)/1.e3
   dz1 = f1->hybrid_ht(kz+1)/1.e3 - f1->hybrid_ht(kz-1)/1.e3
   sh1 = dim_avg( f1->u(4:nmon-1,kz+1,{-lat0:lat0},0) - \
                  f1->u(4:nmon-1,kz-1,{-lat0:lat0},0) )
   ng1 = dim_avg( f1a->field424   (4:nmon-1,kz,{-lat0:lat0},0) )*86400.
   cg1 = dim_avg( f1b->unspecified(4:nmon-1,kz,{-lat0:lat0},0) )*86400.

   p1 = 1.e3*exp(-z1/7.)
print(z1+"  "+p1)

   fnam2 = "../dat/u_"+expname(1)+".nc"
   fnam2a = "../dat/ng_"+expname(1)+".nc"

   f2 = addfile(fnam2 ,"r")
   f2a= addfile(fnam2a,"r")
   u2 = dim_avg( f2->u(4:nmon-1,kz,{-lat0:lat0},0) )
   sh2 = dim_avg( f2->u(4:nmon-1,kz+1,{-lat0:lat0},0) - \
                  f2->u(4:nmon-1,kz-1,{-lat0:lat0},0) )
   ng2 = dim_avg( f2a->field424(4:nmon-1,kz,{-lat0:lat0},0) )*86400.

   nz = dimsizes(p1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,kz0,{-lat0:lat0}) )
   sh3 = dim_avg( f3->u(m1:m1+nmon0-1,kz0+1,{-lat0:lat0}) - \
                  f3->u(m1:m1+nmon0-1,kz0-1,{-lat0:lat0}) )
   p  = f3->levels(kz0)
   z3 = 7.0*log(1.e3/p)
   dz3 = 7.0*log(f3->levels(kz0-1)/f3->levels(kz0+1))
print(p)

   nt3 = nmon0

   u  = new((/2,nz,nt/),float)
   ut = new((/2,nz,nt/),float)
   sh = new((/2,nz,nt/),float)
   gw = new((/4,nz,nt/),float)
   do n=0, nt-1
     u (:,:,n) = (/u2 (n,:),u1 (n,:)/)
     sh(:,:,n) = (/ sh2(n,:)/dz1(:), sh1(n,:)/dz1(:) /)
     gw(:,:,n) = (/ng2(n,:),ng1(n,:)+cg1(n,:),ng1(n,:),cg1(n,:)/)
   end do
   do n=1, nt-2
     ut(:,:,n) = (u(:,:,n+1)-u(:,:,n-1))/60.
   end do
   ut(:,:,0   ) = (u(:,:,1   )-u(:,:,0   ))/30.
   ut(:,:,nt-1) = (u(:,:,nt-1)-u(:,:,nt-2))/30.

   ur  = new((/nz,nt3/),float)
   utr = new((/nz,nt3/),float)
   shr = new((/nz,nt3/),float)
   do n=0, nt3-1
     ur (:,n) = u3(n,:)
     shr(:,n) = sh3(n,:)/dz3(:)
   end do
   do n=1, nt3-2
     utr(:,n) = (ur(:,n+1)-ur(:,n-1))/60.
   end do
   utr(:,0    ) = (ur(:,1    )-ur(:,0    ))/30.
   utr(:,nt3-1) = (ur(:,nt3-1)-ur(:,nt3-2))/30.

   t = fspan(1,nt,nt)
   t3 = fspan(1,nt3,nt3)

;   uss = new((/3,nz,12/),float)
;   do n=0, 11
;     uss(0:1,:,n) = dim_avg(u(0:1,:,ispan(n,nt-1,12)))
;     uss(2  ,:,n) = dim_avg(ur(:,ispan(n,nt3-1,12)))
;   end do
;   umean = dim_avg(uss)
;
;   udes  = new(dimsizes(u),float)
;   udes0 = new(dimsizes(ur),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(0:1,:,mod(n,12))
;   end do
;   do n=0, nt3-1
;     udes0(:,n) = ur(:,n) - uss(2,:,mod(n,12))
;   end do
;
;   do k=0, nz-1
;     udes(0,k,:) = udes(0,k,:) + umean(0,k)
;     udes(1,k,:) = udes(1,k,:) + umean(1,k)
;     udes0 (k,:) = udes0 (k,:) + umean(2,k)
;   end do

;--------------------------------------- process 0
 ndist = dimsizes(udist)
 itv = udist(1) - udist(0)
 uitv = (/udist-0.5*itv,udist+0.5*itv/)

 nti = (nt-1)*n_int+1
 ti  = fspan(1,nt,nti)

 ui  = linint1(t,u ,False,ti,0)
 uti = linint1(t,ut,False,ti,0)
 shi = linint1(t,sh,False,ti,0)
 gwi = linint1(t,gw,False,ti,0)

 nt3i = (nt3-1)*n_int+1
 t3i  = fspan(1,nt3,nt3i)

 u3i  = linint1(t3,ur ,False,t3i,0)
 ut3i = linint1(t3,utr,False,t3i,0)
 sh3i = linint1(t3,shr,False,t3i,0)

 dura = new((/2,3,nz,ndist/),float)
 dura = 0.
 do i=0, 1
 do k=0, nz-1
 do n=0, nt-1
   if (ut(i,k,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( u(i,k,n) .ge. uitv(0,ii) .and. u(i,k,n) .lt. uitv(1,ii) ) then
         dura(0,i,k,ii) = dura(0,i,k,ii) + 1.
       end if
     end do
   else
     do ii=0, ndist-1
       if ( u(i,k,n) .ge. uitv(0,ii) .and. u(i,k,n) .lt. uitv(1,ii) ) then
         dura(1,i,k,ii) = dura(1,i,k,ii) + 1.
       end if
     end do
   end if
 end do
 end do
 end do
 do k=0, nz-1
 do n=0, nt3-1
   if (utr(k,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ur(k,n) .ge. uitv(0,ii) .and. ur(k,n) .lt. uitv(1,ii) ) then
         dura(0,2,k,ii) = dura(0,2,k,ii) + 1.
       end if
     end do
   else
     do ii=0, ndist-1
       if ( ur(k,n) .ge. uitv(0,ii) .and. ur(k,n) .lt. uitv(1,ii) ) then
         dura(1,2,k,ii) = dura(1,2,k,ii) + 1.
       end if
     end do
   end if
 end do
 end do
 do is=0, 1
 do i=0, 2
 do k=0, nz-1
   dura(is,i,k,:) = dura(is,i,k,:)/ncyc(i,is,k)
 end do
 end do
 end do
 dura@_FillValue = 0.

 tend  = new((/2,3,nz,ndist/),float,1.e20)
 tend  = 0.
 gwd   = new((/2,4,nz,ndist/),float,1.e20)
 gwd   = 0.
 shear = new((/2,3,nz,ndist/),float,1.e20)
 shear = 0.
 ncom  = new((/2,5,nz,ndist/),float)
 ncom  = 0.
 do i=0, 3
 do k=0, nz-1
 do n=0, nti-1
   i0 = min((/1,i/))
   if (uti(i0,k,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( ui(i0,k,n) .ge. uitv(0,ii) .and. ui(i0,k,n) .lt. uitv(1,ii) ) then
         ncom(0,i,k,ii) = ncom(0,i,k,ii) + 1.
         gwd (0,i,k,ii) = gwd (0,i,k,ii) + gwi(i,k,n)
         if (i .lt. 2) then
           tend (0,i,k,ii) = tend (0,i,k,ii) + uti(i,k,n)
           shear(0,i,k,ii) = shear(0,i,k,ii) + shi(i,k,n)
         end if
       end if
     end do
   else
     do ii=0, ndist-1
       if ( ui(i0,k,n) .ge. uitv(0,ii) .and. ui(i0,k,n) .lt. uitv(1,ii) ) then
         ncom(1,i,k,ii) = ncom(1,i,k,ii) + 1.
         gwd (1,i,k,ii) = gwd (1,i,k,ii) + gwi(i,k,n)
         if (i .lt. 2) then
           tend (1,i,k,ii) = tend (1,i,k,ii) + uti(i,k,n)
           shear(1,i,k,ii) = shear(1,i,k,ii) + shi(i,k,n)
         end if
       end if
     end do
   end if
 end do
 end do
 end do

 do k=0, nz-1
 do n=0, nt3i-1
   if (ut3i(k,n) .gt. 0.0) then
     do ii=0, ndist-1
       if ( u3i(k,n) .ge. uitv(0,ii) .and. u3i(k,n) .lt. uitv(1,ii) ) then
         ncom (0,4,k,ii) = ncom (0,4,k,ii) + 1.
         tend (0,2,k,ii) = tend (0,2,k,ii) + ut3i(k,n)
         shear(0,2,k,ii) = shear(0,2,k,ii) + sh3i(k,n)
       end if
     end do
   else
     do ii=0, ndist-1
       if ( u3i(k,n) .ge. uitv(0,ii) .and. u3i(k,n) .lt. uitv(1,ii) ) then
         ncom (1,4,k,ii) = ncom (1,4,k,ii) + 1.
         tend (1,2,k,ii) = tend (1,2,k,ii) + ut3i(k,n)
         shear(1,2,k,ii) = shear(1,2,k,ii) + sh3i(k,n)
       end if
     end do
   end if
 end do
 end do

 ncom@_FillValue = 0.

 tend (:,0:1,:,:) = tend (:,0:1,:,:)/ncom(:,0:1,:,:)
 gwd  (:,:  ,:,:) = gwd  (:,:  ,:,:)/ncom(:,0:3,:,:)
 shear(:,0:1,:,:) = shear(:,0:1,:,:)/ncom(:,0:1,:,:)

 tend (:,2,:,:) = tend (:,2,:,:)/ncom(:,4,:,:)
 shear(:,2,:,:) = shear(:,2,:,:)/ncom(:,4,:,:)

;--------------------------------------- plot
 ymin = (/-0.1,0.1/)
 ymax = (/1.2,-1.2/)
 xmin = (/min(uitv(0,:)),max(uitv(1,:))/)
 xmax = (/max(uitv(1,:)),min(uitv(0,:))/)
 wi = 32
 he = 24
 xf = 10 + (/0,1,0/)*39.5
 yf = 10 + (/0,1,2/)*28

 mstr = new((/nz,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"Positive","Negative"/)+"-tendency phase"
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "U (m s:S:-1:N:)"
 ystr = (/"GWD (m s:S:-1:N: / day)",""/)
 lbstr = (/"BGW (CTL)","BGW (CC)","CGW (CC)"/)

 i2 = (/0,2,3/)
 fcol = (/6,4,2/)
 xx0 = (/ (/-1,0,0/), (/0,-1,-1/) /)*itv*0.3

 do k=1, 1 ;0, nz-1
 do is=0, 1

   var0 = array_append_record( tend(is,0:2,k,:), gwd(is,0:1,k,:), 0 )
   var1 = array_append_record( dura(is,0:2,k,:)*(-1)^is * 0.05, shear(is,0:2,k,:)*0.05, 0 )
   var  = array_append_record( var0, var1, 0 )
   plot = xylineplot(wks,udist,var,"D",(/16,0,0,0,0,16,0,0,16,0,0/),"T",1.0,"C",(/1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1/))
   xyaddline(plot,udist*2,0.,"D",0,"T",0.5,"C",colind)

   dum1a = new((/3,ndist/),graphic)
   dum1b = new((/3,ndist/),graphic)
   do ii=0, ndist-1
     yyp = (/0.,0.,0.,0.,0./)
     do i=0, 2
       if ( .not. ismissing(gwd(is,i2(i),k,ii)) ) then
         tmp = 0.
         if (i .eq. 2) then
           if ( gwd(is,i2(i),k,ii) .ge. 0. ) then
             tmp = max(yyp)
           else
             tmp = min(yyp)
           end if
         end if
         xxp = xx0(is,i) + (/0,0,1,1,0/)*(itv*0.3) + udist(ii)
         yyp = tmp + (/0,1,1,0,0/)*gwd(is,i2(i),k,ii)
         pgres                  = True
         pgres@gsLineThicknessF = 0.8
         pgres@gsFillColor      = fcol(i)
         dum1a(i,ii) = gsn_add_polygon (wks,plot,xxp,yyp,pgres)
         dum1b(i,ii) = gsn_add_polyline(wks,plot,xxp,yyp,pgres)
       end if
     end do
   end do

   set_size(plot,xf(is),yf(k),wi,he)
   axis_range(plot,xmin(is),xmax(is),ymin(is),ymax(is),False,False)
   axis_str(plot,mstr(0,is),xstr(2),ystr(is))
   setvalues plot
    "tmXBMode" : "Explicit"
    "tmXBValues" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
    "tmXBLabels" : array_append_record( uitv(0,:), max(uitv(1,:)), 0 )
   end setvalues

;   dr_txt(wks,p(k)+" hPa",xf(0)+8,yf(k)+3,22,8,0.)

   draw(plot)

   delete(dum1a)
   delete(dum1b)

 end do
 end do

 do i=0, 2
   xxp = ( xf(0)+28 + (/0,0,1,1,0/)*wi(0)*0.5 )/p_xpg
   yyp = ( yf(1)+94+i*6 + (/0,1,1,0,0/)*he*0.1    )/p_ypg
   pgres                  = True
   pgres@gsLineThicknessF = 0.8
   pgres@gsFillColor      = fcol(i)
   gsn_polygon_ndc (wks,xxp,1.-yyp,pgres)
   gsn_polyline_ndc(wks,xxp,1.-yyp,pgres)
   txtres               = True
   txtres@txFont        = 21
   txtres@txFontHeightF = 0.005
   txtres@txJust        = "CenterLeft"
   gsn_text_ndc(wks,lbstr(i),xxp(2)+0.005,1.-0.5*(yyp(0)+yyp(1)),txtres)
 end do

 lg = lgplot(wks,(/"dU/dt (CTL)","dU/dt (CC)"/),"D",(/16,0,0,0,0,16,0,0,16,0,0/),"T",1.0,"C",(/1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1/),2)
 set_size_lg(lg,xf(0)+17.0,yf(1)+0.7,wi*0.28,5.2)
 setvalues lg
  "lgLabelOffsetF" : 0.2
  "lgLineDashSegLenF"  : 0.030
  "lgLabelFontHeightF" : 0.005
 end setvalues
 draw(lg)


 frame(wks)

end



