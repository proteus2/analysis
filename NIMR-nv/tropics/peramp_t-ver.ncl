load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"sonde","era-int","merra"/)

 lat0 = 5.
 p = (/20,30,50/)

 wks = new(2,graphic)
 wks(0) = cr_wks("fig/peramp1_t-ver",600,0)
 wks(1) = cr_wks("fig/peramp2_t-ver",600,0)
 colind = NhlNewColor(wks(0),0.4,0.4,0.4)
 colind = NhlNewColor(wks(1),0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)

 fdir  = "/data18/kyh/qbo_obs"
 fname1 = fdir+"/qbo_"+expname(0)+".nc"
 f1 = addfiles(fname1,"r")

 u1 = f1[0]->u({p},:)

 t1 = f1[0]->t

 nt1 = dimsizes(t1)

 fdir  = "../dat/L60CGW"
 fname3 = fdir+"/"+expname(1:2)+".u_ypt.1979-2010.nc"
 f3 = addfiles(fname3,"r")

 u3 = transpose( doubletofloat(dim_avg( f3[0]->u(:,{p},{-lat0:lat0}) )) )
 u4 = transpose( dim_avg( f3[1]->u(:,{p},{-lat0:lat0}) ) )

 nt3 = (2010-1979+1)*12
 t3 = 1979. + (fspan(1,nt3,nt3)-0.5)/12.
 nt4 = nt3
 t4 = t3

 t_stt = min( (/t1(0),t3(0),t4(0)/) )
 t_end = max( (/t1(nt1-1),t3(nt3-1),t4(nt4-1)/) )
 nt = floattointeger( (t_end - t_stt)*12. ) + 1
 t = fspan(t_stt,t_end,nt)

 nz = dimsizes(p)

 u = new((/ne,nz,nt/),float)
 n = floattointeger( (t1(0) - t_stt)*12. )
 u(0,:,n:n+nt1-1) = (/ u1(:,:) /)
 n = floattointeger( (t3(0) - t_stt)*12. )
 u(1,:,n:n+nt3-1) = (/ u3(:,:) /)
 n = floattointeger( (t4(0) - t_stt)*12. )
 u(2,:,n:n+nt4-1) = (/ u4(:,:) /)

 u!2 = "t"
 u&t = t

; uu = u
 uu = runave_Wrap(u,5,0)

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

 nd = new((/ne,nz/),integer)
 nd = 0
 t0  = new((/ne,nz,nt/),float,1.e20)
 peramp = new((/2,3,ne,nz,nt/),float,1.e20)

 do ie=0, ne-1
 do k=0, nz-1
   tmpi = 0
   do n=1, nt-1
   if ( .not. ( ismissing(uu(ie,k,n-1)) .or. ismissing(uu(ie,k,n)) ) ) then
     if ( uu(ie,k,n-1) .ge. 0. .and. uu(ie,k,n) .lt. 0. .and. tmpi .ne. 0 ) then
       tmp1 = 0.5*(t(n-1) + t(n))
       t0(ie,k,nd(ie,k)) = tmp1
       peramp(0,0,ie,k,nd(ie,k)) = (tmp1 - tmp0)*12.
       peramp(1,0,ie,k,nd(ie,k)) = max(abs( uu(ie,k,{tmp0:tmp1}) ))
     end if
     if ( uu(ie,k,n-1) .lt. 0. .and. uu(ie,k,n) .ge. 0. ) then
       tmp0 = 0.5*(t(n-1) + t(n))
       if ( tmpi .ne. 0 ) then
         peramp(0,1,ie,k,nd(ie,k)) = (tmp0 - tmp1)*12
         peramp(1,1,ie,k,nd(ie,k)) = max(abs( uu(ie,k,{tmp1:tmp0}) ))
         peramp(0,2,ie,k,nd(ie,k)) = sum( peramp(0,0:1,ie,k,nd(ie,k)) )
         peramp(1,2,ie,k,nd(ie,k)) = 0.5*sum( peramp(1,0:1,ie,k,nd(ie,k)) )
         nd(ie,k) = nd(ie,k) + 1
       end if
       tmpi = 1
     end if
   end if
   end do
 end do
 end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 lcol = (/3,1,2/)
 ldsh = (/0,1,2/)

;--------------------------------------- plot
 wi = 150
 he = 40
 xf = 10
 yf = 10 + (/0,1,2/)*(he+15)

 mstr = (/ (/"Duration of westerly phases","Duration of easterly phases","QBO period"/), \
           (/"Magnitude of westerly phases","Magnitude of easterly phases","QBO magnitude"/) /)
 xstr = (/"","","year"/)
 ystr = (/"(month)","(m s:S:-1:N:)"/)

 ymin = (/ (/6 ,6 ,19/), (/0 ,10,5 /) /)
 ymax = (/ (/25,25,36/), (/30,40,35/) /)

do iw=0, 1

do iv=0, 2
do ie=0, ne-1

 plot = xylineplot(wks(iw),t0(ie,:,:),peramp(iw,iv,ie,:,:),"D",ldsh(ie),"T",0.8,"C",lcol)
; xyaddline       (plot,t0(ie,:,:),peramp(iv,ie,:,:),"D",ldsh(ie),"T",0.8,"C",lcol)
 xyaddmarker      (plot,1,"S",1.,"T",1.,"C",1)
 getvalues plot
  "xyCoordData" : xydat
 end getvalues
 setvalues xydat
  "caYMissingV" : 1.e20
 end setvalues

 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,1953,2013,ymin(iw,iv),ymax(iw,iv),False,False)
 axis_str(plot,mstr(iw,iv),xstr(iv),ystr(iw))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 5
  "tmXBTickStartF"    : 1950
  "tmXBMinorPerMajor" : 4
 end setvalues

 draw(plot)

end do
end do

 lg = lgplot(wks(iw),p+" hPa","D",0,"T",0.8,"C",lcol,2)
 set_size_lg(lg,xf,yf(0),15,12)
 draw(lg)

 frame(wks(iw))

end do

end

