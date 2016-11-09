load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.25

begin
;--------------------------------------- parameters
 cexp    = "CGW"
 expname = "tc"
 ppid    = "j"
 varname = "u"

 z = (/26.9,21.0/)  ; (/10,15,20,30,40,50/) ~ (/31.9,29.3,26.9,24.5,22.1,21.0/)
 p = (/20,50/)

 y1 = 1950 ;2058 ;1960
 y2 = 2100 ;2097 ;1999
 lat0 = 5.
 opt_ss = 3  ; 0, 1, 2, 3 : no, mean, linear, and quadratic seasonal cycle

 wks = new(2,graphic)
 wks(0) = cr_wks("fig/peramp1_t0-z1-"+expname,600,0)
 wks(1) = cr_wks("fig/peramp2_t0-z1-"+expname,600,0)

;--------------------------------------- read data

 fdir  = "$DATD/AOL60"+cexp
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$
; z = tmp&$tmp!1$ / 1.e3
 nz = dimsizes(z)
 t = tmp&$tmp!0$
 t_char = tochar(t@units)
 if ( t_char(0:10) .ne. "days since " ) then
   print(t@units)
   print("Check this. It should be 'days since ...'.")
   exit
 end if
 y0i = toint(tostring(t_char(11:14)))
 m0i = toint(tostring(t_char(16:17)))
 d0i = toint(tostring(t_char(19:20)))
 print("date0 : "+y0i+" "+m0i+" "+d0i)

 t = y0i + (m0i-1)/12. + (d0i-1)/360. + t/360.
 nt = dimsizes(t)
 print("data : "+t(0)+" - "+t(nt-1))

 u = transpose( dim_avg( tmp(:,{z*1.e3},{-lat0:lat0}) ) )

 delete(tmp)

 u!1 = "t"
 u&t = t

;--------------------------------------- process 0
 if (opt_ss .ne. 0) then
   udes = u
   if (opt_ss .eq. 1) then  ; mean seasonal cycle
     uss = new((/nz,12/),float)
     um = dim_avg(u)
     do n=0, 11
       itmps = ispan(n,nt-1,12)
       uss(:,n) = dim_avg(u(:,itmps))
       udes(:,itmps) = u(:,itmps) - conform_dims(dimsizes(u(:,itmps)),uss(:,n)-um,0)
       delete(itmps)
     end do
   else
     if (opt_ss .eq. 2) then  ; seasonal cycle having a linear trend
       do n=0, 11
         itmps = ispan(n,nt-1,12)
         udes(:,itmps) = dtrend(u(:,itmps),False)
         delete(itmps)
       end do
       uss = u - udes
       utmp = uss - dtrend(uss,False)
     else if (opt_ss .eq. 3) then  ; seasonal cycle having a quadratic trend
       do n=0, 11
         itmps = ispan(n,nt-1,12)
         udes(:,itmps) = dtrend_quadratic(u(:,itmps),0)
         delete(itmps)
       end do
       uss = u - udes
       utmp = uss - dtrend_quadratic(uss,0)
     end if
     end if
     if (opt_ss .eq. 4) then
       nra = 41
       uss = u
       utmp = runave(u,nra,0)
       do n=0, 11
         itmps = ispan(n,nt-1,12)
         uss(:,itmps) = runave(u(:,itmps),nra,0)
;         uss(:,itmps(0:
         delete(itmps)
       end do
       udes = u - uss
     end if
     udes = udes + utmp
   end if 

   u = udes

 end if


; uu = u
 uu = runave_Wrap(u,3,0)

 nd = new(nz,integer)
 nd = 0
 t0 = new((/nz,nt/),float,1.e20)
 peramp = new((/2,3,nz,nt/),float,1.e20)

 do k=0, nz-1
   tmpi = 0
   do n=1, nt-1
   if ( .not. ( ismissing(uu(k,n-1)) .or. ismissing(uu(k,n)) ) ) then
     if ( uu(k,n-1) .ge. 0. .and. uu(k,n) .lt. 0. .and. tmpi .ne. 0 ) then
       tmp1 = 0.5*(t(n-1) + t(n))
       t0(k,nd(k)) = tmp1
       peramp(0,0,k,nd(k)) = (tmp1 - tmp0)*12.
       peramp(1,0,k,nd(k)) = max(abs( uu(k,{tmp0:tmp1}) ))
     end if
     if ( uu(k,n-1) .lt. 0. .and. uu(k,n) .ge. 0. ) then
       tmp0 = 0.5*(t(n-1) + t(n))
       if ( tmpi .ne. 0 ) then
         peramp(0,1,k,nd(k)) = (tmp0 - tmp1)*12
         peramp(1,1,k,nd(k)) = max(abs( uu(k,{tmp1:tmp0}) ))
         peramp(0,2,k,nd(k)) = sum( peramp(0,0:1,k,nd(k)) )
         peramp(1,2,k,nd(k)) = 0.5*sum( peramp(1,0:1,k,nd(k)) )
         nd(k) = nd(k) + 1
       end if
       tmpi = 1
     end if
   end if
   end do
 end do

;--------------------------------------- cn level
 lcol = (/1,2/)
 ldsh = 0
 ltck = 1.8

;--------------------------------------- plot
 wi = 150
 he = 40
 xf = 10
 yf = 10 + (/0,1,2/)*(he+15)

 mstr = (/ (/"Duration of westerly phases","Duration of easterly phases","QBO period"/), \
           (/"Magnitude of westerly phases","Magnitude of easterly phases","QBO magnitude"/) /)
 xstr = (/"","","year"/)
 ystr = (/"length (month)","(m s:S:-1:N:)"/)

 ymin = (/ (/0 ,6 ,10/), (/0 , 5,0/) /)
 ymax = (/ (/18,18,32/), (/25,35,28/) /)

do iw=0, 1

do iv=0, 2

 plot = xylineplot(wks(iw),t0(0,:),peramp(iw,iv,0,:),"D",ldsh,"T",ltck,"C",-1)
 xyaddline(plot,t0(0,:),runave(peramp(iw,iv,0,:),7,0),"D",ldsh,"T",ltck,"C",lcol(0))
 xyaddmarker(plot,(/1,-1/),"S",3.,"T",1.,"C",lcol)
 getvalues plot
  "xyCoordData" : xydat
 end getvalues
 setvalues xydat
  "caYMissingV" : 1.e20
 end setvalues

 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,y1,y2,ymin(iw,iv),ymax(iw,iv),False,False)
 axis_str(plot,mstr(iw,iv),xstr(iv),ystr(iw))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 10
  "tmXBTickStartF"    : 1950
  "tmXBMinorPerMajor" : 4
 end setvalues

 draw(plot)

 if (min(peramp(iw,iv,0,:)) .lt. ymin(iw,iv)) then
   print("Caution. min( peramp("+iw+","+iv+",0,:) ) = "+  \
         min(peramp(iw,iv,0,:))+" < "+ymin(iw,iv)+"  !!")
 end if
 if (max(peramp(iw,iv,0,:)) .gt. ymax(iw,iv)) then
   print("Caution. max( peramp("+iw+","+iv+",0,:) ) = "+  \
         max(peramp(iw,iv,0,:))+" > "+ymax(iw,iv)+"  !!")
 end if

end do

 lg = lgplot(wks(iw),p+" hPa","D",0,"T",1.8,"C",lcol(0),2)
 set_size_lg(lg,xf,yf(0),15,12)
; draw(lg)

 frame(wks(iw))

end do

end

