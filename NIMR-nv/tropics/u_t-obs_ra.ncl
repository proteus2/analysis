load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"sonde","profiler","era-int","merra"/)

 lat0 = 5.
 p = (/50,30,20,10/)

 wks = cr_wks("fig/u_t-obs_ra",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)

 fdir  = "/data15/kyh/qbo/obs"
 fname1 = fdir+"/qbo_"+expname(0:1)+".nc"
 f1 = addfiles(fname1,"r")

 u1 = f1[0]->u({p},:)
 u2 = f1[1]->u({p},:)

 t1 = f1[0]->t
 t2 = f1[1]->t

 nt1 = dimsizes(t1)
 nt2 = dimsizes(t2)

 fdir  = "../dat"
 fname3 = fdir+"/"+expname(2:3)+".u_ypt.1979-2010.nc"
 f3 = addfiles(fname3,"r")

 u3 = transpose( doubletofloat(dim_avg( f3[0]->u(:,{p},{-lat0:lat0}) )) )
 u4 = transpose( dim_avg( f3[1]->u(:,{p},{-lat0:lat0}) ) )

 nt3 = (2010-1979+1)*12
 t3 = 1979. + (fspan(1,nt3,nt3)-0.5)/12.
 nt4 = nt3
 t4 = t3

 t_stt = min( (/t1(0),t2(0),t3(0),t4(0)/) )
 t_end = max( (/t1(nt1-1),t2(nt2-1),t3(nt3-1),t4(nt4-1)/) )
 nt = floattointeger( (t_end - t_stt)*12. ) + 1
 t = fspan(t_stt,t_end,nt)

 nz = dimsizes(p)

 u = new((/ne,nz,nt/),float,1.e20)
 n = floattointeger( (t1(0) - t_stt)*12. )
 u(0,:,n:n+nt1-1) = (/ u1(:,:) /)
 n = floattointeger( (t2(0) - t_stt)*12. )
 u(1,:,n:n+nt2-1) = (/ u2(:,:) /)
 n = floattointeger( (t3(0) - t_stt)*12. )
 u(2,:,n:n+nt3-1) = (/ u3(:,:) /)
 n = floattointeger( (t4(0) - t_stt)*12. )
 u(3,:,n:n+nt4-1) = (/ u4(:,:) /)

 uu = runave(u,5,0)

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
 per = new((/3,ne,nz,nt/),float,1.e20)

 do ie=0, ne-1
 do k=0, nz-1
   tmpi = 0
   do n=1, nt-1
   if ( .not. ( ismissing(u(ie,k,n-1)) .or. ismissing(u(ie,k,n)) ) ) then
     if ( u(ie,k,n-1) .lt. 0. .and. u(ie,k,n) .ge. 0. ) then
       tmp0 = 0.5*(t(n-1) + t(n))
       if ( tmpi .ne. 0 ) then
         per(1,ie,k,nd(ie,k)) = (tmp0 - tmp1)*12
         per(2,ie,k,nd(ie,k)) = sum( per(0:1,ie,k,nd(ie,k)) )
         nd(ie,k) = nd(ie,k) + 1
       end if
       tmpi = 1
     end if
     if ( u(ie,k,n-1) .ge. 0. .and. u(ie,k,n) .lt. 0. .and. tmpi .ne. 0 ) then
       tmp1 = 0.5*(t(n-1) + t(n))
       t0(ie,k,nd(ie,k)) = tmp1
       per(0,ie,k,nd(ie,k)) = (tmp1 - tmp0)*12.
     end if
   end if
   end do
 end do
 end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 lcol = (/1,1,3,4/)

;--------------------------------------- plot
 wi = 150
 he = 30
 xf = 10
 yf = 10

 mstr  = new(ne*3,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)
 lbopt = new(ne,string)

 mstr = ""
; mstr(0:3) = (/"Sonde data (three stations)","","",""/)
; mstr(4:7) = mstr(0:3)+(/", ERA-Interim, MERRA","","",""/)
; mstr(8:11) = mstr(0:3)+(/", Profiler data (Singapore), ERA-Interim, MERRA","","",""/)
 xstr = ""
 ystr = "(month)"

do k=0, nz-1

 plot = xylineplot(wks,t,u(:,k,:),"D",0,"T",0.5,"C",lcol)
 xyaddline       (plot,t,uu(:,k,:),"D",0,"T",0.5,"C",lcol)
 xyaddline       (plot,t,0,"D",0,"T",0.5,"C",colind)
 getvalues plot
  "xyCoordData" : xydat
 end getvalues
 setvalues xydat
  "caYMissingV" : 1.e20
 end setvalues

 set_size(plot,xf,yf+k*(he+10),wi,he)
 axis_range(plot,1953,2013,-999,999,False,False)
 axis_str(plot,mstr,xstr,ystr)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 5
  "tmXBTickStartF"    : 1950
  "tmXBMinorPerMajor" : 4
 end setvalues

 draw(plot)

end do


 frame(wks)

end

