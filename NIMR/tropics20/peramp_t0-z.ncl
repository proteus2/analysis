load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 ppid    = (/"j" /)
 varname = "u"
 expname2= "sonde"

 lat0 = 15
 z = (/26.9,23.3,21.0/)  ; (/10,20,30,50/) ~ (/30.6,26.9,23.3,21.0/)
 p = (/20,30,50/)
 u0 = 6.

 y1 = 1952
 y2 = 2013

 wks = new(2,graphic)
 wks(0) = cr_wks("fig/peramp1_t0-z-"+expname(0)+"_"+lat0,600,0)
 wks(1) = cr_wks("fig/peramp2_t0-z-"+expname(0)+"_"+lat0,600,0)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new((/ne,nt/),string)
 fdir  = "../dat/xm"
 do ie=0, ne-1
   n = 0
   do iy=y1, y2
   do im=1, 12
     fname(ie,n) = fdir+"/"+expname(ie)+"/u_xm.uan"+expname(ie)+"a.p"+ppid(ie)+"_"+iy+sprinti("%0.2i",im)+".nc"
     n = n + 1
   end do
   end do
 end do
 fexist = isfilepresent(fname)

 n1 = new(ne,integer)
 n2 = new(ne,integer)
 n1 = -1
 n2 = nt-1
 do ie=0, ne-1
   do n=0, nt-1
     if ( fexist(ie,n) ) then
       if (n1(ie) .eq. -1) then
         n1(ie) = n
       end if
     else
       if ( n1(ie) .ne. -1 .and. n2(ie) .eq. nt-1 ) then
         n2(ie) = n - 1
       end if
     end if
   end do
   if (n1(ie) .eq. -1) then
     n2(ie) = -1
   end if
 end do
print(n1/12.+"  "+n2/12.)

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 nz = dimsizes(z)

 u = new((/ne+1,nz,nt/),float,1.e20)
 do ie=0, ne-1
   if (n1(ie) .ne. -1) then
print(expname(ie)+"")
     f1 = addfiles(fname(ie,n1(ie):n2(ie)),"r")
     ListSetType(f1,"cat")
     tmp = dim_avg( f1[:]->$varname(ie)$(:,{z*1.e3},{-lat0:lat0}) )
     u(ie,:,n1(ie):n2(ie)) = transpose( tmp )
     delete(tmp)
   end if
 end do

 fdir  = "/data17/kyh/qbo_obs"
 fname2 = fdir+"/qbo_"+expname2+".nc"
 f2 = addfile(fname2,"r")

 uobs = f2->u({p},:)
 tobs = f2->t
 ntobs = dimsizes(tobs)
 n = floattointeger( (tobs(0) - t(0))*12. )
 u(ne,:,n:n+ntobs-1) = (/ uobs(:,:) /)

 u = u + u0

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

 nd = new((/ne+1,nz/),integer)
 nd = 0
 t0  = new((/ne+1,nz,nt/),float,1.e20)
 peramp = new((/2,3,ne+1,nz,nt/),float,1.e20)

 do ie=0, ne
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
 ldsh = (/0,2/)
 ltck = (/1.8,0.8/)

;--------------------------------------- plot
 wi = 150
 he = 40
 xf = 10
 yf = 10 + (/0,1,2/)*(he+15)

 mstr = (/ (/"Duration of westerly phases","Duration of easterly phases","QBO period"/), \
           (/"Magnitude of westerly phases","Magnitude of easterly phases","QBO magnitude"/) /)
 xstr = (/"","","year"/)
 ystr = (/"(month)","(m s:S:-1:N:)"/)
 estr =  (/cexp(0),"Obs"/)

 ymin = (/ (/6 ,6 ,19/), (/0 , 0,0 /) /)
 ymax = (/ (/25,25,36/), (/40,40,40/) /)

do iw=0, 1

do iv=0, 2
do ie=0, ne

 plot = xylineplot(wks(iw),t0(ie,:,:),peramp(iw,iv,ie,:,:),"D",ldsh(ie),"T",ltck(ie),"C",lcol)
; xyaddline             (plot,t0(ie,:,:),peramp(iv,ie,:,:),"D",ldsh(ie),"T",ltck(ie),"C",lcol)
 xyaddmarker            (plot,1,"S",1.,"T",1.,"C",1)
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

 lg = lgplot(wks(iw),p+" hPa","D",0,"T",1.8,"C",lcol,2)
 set_size_lg(lg,xf,yf(0),15,12)
 draw(lg)

 frame(wks(iw))

end do

end

