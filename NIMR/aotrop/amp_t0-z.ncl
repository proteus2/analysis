load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = (/"CTL"/)
 expname = (/"tk"/)
 ppid    = (/"j" /)
 varname = "u"
 expname2= "sonde"

 lat0 = 5.
 z = (/26.9,23.3,21.0,18.59/)  ; (/10,20,30,50/) ~ (/30.6,26.9,23.3,21.0/)
 p = (/20,30,50,70/)

 y1 = 1950
 y2 = 2014

 wks = cr_wks("fig/amp_t0-z-"+expname(0),600,0)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new((/ne,nt/),string)
 fdir  = "../dat/AOL60"+cexp+"/xm"
 do ie=0, ne-1
   n = 0
   do iy=y1, y2
   do im=1, 12
     fname(ie,n) = fdir(ie)+"/xm.uan"+expname(ie)+"a.p"+ppid(ie)+"_"+iy+sprinti("%0.2i",im)+".nc"
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

 fdir  = "/data18/kyh/qbo_obs"
 fname2 = fdir+"/qbo_"+expname2+".nc"
 f2 = addfile(fname2,"r")

 uobs = f2->u({p},:)
 tobs = f2->t
 ntobs = dimsizes(tobs)
 n = floattointeger( (tobs(0) - t(0))*12. )
 u(ne,:,n:n+ntobs-1) = (/ uobs(:,:) /)

 u!2 = "t"
 u&t = t

;--------------------------------------- process 0
   uss = new((/2,nz,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n+(1953-y1)*12,(2010-y1)*12-1,12)))
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n+(1953-y1)*12,(2013-y1)*12-1,12)))  ; 1953.01 - 2012.12
   end do

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

 udes!2 = "t"
 udes&t = t

; uu = u
; uu = runave_Wrap(u,5,0)
 uu = runave_Wrap(udes,5,0)

 nd = new((/ne+1,nz/),integer)
 nd = 0
 t0  = new((/ne+1,nz,nt/),float,1.e20)
 amp3 = new((/ne+1,nz,nt/),float,1.e20)
 amp = new((/2,ne+1,nz,nt/),float,1.e20)
 avp = new((/4,ne+1,nz,nt/),float,1.e20)
 amp3 = 0.
 amp = 0.
 avp = 0.

 k0 = ind(p .eq. 30)
print(p(k0))

 do ie=0, ne
 do k=0, nz-1
   tmpi = 0
;   kk = k
   kk = k0
   do n=1, nt-1
   if ( .not. ( ismissing(uu(ie,kk,n-1)) .or. ismissing(uu(ie,kk,n)) ) ) then
     if ( uu(ie,kk,n-1) .lt. 0. .and. uu(ie,kk,n) .ge. 0. .and. tmpi .ne. 0 ) then
       tmp1 = 0.5*(t(n-1) + t(n))
       t0(ie,k,nd(ie,k)) = tmp1
       amp(0,ie,k,nd(ie,k)) = sum(uu(ie,k,{tmp0:tmp1})^2)
       avp(0,ie,k,nd(ie,k)) = (tmp1-tmp0)*12.
       avp(2,ie,k,nd(ie,k)) = sum(uu(ie,k,{tmp0:tmp1}))
     end if
     if ( uu(ie,kk,n-1) .ge. 0. .and. uu(ie,kk,n) .lt. 0. ) then
       tmp0 = 0.5*(t(n-1) + t(n))
       if ( tmpi .ne. 0 ) then
         amp(1,ie,k,nd(ie,k)) = sum(uu(ie,k,{tmp1:tmp0})^2)
         avp(1,ie,k,nd(ie,k)) = (tmp0-tmp1)*12.
         avp(3,ie,k,nd(ie,k)) = sum(uu(ie,k,{tmp1:tmp0}))
         nd(ie,k) = nd(ie,k) + 1
       end if
       tmpi = 1
     end if
   end if
   end do
 end do
 end do

 do ie=0, ne
 do k=0, nz-1
   do n=1, nd(ie,k)-2
     amp3(ie,k,n) = sqrt( sum(amp(0:1,ie,k,n-1:n+1))/sum(avp(0:1,ie,k,n-1:n+1)) - \
                          (sum(avp(2:3,ie,k,n-1:n+1))/sum(avp(0:1,ie,k,n-1:n+1)))^2 )
   end do
 end do
 end do
 amp3 = amp3*sqrt(2.)

 amp3 = where( amp3 .eq. 0., 1.e20, amp3 )

;--------------------------------------- cn level
 lcol = (/3,1,4,2/)
 ldsh = (/0,2/)
 ltck = (/1.8,1.2/)

;--------------------------------------- plot
 wi = 63
 he = 90
 xf = 10
 yf = 10 + (/0,1,2/)*(he+13)

 mstr = "" ;(/"amplitude"/)
 xstr = (/"","","year"/)
 ystr = "amplitude (m s:S:-1:N:)"
 estr =  (/cexp(0),"Obs"/)

 ymin = (/3 , 0,0 /)
 ymax = (/29,40,40/)

do iv=0, 0 ; dummy
do ie=0, ne

 plot = xylineplot(wks,t0(ie,:,:),amp3(ie,:,:),"D",ldsh(ie),"T",ltck(ie),"C",lcol)
; xyaddline             (plot,t0(ie,:,:),amp3(ie,:,:),"D",ldsh(ie),"T",ltck(ie),"C",lcol)
 xyaddmarker            (plot,1,"S",1.,"T",1.,"C",1)
 getvalues plot
  "xyCoordData" : xydat
 end getvalues
 setvalues xydat
  "caYMissingV" : 1.e20
 end setvalues

 set_size(plot,xf,yf(iv),wi,he)
 axis_range(plot,1950,2013,ymin(iv),ymax(iv),False,False)
 axis_str(plot,mstr(iv),xstr(iv),ystr)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 10
  "tmXBTickStartF"    : 1950
  "tmXBMinorPerMajor" : 9
 end setvalues

 draw(plot)

end do
end do

do k=0, nz-1
 lg = lgplot(wks,p(k)+" hPa","D",0,"T",1.8,"C",lcol(k),2)
 set_size_lg(lg,xf+1+k*16.2,yf(0)+1,6,4)
 draw(lg)
end do

 frame(wks)

end

