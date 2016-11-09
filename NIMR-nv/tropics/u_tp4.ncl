load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"tr","ts","tt","ui","uj"/)
 ppid    = (/"j" ,"j" ,"j" ,"j" ,"j" /)
 varname = "u"+(/"_2","_2","_2","_2","_2"/)

 y1 = 1942
 y2 = 1965
 lat0 = 5.

 wks = cr_wks("fig/u_tp4-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fname = new((/ne,nt/),string)
 fdir  = "../dat/xm"
 do ie=0, ne-1
   n = 0
   do iy=y1, y2
   do im=1, 12
     umtail5 = systemfunc("ex_date_umstr.sh "+iy+" "+im)
     fname(ie,n) = fdir+"/"+expname(ie)+"/u-p_xm.uan"+expname(ie)+"a.p"+ppid(ie)+umtail5+".nc"
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

 ; get z
 ie0 = maxind(n1(:ne-3))
 f0 = addfile(fname(ie0,n1(ie0)),"r")
 tmp = f0->$varname(ie0)$
 zname = tmp!1
 z = 7.0*log(1.e3/tmp&$zname$)
 nz = dimsizes(z)
 delete(tmp)

 ; get z9
 ie0 = ne-2 + maxind(n1(ne-2:))
 f0 = addfile(fname(ie0,n1(ie0)),"r")
 tmp = f0->$varname(ie0)$
 zname = tmp!1
 z9 = 7.0*log(1.e3/tmp&$zname$)
 nz9 = dimsizes(z9)
 delete(tmp)

 u = new((/ne-2,nz,nt/),float,1.e20)
 u9 = new((/2,nz9,nt/),float,1.e20)
 do ie=0, ne-1
   if (n1(ie) .ne. -1) then
print(expname(ie)+"")
     f1 = addfiles(fname(ie,n1(ie):n2(ie)),"r")
     ListSetType(f1,"cat")
     tmp = dim_avg( f1[:]->$varname(ie)$(:,:,{-lat0:lat0}) )
     if (ie .lt. ne-2) then
       u(ie,:,n1(ie):n2(ie)) = transpose( tmp )
     else
       u9(ie-ne+2,:,n1(ie):n2(ie)) = transpose( tmp )
     end if
     delete(tmp)
   end if
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72*2
 he = 20
 xf = 10 + (/0,0,0,0,0,0,0,0/)*78
 yf = 10 + (/1,2,3,0,4/)*32

 mstr  = new(ne,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)
 lbopt = new(ne,string)

 mstr = (/"N1 from 1952/CTL","N2","N3","CTL","CTL99"/)
 xstr(:) = ""
 xstr(ne-1) = "year"
 ystr(:) = "height (km)"
 lbopt(:) = ""
 lbopt(ne-1) = "(m s:S:-1:N:)"

do ie=0, ne-1

 if ( ie .lt. ne-2 ) then
   contour1 = cnshadeplot(wks,t,z ,u (ie,:,:),cnlev,cnfil,lbopt(ie))
 else
   contour1 = cnshadeplot(wks,t,z9,u9(ie-ne+2,:,:),cnlev,cnfil,lbopt(ie))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(ie),yf(ie),wi,he)
 axis_range(contour1,1942,1966,16,45,False,False)
 axis_str(contour1,mstr(ie),xstr(ie),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 1942
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

end do


 frame(wks)

end



