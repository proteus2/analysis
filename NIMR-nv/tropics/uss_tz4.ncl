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
 varname = "u_"+(/1,1,1,1,1/)

 y1 = 1942
 y2 = 1965
 lat0 = 5.

 wks = cr_wks("fig/uss_tz4-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 wks2 = cr_wks("fig/udes_tz4-"+expname(0),600,0)
 gsn_define_colormap(wks2,"StepSeq25")
 colind = NhlNewColor(wks2,0.4,0.4,0.4)

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
     fname(ie,n) = fdir+"/"+expname(ie)+"/u_xm.uan"+expname(ie)+"a.p"+ppid(ie)+umtail5+".nc"
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
 ie0 = maxind(n1(:ne-2))
 f0 = addfile(fname(ie0,n1(ie0)),"r")
 tmp = f0->$varname(ie0)$
 zname = tmp!1
 z = tmp&$zname$/1.e3
 nz = dimsizes(z)
 delete(tmp)

 ; get z9
 ie0 = ne-1
 f0 = addfile(fname(ie0,n1(ie0)),"r")
 tmp = f0->$varname(ie0)$
 zname = tmp!1
 z9 = tmp&$zname$/1.e3
 nz9 = dimsizes(z9)
 delete(tmp)

 u = new((/ne-1,nz,nt/),float,1.e20)
 u9 = new((/1,nz9,nt/),float,1.e20)
 do ie=0, ne-1
   if (n1(ie) .ne. -1) then
print(expname(ie)+"")
     f1 = addfiles(fname(ie,n1(ie):n2(ie)),"r")
     ListSetType(f1,"cat")
     tmp = dim_avg( f1[:]->$varname(ie)$(:,:,{-lat0:lat0}) )
     if (ie .lt. ne-1) then
       u(ie,:,n1(ie):n2(ie)) = transpose( tmp )
     else
       u9(ie-ne+1,:,n1(ie):n2(ie)) = transpose( tmp )
     end if
     delete(tmp)
   end if
 end do

 t = 1942 + (fspan(1,nt,nt)-0.5)/12.

 uss = new((/ne-1,nz,14/),float)
 uss9 = new((/1,nz9,14/),float)
 udes = u
 udes9 = u9
 do ie=0, ne-1
   if ( n1(ie) .ne. -1 .and. n2(ie)-n1(ie)+1 .ge. 24 ) then
     if (ie .lt. ne-1) then
       do n=1, 12
         uss(ie,:,n) = dim_avg(u(ie,:,ispan(n1(ie)+n-1,n2(ie),12)))
       end do
       uss(ie,:,0 ) = uss(ie,:,12)
       uss(ie,:,13) = uss(ie,:,1 )
       do n=n1(ie), n2(ie)
         udes(ie,:,n) = udes(ie,:,n) - uss(ie,:,mod(n-n1(ie)+1,12))
       end do
     else
       do n=1, 12
         uss9(ie-ne+1,:,n) = dim_avg(u9(ie-ne+1,:,ispan(n1(ie)+n-1,n2(ie),12)))
       end do
       uss9(ie-ne+1,:,0 ) = uss9(ie-ne+1,:,12)
       uss9(ie-ne+1,:,13) = uss9(ie-ne+1,:,1 )
       do n=n1(ie), n2(ie)
         udes9(ie-ne+1,:,n) = udes9(ie-ne+1,:,n) - uss9(ie-ne+1,:,mod(n-n1(ie)+1,12))
       end do
     end if
   end if
 end do

 tss = ispan(0,13,1)*1.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,40,19)
 cnlev2 = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72*2
 he = 20
 xf = 10 + (/0,0,0,0,0/)*78
 yf = 10 + (/1,2,3,0,4/)*32

 mstr  = new(ne,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)
 lbopt = new(ne,string)

 mstr = (/"N1 from 1952/CTL","N2","N3","CTL","CTL99"/)
 xstr(:) = ""
 xstr(ne-1) = "month"
 ystr(:) = "height (km)"
 lbopt(:) = ""
 lbopt(ne-1) = "(m s:S:-1:N:)"
 xstr2 = xstr
 xstr2(ne-1) = "year"

do ie=0, ne-1

 if ( ie .lt. ne-1 ) then
   contour1 = cnshadeplot(wks ,tss,z ,uss  (ie,:,:),cnlev ,cnfil,lbopt(ie))
   contour2 = cnshadeplot(wks2,t  ,z ,udes (ie,:,:),cnlev2,cnfil,lbopt(ie))
 else
   contour1 = cnshadeplot(wks ,tss,z9,uss9 (0 ,:,:),cnlev ,cnfil,lbopt(ie))
   contour2 = cnshadeplot(wks2,t  ,z9,udes9(0 ,:,:),cnlev2,cnfil,lbopt(ie))
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
 setvalues contour2
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

 set_size(contour1,xf(ie),yf(ie),wi/4.,he)
 axis_range(contour1,0.5,12.5,35,65,False,False)
 axis_str(contour1,mstr(ie),xstr(ie),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBTickStartF"    : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 set_size(contour2,xf(ie),yf(ie),wi,he)
 axis_range(contour2,1942,1966,35,65,False,False)
 axis_str(contour2,mstr(ie),xstr2(ie),ystr(ie))
 setvalues contour2
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 1942
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks )
 cnfinalize(contour2,wks2)

 draw(contour1)
 draw(contour2)

end do


 frame(wks )
 frame(wks2)

end



