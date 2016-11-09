load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
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

 y1 = 1953
 y2 = 2015
 lat0 = 20

 z = (/21.0,23.3,26.9/)  ; (/10,20,30,50/) ~ (/30.6,26.9,23.3,21.0/)
 p = (/50,30,20/)

 wks = cr_wks("fig/u_t0-z-"+expname(0)+"_"+lat0,600,0)
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

 nz = dimsizes(z)

 u = new((/ne,nz,nt/),float,1.e20)
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
 lcol = (/2,1,3/)

;--------------------------------------- contour
 nl = 3
 ys = (/1953,1971,1989/)
 yl = 18
 yls = 2

 wi = 126
 he = 20
 xf = 34
 yf = 10 + (/0,1,2/)*28

 wi2 = tofloat(wi)/yl*yls

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = (/cexp(0),"",""/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "U (m s:S:-1:N:)"

do il=0, nl-1

 plot = xylineplot(wks,t,-6.,"D",0,"T",0.5,"C",1)
 do k=0, nz-1
   xyaddline(plot,t,u(0,k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,ys(il),ys(il)+yl,-30,15,False,False)
 axis_str(plot,mstr(il),xstr(il),ystr(il))
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

if (il .gt. 0) then
  setvalues plot
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  plot2 = xylineplot(wks,t,-6.,"D",0,"T",0.5,"C",1)
  do k=0, nz-1
    xyaddline(plot2,t,u(0,k,:),"D",0,"T",1.0,"C",lcol(k))
  end do
  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),-30,15,False,False)
  axis_str(plot2,"","",ystr(il))
  setvalues plot2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  draw(plot2)
end if

 draw(plot)

end do


 frame(wks)

end



