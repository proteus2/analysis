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
 y2 = 2006
 lat0 = 5.

 z = 26.9  ; (/10,20,30,50/) ~ (/30.6,26.9,23.3,21.0/)
 p = 50

 wks = cr_wks("fig/fu_t0-z"+toint(z)+"-"+expname(0),600,0)
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
     u = new((/ne,nt/),float,1.e20)
   end if
   u(ie,:) = dim_avg( f1->$varname$(:,{z*1.e3},{-lat0:lat0}) )
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("../dat/tem/uan"+expname(0)+".tem2_yzt.nc","r")
 utend = dim_avg( f2->u_tend(:,{z*1.e3},{-lat0:lat0}) )
 uforce = dim_avg( f2->u_force(:,{z*1.e3},{-lat0:lat0}) )
 epd    = dim_avg( f2->epd    (:,{z*1.e3},{-lat0:lat0}) )

 f3 = addfile("../dat/du_xm.uan"+expname(0)+"a.pe.nc","r")
 du_gw = dim_avg( f3->du_gw(:,{z*1.e3},{-lat0:lat0}) )*72.
 if (cexp(0) .ne. "CTL") then
   du_gwc = dim_avg( f3->dudt_gwconv(:,{z*1.e3},{-lat0:lat0}) )*86400.
   du_gw = du_gw + du_gwc
 end if

;--------------------------------------- process 0
 var = (/utend,epd   ,du_gw,uforce+du_gw/) ; ,u(0,:)/100./)

 nv = 4

;--------------------------------------- cn level
 lcol = (/1,4,2,3,6/)

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

 mstr = (/"","",""/);(/cexp(0),"",""/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "force (m s:S:-1:N: day:S:-1:N:)"

do il=0, nl-1

 plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
 do k=0, nv-1
   xyaddline(plot,t,var(k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,ys(il),ys(il)+yl,-0.6,0.6,False,False)
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

  plot2 = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
  do k=0, nv-1
    xyaddline(plot2,t,var(k,:),"D",0,"T",1.0,"C",lcol(k))
  end do
  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),-0.3,0.3,False,False)
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



