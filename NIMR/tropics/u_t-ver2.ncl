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
 p = (/50,30,20/)

 wks = cr_wks("fig/u_t-ver2",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nz = dimsizes(p)

 fdir  = "/data15/kyh/qbo/obs"
 fname1 = fdir+"/qbo_"+expname(0)+".nc"
 f1 = addfiles(fname1,"r")

 u1 = f1[0]->u({p},:)

 t1 = f1[0]->t


 fdir  = "../dat"
 fname3 = fdir+"/"+expname(1:2)+".u_ypt.1979-2010.nc"
 f3 = addfiles(fname3,"r")

 u3 = transpose( doubletofloat(dim_avg( f3[0]->u(:,{p},{-lat0:lat0}) )) )
 u4 = transpose( dim_avg( f3[1]->u(:,{p},{-lat0:lat0}) ) )

 nt3 = (2010-1979+1)*12
 t3 = 1979. + (fspan(1,nt3,nt3)-0.5)/12.
 t4 = t3

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
 wi = 72*2
 he = 20
 xf = 10
 yf = 10 + (/0,9,9,1,3,5,2,4,6/)*28 + (/0,1,2,0,1,2,0,1,2/)*10

 mstr  = new(ne*3,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)

 mstr(:) = ""
 mstr(0) = "Station data"
 mstr(4) = "ERA-Interim"
 mstr(5) = "MERRA"
 xstr(:) = ""
 xstr(ne-1) = "year"
 ystr(:) = "U (m s:S:-1:N:)"

 plot = new(3,graphic)
 y1 = (/1953,1973,1993/)
do ie=0, ne-1

do i=0, 2
 if ( ie .eq. 0 ) then
   plot(i) = xylineplot(wks,t1,0.,"D",0,"T",0.5,"C",1)
   do k=0, nz-1
     xyaddline(plot(i),t1,u1(k,:),"D",0,"T",1.0,"C",lcol(k))
   end do
 end if
 if ( ie .eq. 1 ) then
   plot(i) = xylineplot(wks,t3,0.,"D",0,"T",0.5,"C",1)
   do k=0, nz-1
     xyaddline(plot(i),t3,u3(k,:),"D",0,"T",1.0,"C",lcol(k))
   end do
 end if
 if ( ie .eq. 2 ) then
   plot(i) = xylineplot(wks,t4,0.,"D",0,"T",0.5,"C",1)
   do k=0, nz-1
     xyaddline(plot(i),t4,u4(k,:),"D",0,"T",1.0,"C",lcol(k))
   end do
 end if

 set_size(plot(i),xf,yf(ie+i*ne),wi,he)
 axis_range(plot(i),y1(i),y1(i)+22,-40,30,False,False)
 axis_str(plot(i),mstr(ie+i*ne),"",ystr(ie))
 setvalues plot(i)
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : y1(i)
  "tmXBMinorPerMajor" : 0
 end setvalues

 if ( ie .eq. 0 .or. i .ge. 1 ) then
   draw(plot(i))
 end if
end do

end do


 frame(wks)

end

