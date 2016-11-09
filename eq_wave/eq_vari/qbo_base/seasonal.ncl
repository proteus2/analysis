load "$UTIL/ncl/header.ncl"

xspace = 160  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 load "./ftn_seasonal.ncl"
 load "./param0.ncl"

 fname = param0_fname
 vname = param0_vname

 t1 = param0_t1
 t2 = param0_t2

 t1_1 = param0_t1_1
 t2_1 = param0_t2_1
 t1_2 = param0_t1_2
 t2_2 = param0_t2_2
 t1_3 = param0_t1_3
 t2_3 = param0_t2_3

 if ( tofloat(toint(t1)) .ne. t1 ) then
   print("t1 should be integer.")
   exit
 end if

 ; for the title
 t2t = t2
 if ( tofloat(toint(t2)) .eq. t2 ) then
   l_todec = True
   t2t = t2t - 1
 end if

 wks = cr_wks("fig/seasonal_"+t1+"-"+t2t,600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 f1 = addfile(fname,"r")
 u0 = f1->$vname$
 tmpstr = f1->$vname$!0
 p0 = f1->$tmpstr$
 tmpstr = f1->$vname$!1
 t0 = f1->$tmpstr$

 u0@_FillValue = -999.

 p = p0
 nz = dimsizes(p)

;--------------------------------------- process 0
 u   = u0
 u_1 = u0
 u_2 = u0
 u_3 = u0
 do k=0, nz-1
   u  (k,:) = where( t0 .lt. t1   .or. t0 .gt. t2   , u0@_FillValue, u  (k,:) )
   u_1(k,:) = where( t0 .lt. t1_1 .or. t0 .gt. t2_1 , u0@_FillValue, u_1(k,:) )
   u_2(k,:) = where( t0 .lt. t1_2 .or. t0 .gt. t2_2 , u0@_FillValue, u_2(k,:) )
   u_3(k,:) = where( t0 .lt. t1_3 .or. t0 .gt. t2_3 , u0@_FillValue, u_3(k,:) )
 end do

 uss = new((/4,3,nz,12/),float)
 do ih=0, 2
   uss(0,ih,:,:)  = seasonal(u  ,ih)
   uss(1,ih,:,:)  = seasonal(u_1,ih)
   uss(2,ih,:,:)  = seasonal(u_2,ih)
   uss(3,ih,:,:)  = seasonal(u_3,ih)
 end do

 t = ispan(0,13,1)
 ua = new((/4,3,nz,14/),float)
 ua(:,:,:,1:12) = uss
 ua(:,:,:,0) = ua(:,:,:,12)
 ua(:,:,:,13) = ua(:,:,:,1)

;--------------------------------------- cn level
 lcol = (/1,2,4,3/)
 ltck = (/1.5,0.7,0.7,0.7/)

;--------------------------------------- contour
 nl = nz

 wi = 30
 he = 20
 xf = 10 + (/0,1,2,3/)*(wi+10) - (/0,0,3,6/)
 yf = 10 + ispan(nl-1,0,1)*(he+5)

 mstr0 = (/"All stations","Seasonal cycle","1st harmonic","Harmonics 1+2"/)
 mstr = new(nl,string)
 mstr = ""
 xstr = new(nl,string)
 xstr = ""
 xstr(0) = "month"
 ystr = "U [m s:S:-1:N:]"

do il=0, nl-1

 mstr(nl-1) = mstr0(0)
 plot = xylineplot(wks,t,ua(0,:,il,:),"D",0,"T",ltck,"C",lcol)
 xyaddmarker(plot,16,"S",0.8,"T",1.0,"C",1)

 set_size(plot,xf(0),yf(il),wi,he)
 axis_range(plot,0.5,12.5,-999,999,False,False)
 axis_str(plot,mstr(il),xstr(il),ystr)
 setvalues plot
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : ispan(1,12,1)
 end setvalues

 if ( p(il) .gt. 20 .and. p(il) .lt. 50 ) then
   dr_txt(wks,p(il)+" hPa",xf(0)+5,yf(il)+he-2,22,5,0.)
 else
   dr_txt(wks,p(il)+" hPa",xf(0)+5,yf(il)+2,22,5,0.)
 end if

 draw(plot)

 lg = lgplot(wks,(/"total","1st harmonic","harmonics 1+2"/),"D",0,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(0),yf(0)+he+8,wi*0.5,12)
 draw(lg)


do ih=0, 2

 mstr(nl-1) = mstr0(ih+1)
 plot = xylineplot(wks,t,ua(:,ih,il,:),"D",0,"T",ltck,"C",lcol)
 xyaddmarker(plot,16,"S",0.8,"T",1.0,"C",1)

 set_size(plot,xf(ih+1),yf(il),wi,he)
 axis_range(plot,0.5,12.5,-999,999,False,False)
 axis_str(plot,mstr(il),xstr(il),"")
 setvalues plot
  "tmXBMode" : "Explicit"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : ispan(1,12,1)
 end setvalues

 draw(plot)

end do  ; ih

 s_stations = (/"all stations","Canton Island","Gan/Maledive Islands","Singapore"/)
 if ( l_todec ) then
   s_tt = (/"Jan. ","Dec. "/)+(/t1,t2t/)
 else
   s_tt = (/"",""/)
 end if
 s_periods = (/""," ("+s_tt(0)+" ~ Aug. 1967)"," (Sep. 1967 ~ Dec. 1975)"," (Jan. 1976 ~ "+s_tt(1)+")"/)
 lg = lgplot(wks,s_stations+s_periods,"D",0,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(1),yf(0)+he+8,wi*1.0,16)
 draw(lg)

end do  ; il

 frame(wks)

end

