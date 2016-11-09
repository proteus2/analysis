load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 u_expname = (/"era-int","merra","cfsr","jra55"/);,"era-int_f"/)
 u_datname0= (/"ERA-I"  ,"MERRA","CFSR","JRA55"/);,"ERA-I-nr"/)
; ie = 0

 y0 = 1979
 feb_1st = True
 y1 = 1979
 y2 = 2010
 lat0u = 5
 lata = lat0u + 2.5
 p = 30
 p_em = 19  ; (50,30,20,10) : (22,19,17,14)

 wks = cr_wks("fig/u_t-p"+p,600,0)
; gsn_define_colormap(wks,"amwg")
; gsn_define_colormap(wks,"cosam")
; gsn_define_colormap(wks,"radar")
 colind = NhlNewColor(wks,0.70,0.70,0.70)

;--------------------------------------- read data
 ne = dimsizes(u_expname)
 nw = 8
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12
 ntr = nt
 if ( feb_1st .and. ( n0 .eq. 0 ) ) then
   ntr = nt - 1
 end if

 u_vname = (/"U_GDS0_ISBL","u","UGRD","U_GRD_GDS0_ISBL"/);,"u"/)
 u_fnam0 = "$DATD/clim/"+u_datname0+"/"+u_expname+".u_ypt.nc"
 uc = new(nt,float)
 uc = 0.
 do ie=0, dimsizes(u_expname)-1
   f0 = addfile(u_fnam0(ie),"r")
   s_lat = f0->$u_vname(ie)$!2
   lat = tofloat( f0->$s_lat$({-lata:lata}) )
   u0 = new(nt,float)
;   if (u_expname(ie) .eq. "era-int_f") then
;     do n=0, nt-1
;       u0(n) = ftcurvi(-lat0u,lat0u,lat, \
;               tofloat( f0->$u_vname(ie)$(n0+n,{p_em},{-lata:lata}) ))
;     end do
;   else
     do n=0, nt-1
       u0(n) = ftcurvi(-lat0u,lat0u,lat, \
               tofloat( f0->$u_vname(ie)$(n0+n,{p},{-lata:lata}) ))
     end do
;   end if
   u0 = u0 / (lat0u+lat0u)
   uc = uc + u0/tofloat(dimsizes(u_expname))
   delete(lat)
 end do

   t = y1 + (fspan(1,nt,nt)-0.5)/12.
;   tt = y1 + (fspan(1,nt-1,nt-1))/12.

;--------------------------------------- process 0
 zeropt = new(nt,integer)
 zeropt = -999

 tag = 0
 nlast = -999
 nlast2 = -999
 cnt = 0
 do n=1, nt-1
   if ( uc(n)*uc(n-1) .le. 0. ) then
     zeropt(cnt) = n
     if (n-nlast2 .lt. 6) then
       zeropt(cnt-1:cnt) = -999
       cnt = cnt - 1
     else
       cnt = cnt + 1
     end if
     nlast2 = nlast
     nlast = n
   end if
 end do
 nzero = cnt

 ;yh +1 for 30 hPa
 nzero = cnt + 1
 zeropt(nzero-1) = nt-1

 maxmin = new((/2,nzero-1/),integer)
 do cnt=0, nzero-2
   if ( sum(uc(zeropt(cnt):zeropt(cnt+1))) .gt. 0. ) then
     maxmin(0,cnt) = maxind(uc(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = 1
   else
     maxmin(0,cnt) = minind(uc(zeropt(cnt):zeropt(cnt+1))) + zeropt(cnt)
     maxmin(1,cnt) = -1
   end if
 end do

;--------------------------------------- cn level
 lcol = (/2,3,13,4,1,7/)  ; default
; lcol = (/15,9,13,3,1,17/)  ; amwg
; lcol = (/11,7,10,4,1,17/)  ; cosam
; lcol = (/11,7,9,4,1,17/)  ; radar
 ltck = (/1,1,1,1,1/)*0.8
; ltck = 1.2 + 0.2*fspan(0,-6,7)
 imrk = (/1,1,1,1,1/)

;--------------------------------------- contour
 wi = 140
 he = 15
 xf = 15
 yf = 10

 xstr = "year"
 ystr = "(m s:S:-1:N:)"
 txtstr = "U"

 plot = xylineplot(wks,t,uc,"D",0,"T",0.5,"C",1)
; xyaddmarker(plot,imrk,"S",0.8,"T",10.0,"C",lcol)

 xyaddline(plot,t,0.,"D",0,"T",0.6,"C",1)

 dumpb = new(nzero,graphic)
 respg = True
 respg@gsFillColor = colind
 do cnt=0, nzero-3
   if (maxmin(1,cnt) .eq. -1) then
     polybx = (/ t(maxmin(0,cnt)), t(maxmin(0,cnt+1)), t(maxmin(0,cnt+1)), \
                 t(maxmin(0,cnt)), t(maxmin(0,cnt)) /)
     polyby = (/-1,-1,1,1,-1/)*999
     dumpb(cnt) = gsn_add_polygon(wks,plot,polybx,polyby,respg)
   end if
 end do

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,y1,y2+1,-999,999,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 draw(plot)

 dr_txt2(wks,txtstr,xf+1,yf+he-2,22,6,0.,"CenterLeft",1)

 frame(wks)

end

