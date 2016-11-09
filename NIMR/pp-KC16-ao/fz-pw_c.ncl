load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"tc"/)
 expdir  = (/"AOL60CGW"/)
 varname = "unspecified"

 lat0 = 10.
 z0   = 17.41  ; 18.59  ; 17.41
 y1 = 1955 ;2058 ;1960
 y2 = 1994 ;2097 ;1999
 y0 = 1951
 ydif = 100
 nc = 61
 dc = 2.
 ang = 45.

 wks = cr_wks("fig/fz-pw_c-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
   fnam1 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.ph.nc"
;   fnam2 = "$DATD/"+expdir(0)+"/mmxm.uan"+expname(0)+"a.pg.nc"
   fnam3 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"
print(fnam1+"")
;print(fnam2+"")

   f1 = addfile(fnam1,"r")
;   f2 = addfile(fnam2,"r")
   f3 = addfile(fnam3,"r")
print(f1->$f1->unspecified!1$({z0*1.e3-200})+"  "+(z0*1.e3)+"  "+  \
      f1->$f1->unspecified!1$({z0*1.e3+200}))

   nt = (y2 - y1 + 1)*12
   n1 = (y1 - y0)*12
   n2 = n1 + nt - 1
   n1f = (y1+ydif - y0)*12
   n2f = n1f + nt - 1

   nw = 3

   rw = new((/nw,2,nc/),float)
do ic=0, nc-1
   if (ic .eq. 0) then
     varnamei0 = varname
   else
     varnamei0 = varname+"_"+ic
   end if
   varnamei1 = varname+"_"+(ic+nc)
   rw(0,0,ic) = dim_avg_n( 0.5*  \
          dim_avg( f1->$varnamei0$(n1:n2,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei0$(n1:n2,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(0,1,ic) = dim_avg_n( (-0.5)*  \
          dim_avg( f1->$varnamei1$(n1:n2,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei1$(n1:n2,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(1,0,ic) = dim_avg_n( 0.5*  \
          dim_avg( f1->$varnamei0$(n1f:n2f,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei0$(n1f:n2f,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
   rw(1,1,ic) = dim_avg_n( (-0.5)*  \
          dim_avg( f1->$varnamei1$(n1f:n2f,{z0*1.e3-200},{-lat0:lat0}) +  \
                   f1->$varnamei1$(n1f:n2f,{z0*1.e3+200},{-lat0:lat0}) ), 0)/dc
end do
   rw = abs(rw)

;   rw(2,:,:) = 0.5*(rw(0,:,:)+rw(1,:,:))
   rw(2,:,:) = rw(0,:,:)

   rw@_FillValue = -999.
   rw(2,:,:) = where(rw(2,:,:) .eq. 0., rw@_FillValue, rw(2,:,:))

;   rw(2,0,:) = rw(0,0,:) + rw(1,0,:)
;   rw(2,1,:) = rw(0,1,:) + rw(1,1,:)

   um = dim_avg( f3->u(n1:n2,{4.2e3:z0*1.e3},{-lat0:lat0}) )
   ud = new((/nw,2,nt/),float)
   do n=0, nt-1
     k1 = maxind( um(n,:) )
     k0 = minind( um(n,:) )
;     ud(0,0,n) = um(n,k1) - min( um(n,:k1) )
;     ud(0,1,n) = um(n,k0) - max( um(n,:k0) )
     ud(1,0,n) = max( um(n,1:) ) ;- um(n,0)
     ud(1,1,n) = min( um(n,1:) ) ;- um(n,0)
   end do
   ud(0,:,:) = ud(1,:,:)

   t = ispan(y1,y2,1)

   wss = new((/nw,2,nc/),float)
   wss = rw*1.e3

   uss = dim_avg( ud )

;--------------------------------------- process 0
 rwrdif = (wss(1,:,:) - wss(0,:,:))/wss(2,:,:)
; rwrdif = where(wss(0,:,:) .lt. 0.5e-4, rwrdif@_FillValue, rwrdif)

; wss = runave(wss,5,0)
; uss = runave(uss,5,0)

 c = ispan(-nc/2,nc/2,1)*dc*cos(ang*3.141592/180.)

 wss = where(wss .eq. 0., 1.e-10, wss)

;--------------------------------------- cn level
; cnlev = (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3/)
; cnfil = (/0 ,10,9 ,8 ,7  ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,20 ,19 ,18/)
; cnlev = (/0.05,0.1,0.25,0.5,1,2.5,3.5,5,6,7,8,9,10,15,20/)*0.2
 cnlev = (/0.1,0.2,0.4,0.7,1,2,4,7,10,20,30,40/)*0.1
 cnfil = (/0,5,6,7,8,12,14,15,16,17,19,20,21/)

 ldsh = (/0,16,0/)
 ltck = (/1.2,1.2,2.0/)

;--------------------------------------- contour
 wname = (/ (/"E-CGW","E-BGW","E-GW"/), (/"W-CGW","W-BGW","W-GW"/) /)
 wi = 30
 he = 20
 he2 = 10
 xf = 12 + (/0,1,2,3/)*36
 yf = 50 + (/0,1,2,3/)*28 + (/0,0,0,0/)*12

 mstr = ""
 xstr = "phase speed [m s:S:-1:N:]"
 ystr = "flux [mPa / m s:S:-1:N:]"
 ystr2 = "change [%]"
; txtstr = wname
 txtstr = (/"CGW","BGW"/)
; txtcol = (/ (/0,0,0/), (/0,1,1/) /)

; wpl = (/ wss(0,0,:), wss(0,1,:), wss(1,0,:), wss(1,1,:) /)
; plot = xylineplot(wks,c,wpl,"D",(/0,0,14,14/),"T",1.2,"C",(/1,2,1,2/))
 wpl = (/ wss(0,0,:), wss(0,1,:) /)
 plot = xylineplot(wks,c,wpl,"D",(/0,14/),"T",1.2,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-40,40,1.e-4,1.e-1,False,True)
 axis_str(plot,mstr,xstr,ystr)

 draw(plot)

 dpl = (/ rwrdif(0,:), rwrdif(1,:) /)*100.
 plot2 = xylineplot(wks,c,dpl,"D",(/0,14/),"T",1.2,"C",1)
 xyaddline(plot2,c,0,"D",0,"T",0.6,"C",colind)

 set_size(plot2,xf,yf-he2-3,wi,he2)
 axis_range(plot2,-40,40,-30,60,False,False)
 axis_str(plot2,"","",ystr2)
 setvalues plot2
  "tmXBLabelsOn" : False
 end setvalues

 draw(plot2)

; dr_txt2(wks,txtstr(iv),xf+0.6,yf(iv)+he-2,22,5,0.,"CenterLeft",0)


 frame(wks)

end

