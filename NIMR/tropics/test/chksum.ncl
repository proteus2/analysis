load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "uj"

 y1 = 1953
 y2 = 2006
 lat0 = 5.

 za = (/36.0,30.6,26.9,23.3,21.0,18.6/)

 wks = cr_wks("fig/chksum-za.ps",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

do kkk=0, 5

; z = 26.9  ; (/10,20,30,50/) ~ (/30.6,26.9,23.3,21.0/)
 z = za(kkk)

; wks = cr_wks("fig/chksum-z"+numeric2int(z,1),600,0)
; colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = 1
 nt = (y2-y1+1)*12 - 1

 fname = "$DATD/L60CGW-t/xm.uan"+expname+"a.pj.nc"
 f1 = addfile(fname,"r")
 u = dim_avg( f1->u(:nt-1,{z*1.e3},{-lat0:lat0}) ) *1.e-2

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 f2 = addfile("$DATD/L60CGW-t/tem/uan"+expname+".tem_yzt.nc","r")
 epd = dim_avg( f2->epd(:nt-1,{z*1.e3},{-lat0:lat0}) )
 epd = (runave_n(epd,3,0,0)*3. + epd )*0.25
; epd = runave_n(epd,3,0,0)

 f3 = addfile("$DATD/L60CGW-t/epf/uan"+expname+".epf_yz_recon0.nc","r")
 epde = dim_avg( f3->epd(:,0,{z*1.e3},{-lat0:lat0}) )
 epdw = dim_avg( f3->epd(:,1,{z*1.e3},{-lat0:lat0}) )
 epde1 = dim_avg( f3->epd(:,4,{z*1.e3},{-lat0:lat0}) ) + \
         dim_avg( f3->epd(:,6,{z*1.e3},{-lat0:lat0}) )
 epdw1 = dim_avg( f3->epd(:,5,{z*1.e3},{-lat0:lat0}) ) + \
         dim_avg( f3->epd(:,7,{z*1.e3},{-lat0:lat0}) )
 epde2 = dim_avg( f3->epd(:,8,{z*1.e3},{-lat0:lat0}) )
 epdw2 = dim_avg( f3->epd(:,9,{z*1.e3},{-lat0:lat0}) )
 epdw3 = dim_avg( f3->epd(:,2,{z*1.e3},{-lat0:lat0}) ) + \
         dim_avg( f3->epd(:,3,{z*1.e3},{-lat0:lat0}) )

 epdea = epde1 + epde2 
 epdwa = epdw1 + epdw2

 f4 = addfile("$DATD/L60CGW-t/epf/uan"+expname+".epfsa0_yz_recon0.nc","r")
 epde1s = dim_avg( f4->epd(:,0,{z*1.e3},{-lat0:lat0}) )
 epdw1s = dim_avg( f4->epd(:,1,{z*1.e3},{-lat0:lat0}) )
 epde1s1 = dim_avg( f4->epd(:,2,{z*1.e3},{-lat0:lat0}) )
 epdw1s1 = dim_avg( f4->epd(:,3,{z*1.e3},{-lat0:lat0}) )
 epde1s2 = dim_avg( f4->epd(:,4,{z*1.e3},{-lat0:lat0}) )
 epdw1s2 = dim_avg( f4->epd(:,5,{z*1.e3},{-lat0:lat0}) )
 epde1s3 = dim_avg( f4->epd(:,6,{z*1.e3},{-lat0:lat0}) )

 f5 = addfile("$DATD/L60CGW-t/epf/uan"+expname+".epfsa1_yz_recon0.nc","r")
 epde2s = dim_avg( f5->epd(:,0,{z*1.e3},{-lat0:lat0}) )
 epdw2s = dim_avg( f5->epd(:,1,{z*1.e3},{-lat0:lat0}) )

 epde1a = epde1s1 + epde1s2
 epdw1a = epdw1s1 + epdw1s2

 epde9 = epde1s + epde2
 epdw9 = epdw1s + epdw2 + epdw3

 epde8 = epde9 + epde2s
 epdw8 = epdw9 + epdw2s

;--------------------------------------- process 0
; var = (/u, epd, epdw, epdea+epdwa+epdw3-epd, epdwa+epdw3-epdw /)
; var = (/u, epde1, epdw1, epde1s, epdw1s /)
; var = (/u, epde1s, epdw1s, epde1a, epdw1a /)
; var = (/u, epde, epdw, epde9, epdw9 /)
 var = (/u, epd, epde+epdw-epd, epde9+epdw9, epde9+epdw9-epd /)

 nv = 5

;--------------------------------------- cn level
 lcol = (/1,4,2,3,6,5,8/)

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

 mstr = (/"","",""/);(/cexp,"",""/)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "force (m s:S:-1:N: day:S:-1:N:)"

do il=0, nl-1

 plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
 do k=0, nv-1
   xyaddline(plot,t,var(k,:),"D",0,"T",1.0,"C",lcol(k))
 end do

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,ys(il),ys(il)+yl,-0.3,0.3,False,False)
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

end do  ; kkk

end



