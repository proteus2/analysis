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

 wks = cr_wks("fig/u_tp-ver",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)

 fdir  = "/data15/kyh/qbo/obs"
 fname1 = fdir+"/qbo_"+expname(0)+".nc"
 f1 = addfiles(fname1,"r")

 u1 = f1[0]->u

 t1 = f1[0]->t

 z1 = 7.0*log(1.e3/f1[0]->p)


 fdir  = "../dat"
 fname3 = fdir+"/"+expname(1:2)+".u_ypt.1979-2010.nc"
 f3 = addfiles(fname3,"r")

 z3 = 7.0*log(1.e3/f3[0]->levelist)
 z4 = 7.0*log(1.e3/f3[1]->levels)

 u3 = transpose( doubletofloat(dim_avg( f3[0]->u(:,:,{-lat0:lat0}) )) )
 u4 = transpose( dim_avg( f3[1]->u(:,:,{-lat0:lat0}) ) )

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
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72*2
 he = 20
 xf = 10
 yf = 10 + (/0,9,9,1,2,3,4,5,6/)*28 + (/0,0,0,1,1,1,2,2,2/)*10

 mstr  = new(ne*3,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)

 mstr(0:2) = (/"Station data","",""/)
 mstr(3:5) = mstr(0:2)+(/", ERA-Interim, MERRA","",""/)
 mstr(6:8) = ""
 xstr(:) = ""
 xstr(ne-1) = "year"
 ystr(:) = "height (km)"
 lbopt = (/"","","(m s:S:-1:N:)"/)

 contour = new(3,graphic)
 y1 = (/1953,1973,1993/)
do ie=0, ne-1

do i=0, 2
 if ( ie .eq. 0 ) then
   contour(i) = cnshadeplot(wks,t1,z1,u1,cnlev,cnfil,0)
 end if
 if ( ie .eq. 1 ) then
   contour(i) = cnshadeplot(wks,t3,z3,u3,cnlev,cnfil,0)
 end if
 if ( ie .eq. 2 ) then
   contour(i) = cnshadeplot(wks,t4,z4,u4,cnlev,cnfil,lbopt(i))
 end if
 setvalues contour(i)
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

 set_size(contour(i),xf,yf(ie+i*ne),wi,he)
 axis_range(contour(i),y1(i),y1(i)+22,16,45,False,False)
 axis_str(contour(i),mstr(ie+i*ne),"",ystr(ie))
 setvalues contour(i)
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : y1(i)
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour(i),wks)

 if ( ie .eq. 0 .or. i .ge. 1 ) then
   draw(contour(i))
 end if
end do

end do


 frame(wks)

end

