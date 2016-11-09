load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
; expname = (/"te","tf","tg","th","ui","up"/)
; ppid    = (/"e" ,"e" ,"e" ,"e" ,"j" ,"j" /)
; varname = "u_"+(/1,1,1,1,1,1/)
 expname = "ui"
 ppid    = "j"
 varname = "u_1"

 y1 = 1942
 y2 = 1969

 wks = cr_wks("fig/spin_u_yz-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*12

 fname = new(nt,string)
 fname2= fname
 fdir  = "../dat/xm"
 n = 0
 do iy=y1, y2, 1
 do im=1, 12
   umtail5 = systemfunc("ex_date_umstr.sh "+iy+" "+im)
   fname(n) = fdir+"/"+expname+"/u_xm.uan"+expname+"a.p"+ppid+umtail5+".nc"
   fname2(n) = fdir+"/"+expname+"/pt_xm.uan"+expname+"a.p"+ppid+umtail5+".nc"
   n = n + 1
 end do
 end do
 fexist = isfilepresent(fname)

 n1 = -1
 n2 = nt-1
 do n=0, nt-1
   if ( fexist(n) ) then
     if (n1 .eq. -1) then
       n1 = n
     end if
   else
     if ( n1 .ne. -1 .and. n2 .eq. nt-1 ) then
       n2 = n - 1
     end if
   end if
 end do
 if (n1 .eq. -1) then
   n2 = -1
 end if
print(n1+"  "+n2)

 ; get z
 f0 = addfile(fname(n1),"r")
 tmp = f0->$varname$
 zname = tmp!1
 z = tmp&$zname$/1.e3
 nz = dimsizes(z)
 delete(tmp)

 lat = f0->latitude
 ny = dimsizes(lat)

 u = new((/2,nt/24,nz,ny/),float,1.e20)
print(expname+"")
 f1 = addfiles(fname(n1:n2),"r")
 ListSetType(f1,"cat")
 f2 = addfiles(fname2(n1:n2),"r")
 ListSetType(f2,"cat")
 tmp = new((/2,nt/24,nz,ny/),float,1.e20)
 do n=0, nt/24-1
   tmp(0,n,:,:) = dim_avg_n( f1[n*24:n*24+23]->$varname$(:,:,:), 0)
   tmp(1,n,:,:) = dim_avg_n( f2[n*24:n*24+23]->theta(:,:,:), 0)
 end do
 do k=0, nz-1
   tmp(1,:,k,:) = tmp(1,:,k,:)*(exp(-z(k)/7.))^0.286 - 273.
 end do
 u(:,:,:,:) = tmp(:,:,:,:)
 do n=0, nt/24-1
   u(:,n,:,:) = tmp(:,n,:,:) - tmp(:,4,:,:)
 end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)/2.
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 30
 he = 30
 xf = 10 + (/0,1,2,3,4,5,6,0,1,2,3,4,5,6/)*32
 yf = 10 + (/0,0,0,0,0,0,0,1,1,1,1,1,1,1/)*90

 mstr = expname
 xstr = "year"
 lbopt = "(m s:S:-1:N:)"

do i=0, 1
do n=0, nt/24-1

 contour1 = cnshadeplot(wks,lat,z,u(i,n,:,:),cnlev/(i+1.),cnfil,lbopt)
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

 set_size(contour1,xf(n),yf(n)+i*40,wi,he)
 axis_range(contour1,-999,999,0,999,False,False)
 axis_str(contour1,mstr,xstr,"")

 cnfinalize(contour1,wks)


 draw(contour1)

end do
end do

 frame(wks)

end



