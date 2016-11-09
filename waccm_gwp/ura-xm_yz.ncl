load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

do loopctl=1, 4

begin

;--------------------------------------- parameters
 hdname  = (/"test1_0.9","test1_0.9"/)
 expname = (/"2"        ,"2"    /)
 varname = (/"U"/)
 titname = (/""/)

 y1 = 1
 y2 = 5
; monselect = 1  ; 0,1,2,3,4,5 - Annual, DJF, MAM, JJA, SON, AMJ
 monselect = loopctl

 zmin = 0
 zmax = 65

;--------------------------------------- workstation
 if (monselect .eq. 0) then
   mon = ispan(1,12,1)
   tstr = "Annual"
 else
   monlist = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/), (/4,5,6/) /)
   tstrlist = (/"DJF","MAM","JJA","SON","AMJ"/)
   mon = monlist(monselect-1,:)
   tstr = tstrlist(monselect-1)
 end if

 wks = cr_wks("fig/ura-xm_yz-"+tstr,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
 nt = (y2-y1+1)*dimsizes(mon)
 nv = dimsizes(varname)

 fname = new((/nt,2/),string)
 n = 0
 do iy=y1, y2
 do im=0, dimsizes(mon)-1
   fname(n,:) = "dat"+expname+"/"+hdname+".cam2.h0."+ \
                sprinti("%0.4i",iy)+"-"+sprinti("%0.2i",mon(im))+".nc"
   n = n + 1
 end do
 end do

 fname4 = "dat/merra.u_ypt."+(/"djf","mam","jja","son"/)+".1979-2012.nc"
 f4 = addfile(fname4(monselect-1),"r")
 u_ra = f4->u(0,:,:)
 y_ra = f4->latitude
 z_ra = 7.0*log(1.e3 / f4->levels)

 ; get axis
 f0 = addfile(fname(0,0),"r")
 p = f0->lev
 y = f0->lat
 nz = dimsizes(p)
 ny = dimsizes(y)
 z = 7.0*log(1.e3/p)

 ; get variables
 var0 = new((/nv,nz,ny/),float,1.e20)
 var0 = 0.
 var = new((/nv,nz,ny/),float,1.e20)
 var = 0.
 do n=0, nt-1
 do iv=0, nv-1
   f1 = addfile(fname(n,0),"r")
   f2 = addfile(fname(n,1),"r")
   tmp0 = dim_avg_n( dim_avg( f2->$varname(iv)$(:,:,:,:) ), 0 )
   tmp  = dim_avg_n( dim_avg( f1->$varname(iv)$(:,:,:,:) ), 0 ) - tmp0
   var0(iv,:,:) = var0(iv,:,:) + tmp0(:,:)
   var (iv,:,:) = var (iv,:,:) + tmp (:,:)
 end do
 end do
 var0(:,:,:) = var0(:,:,:)/tofloat(nt)
 var (:,:,:) = var (:,:,:)/tofloat(nt)

; t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0
 z_r = z
 var_r = var
 var0_r = var0
 do k=0, nz-1
   z_r(k) = z(nz-1-k)
   var_r(:,k,:) = var(:,nz-1-k,:)
   var0_r(:,k,:) = var0(:,nz-1-k,:)
 end do

 var_ra = linint2(y,z_r,var0_r,False,y_ra,z_ra,0)
 var_ra = var_ra(0,:,:) - u_ra

;--------------------------------------- cn level
; cnlev = fspan(-30,30,31)
 cnlev = fspan(-150,150,31)
; cnlev = (/-1000,-500,-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200,500,1000/)*0.1
; cnfil = (/17,19,21,22,23,24,25,26,0,0,6,5,4,3,2,11,9,7/)

;--------------------------------------- contour
 wi = 30
 he = 30*(65./120.)
 xf = 15 + (/0,1,2/)*(wi+8)
 yf = 10 + (/0,0,0/)*(he+8)

 xstr  = new(nv,string,"")
 ystr  = new(nv,string,"")

 xstr(:) = "latitude (deg)"
 ystr(0) = "height (km)"

do iv=0, nv-1

 contour0 = cnshadeplot(wks,y,z,var0(iv,:,:),(/0./),(/colind,0/),"")

 set_size(contour0,xf(iv),yf(iv),wi,he)
 axis_range(contour0,-999,999,zmin,zmax,False,False)
 axis_str(contour0,titname(iv),xstr(iv),ystr(iv))

 cnfinalize(contour0,wks)

; draw(contour0)

 contour1 = cnlineplot(wks,y_ra,z_ra,u_ra,cnlev,0)
 cnline_att(contour1,"D",(/0,0,16/),"T",(/1.,0.4,1./),"C",(/1,colind,1/))

 set_size(contour1,xf(iv),yf(iv),wi,he)
 axis_range(contour1,-999,999,zmin,zmax,False,False)
 axis_str(contour1,titname(iv),xstr(iv),ystr(iv))

 cnfinalize(contour1,wks)

 draw(contour1)

end do


 frame(wks)

end

end do  ; nnn

