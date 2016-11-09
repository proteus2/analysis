load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 jobname = (/"tc","tk"/)
 expname = (/"CC","CTL"/)
 dirname = "AOL60"+(/"CGW","CTL"/)
 varname = (/"w_res","rho"/)

 y0 = 1951
 y1 = 1951
 y2 = 2097
 p = 70.
 z_rng = (/16,21/)*1.e3
 lat0 = 50.
 navg = 30
 mms = "jja"

 wks = cr_wks("fig/w70_t-"+mms+"-"+jobname(dimsizes(jobname)-1),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 ne = dimsizes(jobname)
 nyr = y2 - y1 + 1

 mm3s = (/"djf","mam","jja","son"/)
 mm3i = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
 do is=0, dimsizes(mm3s)-1
   if ( mms .eq. mm3s(is) ) then
     mm = mm3i(is,:)
   end if
 end do
 nm = dimsizes(mm)

do ie=0, ne-1

 fname = "../dat/"+dirname(ie)+"/tem/uan"+jobname(ie)+".tem-wc_yzt.nc"
 fnamer = "../dat/"+dirname(ie)+"/xm.uan"+jobname(ie)+"a.pj.nc"
 f0 = addfile(fname,"r")
 fr = addfile(fnamer,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp = f0->$varname(0)$(n0:n0+nyr*12-1,{z_rng(0):z_rng(1)},{-lat0:lat0})
 if ( .not. isdefined("lat") ) then
  lata = tmp&$tmp!2$
  lata!0 = "lat"
  lata&lat = lata
  lat = lata({-lat0:lat0})
  ny = dimsizes(lat)
  nz = dimsizes(tmp(0,:,0))
  w = new((/ne,nyr,nm,nz,ny/),float)
  lnp = new((/ne,nyr,nm,nz,ny/),float)
 end if
 rho = fr->$varname(1)$(n0:n0+nyr*12-1,{z_rng(0):z_rng(1)},{-lat0:lat0})
 pa = fr->p(n0:n0+nyr*12-1,{z_rng(0):z_rng(1)},{-lat0:lat0})

 do im=0, nm-1
   w(ie,:,im,:,:) = tmp(ispan(mm(im)-1,nyr*12-1,12),:,{-lat0:lat0})* \
                    rho(ispan(mm(im)-1,nyr*12-1,12),:,{-lat0:lat0})
   lnp(ie,:,im,:,:) = log(pa(ispan(mm(im)-1,nyr*12-1,12),:,{-lat0:lat0}))
 end do

 delete(tmp)
 delete(rho)

end do  ; ie

;--------------------------------------- process 0
 do j=0, ny-1
   w(:,:,:,:,j) = w(:,:,:,:,j)*cos(lat(j)*3.141592/180.)
 end do
 dlat = lat(2) - lat(1)
 wgt = (2.*3.141592*6370.e3)^2*(dlat/360.) * 1.e-6
 w = w*wgt

 tmp = exp( linint1_n(lnp*(-1.),log(w),False,log(p*100.)*(-1.),0,3) )
 if ( dimsizes(dimsizes(tmp)) .eq. 5 ) then
   w_p = tmp(:,:,:,0,:)
 else
   w_p = tmp
 end if

 w_p = where(w_p .gt. 0, w_p, 0)
 ws = dim_avg(dim_sum(w_p))

 ws10 = runave(ws,navg,0)

 t = ispan(y1,y2,1)
 if (navg/2*2 .eq. navg) then
   t10 = ispan(y1,y2-1,1)+0.5
 else
   t10 = t
 end if

;--------------------------------------- cn level
 scl = 3
 ldsh = (/0,14/)

;--------------------------------------- contour
 wi = 60
 he = 40
 xf = 20
 yf = 10

 xstr = "year"
 ystr = "mass flux (10:S:"+scl+":N: kt s:S:-1:N:)"

 plot = xylineplot(wks,t,ws/(10.^scl),"D",ldsh,"T",1.0,"C",colind)
 xyaddline(plot,t10,ws10/(10.^scl),"D",ldsh,"T",1.0,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-999,999,-999,999,False,False)
 axis_str(plot,"",xstr,ystr)

 draw(plot)

 if (ne .gt. 1) then
  lg = lgplot(wks,expname,"D",ldsh,"T",1.0,"C",1,2)
  set_size_lg(lg,xf+1,yf,15,10)
  draw(lg)
 end if

 frame(wks)

end

