load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 casecode = "ra"
; casecode = "hg"

 if (casecode .eq. "ra") then
   load "./param-ra-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ; (/"CGW","CTL"/)
 expname = param0_expname  ;(/"tc" ,"tk" /)
 fdir_u0 = param0_fdir_u0  ; "$DATD/AOL60"+cexp
 fname_u0 = param0_fname_u0  ; "$DATD/AOL60"+cexp
 u0varname = param0_u0varname  ;"u_p_uv"

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 nma = param0_nma
 nph = param0_nph

 xmin = param0_umin_p1
 xmax = param0_umax_p1
 ymin = param0_umin_p2
 ymax = param0_umax_p2

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = param0_fname_uv
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 vnamef = "circ"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_yz-"+nma+"m"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
 else
   if (expname(0) .eq. "era-int_f") then
     foutname = vnamef+"_yz-"+nma+"m"+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_yz-"+nma+"m"+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 u = new((/ne,nt/),float,1.e20)
 ush = new((/ne,nt/),float,1.e20)
 do ie=0, ne-1
   f1 = addfile(fname_u0(ie),"r")
   u(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p1},{-lat0:lat0}) ) )
   ush(ie,:) = tofloat( dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) )
   if (dimsizes(p2) .eq. 2) then
     ush(ie,:) = dim_avg( f1->$u0varname$(it0:it9,{p2(0)},{-lat0:lat0}) ) -  \
                 dim_avg( f1->$u0varname$(it0:it9,{p2(1)},{-lat0:lat0}) )
   end if
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 if (casecode .eq. "ra") then
   fname = "$DATD/tem/"+cexp_v(0)+"/2000/"+expname_v(0)+".tem_ypt.2000.01.nc"
 end if
 if (casecode .eq. "hg") then
;   fname = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname(0)+"a.pj_"+y1+"01.nc"
 end if
 ff = addfile(fname,"r")
 lat = ff->lat
 p   = ff->p({500:10})
 z   = 7.0*log(1000./p)
 ny = dimsizes(lat)
 nz = dimsizes(p)

 var = new((/ne,nt,3,nz,ny/),float,1.e20)
 n = 0
 do yy=y1, y2
 do m=1, 12
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)
   if (casecode .eq. "ra") then
     fname = "$DATD/tem/"+cexp_v+"/"+yy+"/"+expname_v+".tem_ypt."+yy+"."+mm+".nc"
     fnameu = "$DATD/clim/"+cexp_v+"/monthly/"+yy+"/"+expname_v+".u_yp."+yyyymm+".nc"
     do ie=0, ne-1
       ff = addfile(fname(ie),"r")
       ffu = addfile(fnameu(ie),"r")
       var(ie,n,0,:,:) = ff->v_res({500:10},:)
       var(ie,n,1,:,:) = ff->w_res({500:10},:)
       var(ie,n,2,:,:) = (/ ffu->$u0varname$(0,{500:10},:) /)
     end do
   end if

   if (casecode .eq. "hg") then
     fname = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
     do ie=0, ne-1
       ff = addfile(fname(ie),"r")
;       var(ie,n,iv,:,:) = ff->$varname(iv)$(0,{p_v},:,:)
     end do
   end if

   n = n + 1
 end do
 end do

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmp = u
   u(:,0) = tmp(:,nt-1)
   u(:,1:) = tmp(:,:nt-2)
   tmp = ush
   ush(:,0) = tmp(:,nt-1)
   ush(:,1:) = tmp(:,:nt-2)
   tmpv = var
   var(:,0,:,:,:) = tmpv(:,nt-1,:,:,:)
   var(:,1:,:,:,:) = tmpv(:,:nt-2,:,:,:)
   delete(tmp)
   delete(tmpv)
 end if

 uc = new((/ne,12/),float,1.e20)
 ushc = new((/ne,12/),float,1.e20)
 varc = new((/ne,12,3,nz,ny/),float,1.e20)
 if (casecode .eq. "hg") then
   ; use only 50 yrs
   do im=0, 11
     uc(:,im) = dim_avg_n(u(:,im:50*12-1:12), 1)
     ushc(:,im) = dim_avg_n(ush(:,im:50*12-1:12), 1)
     varc(:,im,:,:,:) = dim_avg_n(var(:,im:50*12-1:12,:,:,:), 1)
   end do
   if ( nma .eq. 3 .and. nt .gt. 50*12 ) then
     ; exclude the first (latest) Dec.
     uc(:,0) = dim_avg_n(u(:,12:50*12:12), 1)
     ushc(:,0) = dim_avg_n(ush(:,12:50*12:12), 1)
     varc(:,0,:,:,:) = dim_avg_n(var(:,12:50*12:12,:,:,:), 1)
   end if 
 else
   do im=0, 11
     uc(:,im) = dim_avg_n(u(:,im::12), 1)
     ushc(:,im) = dim_avg_n(ush(:,im::12), 1)
     varc(:,im,:,:,:) = dim_avg_n(var(:,im::12,:,:,:), 1)
   end do
 end if

 ua = u
 usha = ush
 vara = var
 do n=0, nt-1, 12
 do im=0, 11
   ua(:,n+im) = u(:,n+im) - uc(:,im)
   usha(:,n+im) = ush(:,n+im) - ushc(:,im)
   vara(:,n+im,:,:,:) = var(:,n+im,:,:,:) - varc(:,im,:,:,:)
 end do
 end do


 x = ua - 0.5*tofloat(xmin+xmax)
 y = (usha - 0.5*(ymin+ymax))*(tofloat(xmax-xmin)/tofloat(ymax-ymin))

 deg = where(x .ne. 0., atan(y/x)*180./3.141592, 90.)
 deg = where(abs(y) .gt. abs(x) .and. x*y .gt. 0.,  \
             90. - atan(x/y)*180./3.141592, deg)
 deg = where(abs(y) .gt. abs(x) .and. x*y .lt. 0.,  \
             -(90. - (atan(abs(x/y))*180./3.141592)), deg)
 deg = where(x .eq. 0. .and. y .lt. 0., -90., deg)
 deg = where(x .lt. 0., deg+180., deg)      ; -90 <= d < 270
 deg = 90. - deg                            ; -180 < d <= 180 (clockwise)
 deg = where(deg .lt. 0., deg+360., deg)    ;   0 <= d < 360

 deg00 = 45.
 deg0 = fspan(deg00,deg00+360.-(360./tofloat(nph)),nph)

 phs = where(deg .lt. deg0(0) .or. deg .ge. deg0(nph-1), 1, 999)
 do i=0, nph-2
   phs = where(deg .ge. deg0(i) .and. deg .lt. deg0(i+1), i+2, phs)
 end do


 var8m = new((/ne,12/nma,nph,3,nz,ny/),float)
 var8m = 0.
 std8m = var8m
 navg = new((/ne,12/nma,nph/),integer)
 navg = 0
 prob = new((/ne,12/nma,nph,3,nz,ny/),float)
 prob = 0.
 tmp = new((/3,nz,ny,nt/12*nma/),float)
 do ie=0, ne-1
 do im=0, 12/nma-1
 do ip=0, nph-1
   tmp = 0.
   do ma=0, nma-1
   do iy=0, nt/12-1
     imon = iy*12+im*nma+ma
     if (phs(ie,imon) .eq. ip+1) then
       tmp(:,:,:,iy*nma+ma) = vara(ie,imon,:,:,:)
       navg(ie,im,ip) = navg(ie,im,ip) + 1
     end if
   end do
   end do
   if (navg(ie,im,ip) .gt. 1) then
     var8m(ie,im,ip,:,:,:) = dim_sum(tmp)/tofloat(navg(ie,im,ip))
     std8m(ie,im,ip,:,:,:) = sqrt( ( dim_sum(tmp^2) -  \
                        tofloat(navg(ie,im,ip))*var8m(ie,im,ip,:,:,:)^2 )/  \
                      tofloat(navg(ie,im,ip)-1) )
     prob(ie,im,ip,:,:,1:ny-2) = (var8m(ie,im,ip,:,:,1:ny-2) - 0.)*  \
                            sqrt(navg(ie,im,ip))/std8m(ie,im,ip,:,:,1:ny-2)
     prob(ie,im,ip,:,:,1:ny-2) = student_t(prob(ie,im,ip,:,:,1:ny-2),navg(ie,im,ip)-1)
   end if
 end do
 end do
 end do

;--------------------------------------- cn level
 cnlev = (/-25,-20,-15,-10,-7,-5,-2,-1,-0.5,-0.2,0.,0.2,0.5,1,2,5,7,10,15,20,25/)

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(nph-1) = "u [m:S:-1:N:]"

 mstr = cexp
;p xstr = (/"","lat [deg]","","lat [deg]"/)
;p ystr = (/"z [km]","z [km]","",""/)
 xstr = (/"lat [deg]","lat [deg]","lat [deg]","lat [deg]"/)
 ystr = (/"z [km]","","",""/)

 wi = 50
 he = 30
 xf = new(nph,float)
 yf = new(nph,float)
 xf = 15
 yf = 10
;p xf(nph/2:) = xf(nph/2:) + (wi+8)
;p yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+12.5)
;p yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+12.5)
 xf = 15 + ispan(0,nph-1,1)*(wi+8)


var8m(0,:,:,0:1,:,:) = var8m(0,:,:,0:1,:,:)*100.
var8m(0,:,:,1,:,:) = var8m(0,:,:,1,:,:)*111.

do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, nph-1

 vw = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : var8m(0,im,ip,0:1,:,::4)
  "vfXArray" : lat(::4)
  "vfYArray" : z
;  "vfMissingUValueV" : 1.e15
;  "vfMissingVValueV" : 1.e15
 end create
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vw
  "vcRefMagnitudeF" : 10.0
  "vcRefLengthF"    : 0.03
  "vcRefAnnoOn"     : True
  "vcGlyphStyle"    : "FillArrow"
  "vcFillArrowHeadYF" : 0.18
  "vcFillArrowHeadXF" : 0.5
  "vcMinMagnitudeF" : 0.25

;  "vcRefAnnoArrowAngleF" : 90.
;  "vcRefAnnoFontHeightF" : 0.024
  "vcRefAnnoString1On"   : False
;  "vcRefAnnoString1"     : "1 mPa"
;  "vcRefAnnoString2On"   : False
  "vcRefAnnoString2"     : "10 cm s:S:-1:N:"
;  "vcRefAnnoPerimSpaceF" : 0.6
;  "vcRefAnnoOrthogonalPosF" : 0.15
;  "vcRefAnnoParallelPosF"   :  1. - nn*0.1

 end create
 set_size(vec1,xf(ip),yf(ip),wi,he)
 axis_range(vec1,-25,80,-999,30,False,False)

 contour1 = cnshadeplot(wks,lat,z,var8m(0,im,ip,2,:,:),  \
            cnlev,cnfil,cnlab(ip))
; cnaddline(contour1,wks,lat,z,prob(0,im,ip,0,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
; cnaddline(contour1,wks,lat,z,prob(0,im,ip,1,:,:)*100.,5.,0,"D",1,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf(ip),wi,he)
 axis_range(contour1,-25,80,-999,30,False,False)
 axis_str(contour1,"phase "+(ip+1),xstr(ip),ystr(ip))

 cnfinalize(contour1,wks)

 overlay(contour1,vec1)

 draw(contour1)
; cnfinalize(vec1,wks)
; draw(vec1)

 dr_txt(wks,navg(0,im,ip)+"",xf(ip)+wi-3,yf(ip)+he-2.5,22,8,0.)

end do

 frame(wks)
 delete(wks)

end do

end

