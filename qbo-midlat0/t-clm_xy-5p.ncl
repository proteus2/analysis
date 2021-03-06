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

 varname = param0_varname_t
 p_v = (/70,100,200,500,850/)
 if (casecode .eq. "hg") then
 if ( param0_lowlev ) then
   varname = param0_varname_ll
   delete(p_v)
;   p_v = param0_z_v_ll
 end if
 end if

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 nma = param0_nma
; nph = param0_nph

 xmin = param0_umin_p1
 xmax = param0_umax_p1
 ymin = param0_umin_p2
 ymax = param0_umax_p2

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = param0_fname_t
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 vnamef = "t-clm"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_xy-5p_"+nma+"m_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   if ( param0_lowlev ) then
     foutname = vnamef+"_xy-z"+(p_v/1.e3)+"_"+nma+"m_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 else
   if (expname(0) .eq. "era-int_f") then
     foutname = vnamef+"_xy-5p_"+nma+"m_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_xy-5p_"+nma+"m_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12
 nz = dimsizes(p_v)
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

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname = "/data11/data-arch/ERA-I/2000/01/era-int."+fname_v+".anal.00.pl.200001.nc"
 end if
 if (casecode .eq. "hg") then
   fname = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname(0)+"a.pj_"+y1+"01.nc"
 end if
 ff = addfile(fname,"r")
 axisname = ff->$varname$!3
 lon = ff->$axisname$
 axisname = ff->$varname$!2
 lat = ff->$axisname$
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 var = new((/ne,nt,nz,ny,nx/),float,1.e20)
 n = 0
 do yy=y1, y2
 do m=1, 12
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)

   if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
     fname_dir = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm
     fname00 = fname_dir+"/"+  \
               expname_v+"."+fname_v+".anal.00.pl."+yyyymm+".mm.nc"
     fname06 = fname_dir+"/"+  \
               expname_v+"."+fname_v+".anal.06.pl."+yyyymm+".mm.nc"
     fname12 = fname_dir+"/"+  \
               expname_v+"."+fname_v+".anal.12.pl."+yyyymm+".mm.nc"
     fname18 = fname_dir+"/"+  \
               expname_v+"."+fname_v+".anal.18.pl."+yyyymm+".mm.nc"
     do ie=0, ne-1
       ff00 = addfile(fname00(ie),"r")
       ff06 = addfile(fname06(ie),"r")
       ff12 = addfile(fname12(ie),"r")
       ff18 = addfile(fname18(ie),"r")
       var(ie,n,:,:,:) = ( ff00->$varname$({p_v},:,:) + \
                         ff06->$varname$({p_v},:,:) + \
                         ff12->$varname$({p_v},:,:) + \
                         ff18->$varname$({p_v},:,:) )*0.25
     end do
   end if

   if (casecode .eq. "hg") then
     fname = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
     do ie=0, ne-1
       ff = addfile(fname(ie),"r")
       var(ie,n,:,:,:) = ff->$varname$(0,{p_v},:,:)
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
 varc = new((/ne,12,nz,ny,nx/),float,1.e20)
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


 var8m = new((/ne,12/nma,nz,ny,nx/),float)
 do ie=0, ne-1
 do im=0, 12/nma-1
   tmp = 0.
   do ma=0, nma-1
   do iy=0, nt/12-1
     imon = iy*12+im*nma+ma
   end do
   end do
   var8m(ie,im,ip,:,:,:) = dim_sum(tmp)/tofloat(navg(ie,im,ip))
 end do
 end do

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*0.5
 if (nph .eq. 4) then
   cnlev = fspan(-10,10,21)*0.2
 end if
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(nph-1) = "T [K]"

 mstr = cexp

 wi = 50
 he = 25
 xf = new(nph,float)
;p yf = new(nph,float)
 yf = new(nz,float)
 xf = 15
 yf = 10
;p xf(nph/2:) = xf(nph/2:) + (wi+8)
;p yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+10)
;p yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+10)
 xf = 15 + ispan(0,nph-1,1)*(wi+4)
 yf = 10 + ispan(0,nz-1,1)*(he+8)


do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, nph-1
do k=0, nz-1

 map = create "mp" mapPlotClass wks
;  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 125.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"
  "mpGeophysicalLineThicknessF" : 0.8

;  "mpShapeMode"     : "FixedAspectFitBB"
  "mpShapeMode"     : "FreeAspect"
  "mpGeophysicalLineColor" : colind ;"black"  ;"Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 0. ;85.
  "mpMaxLonF"      : 360. ;165.
  "mpMinLatF"      : -25. ;0.
  "mpMaxLatF"      : 80. ;50.
 end create
 set_size_mp(map,xf(ip),yf(k),wi,he)

 contour1 = cnshadeplot(wks,lon,lat,var8m(0,im,ip,k,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,prob(0,im,ip,k,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf(k),wi,he)
if (k .eq. 0) then
 axis_str(contour1,"phase "+(ip+1),"","")
end if

 overlay(map,contour1)
; overlay(map,vec1)
; NhlAddOverlay(map,contour1,map)

 draw(map)

if (k .eq. 0) then
 dr_txt(wks,navg(0,im,ip)+"",xf(ip)+wi-3,yf(k)+he-2.5,22,8,0.)
end if

end do
end do

 frame(wks)
 delete(wks)

end do

end

