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

 varname = (/param0_varname_w,param0_varname_uv(0),param0_varname_uv(1),param0_varname_t,param0_varname_z/)
 p_v = (/70,100,200,500,850/)
; p_v = (/70,100,200/)
 if (casecode .eq. "hg") then
 if ( param0_lowlev ) then
   varname(2) = param0_varname_ll
;   delete(p_v)
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
   fname_v = (/param0_fname_w,param0_fname_uv(0),param0_fname_uv(1),param0_fname_t,param0_fname_z/)
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 vnamef = "tuv-clm"
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

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
;   fname1 = "$DATD/tem3d/ERA-I/2000/01/"+expname(0)+"."+fname_v(0)+".2000.01.nc"
   fname0 = "/data11/data-arch/ERA-I/monthly/2000/01/era-int."+fname_v(0)+".anal.pl.200001.mm.nc"
 end if
 if (casecode .eq. "hg") then
   fname0 = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname(0)+"a.pj_"+y1+"01.nc"
 end if
 ff = addfile(fname0,"r")
 axisname = ff->$varname(0)$!2
 lon = ff->$axisname$
 axisname = ff->$varname(0)$!1
 lat = ff->$axisname$({10:85})
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 var = new((/ne,nt,5,nz,ny,nx/),float,1.e20)
 l_firstmissing = False
 n = 0
 do yy=y1, y2
 do m=1, 12
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)
 do iv=0, 4

   if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
     fname_dir = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm
     fname = fname_dir+"/"+expname_v+"."+fname_v(iv)+".anal.pl."+yyyymm+".mm.nc"
     do ie=0, ne-1
       ff = addfile(fname(ie),"r")
       var(ie,n,iv,:,:,:) = (/ ff->$varname(iv)$({p_v},{10:85},:) /)
     end do
   end if

   if (casecode .eq. "hg") then
     fname = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
     do ie=0, ne-1
       ff = addfile(fname(ie),"r")
       var(ie,n,iv,:,:,:) = ff->$varname(iv)$(0,{p_v},{10:85},:)
     end do
   end if

 end do
   n = n + 1
 end do
 end do
 var(:,:,0,:,:,:) = var(:,:,0,:,:,:)*1.e3  ; [mm/s]

 var2 = new((/ne,nt,4,nz,ny,nx/),float,1.e20)

a_earth = 6371229.
deg2rad = tofloat(3.14159265358979323846/180.)
fcor = 2.*7.292116e-5*sin(lat*deg2rad)
itwodx = 1./( (lon(2) - lon(0))*deg2rad*a_earth*cos(lat*deg2rad) )
itwody = 1./( (lat(2) - lat(0))*deg2rad*a_earth )
gradx = new((/ne,nt,nz,ny,nx/),float,1.e20)
grady = new((/ne,nt,nz,ny,nx/),float,1.e20)
itwodx_m = conform_dims((/ne,nt,nz,ny,nx-2/),itwodx,3)
fcor_m = conform_dims((/ne,nt,nz,ny,nx/),fcor,3)
var2(:,:,3,:,:,1:nx-2) = (var(:,:,3,:,:,2:) - var(:,:,3,:,:,:nx-3))*itwodx_m
var2(:,:,3,:,:,0   ) = (var(:,:,3,:,:,1) - var(:,:,3,:,:,nx-1))*itwodx_m(:,:,:,:,0)
var2(:,:,3,:,:,nx-1) = (var(:,:,3,:,:,0) - var(:,:,3,:,:,nx-2))*itwodx_m(:,:,:,:,0)
var2(:,:,3,:,:,:) = var2(:,:,3,:,:,:)/fcor_m
var2(:,:,2,:,1:ny-2,:) = -(var(:,:,3,:,2:,:) - var(:,:,3,:,:ny-3,:))*itwody
var2(:,:,2,:,:,:) = var2(:,:,2,:,:,:)/fcor_m
var2(:,:,2,:,0   ,:) = var2(:,:,2,:,1   ,:)
var2(:,:,2,:,ny-1,:) = var2(:,:,2,:,ny-2,:)

 var2(:,:,0,:,:,:) = var(:,:,1,:,:,:) - var2(:,:,2,:,:,:)  ; u_a
 var2(:,:,1,:,:,:) = var(:,:,2,:,:,:) - var2(:,:,3,:,:,:)  ; v_a

;yh if ( l_firstmissing ) then
   var(:,0,:,:,:,:) = var@_FillValue
   var2(:,0,:,:,:,:) = var2@_FillValue
;yh end if

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmpv = var
   var(:,0,:,:,:,:) = tmpv(:,nt-1,:,:,:,:)
   var(:,1:,:,:,:,:) = tmpv(:,:nt-2,:,:,:,:)
   delete(tmpv)
 end if

 varc = new((/ne,12,5,nz,ny,nx/),float,1.e20)
 if (casecode .eq. "hg") then
   ; use only 50 yrs
   do im=0, 11
     varc(:,im,:,:,:,:) = dim_avg_n(var(:,im:50*12-1:12,:,:,:,:), 1)
   end do
   if ( nma .eq. 3 .and. nt .gt. 50*12 ) then
     ; exclude the first (latest) Dec.
     varc(:,0,:,:,:,:) = dim_avg_n(var(:,12:50*12:12,:,:,:,:), 1)
   end if 
 else
   do im=0, 11
     varc(:,im,:,:,:,:) = dim_avg_n(var(:,im::12,:,:,:,:), 1)
   end do
 end if

 var8m = new((/ne,12/nma,5,nz,ny,nx/),float)
 do im=0, 12/nma-1
   imon = im*nma
   var8m(:,im,:,:,:,:) = dim_avg_n(varc(:,imon:imon+nma-1,:,:,:,:),1)
 end do

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*0.5
; cnlev = fspan(-10,10,21)*0.2
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlevt = ispan(200,300,5)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(0) = "[mm s:S:-1:N:]"
 l_anno = new(nz,logical,False)
 l_anno(nz-1) = True

 mstr = cexp

 wi = 70 ;50
 he = 25
; xf = new(nph,float)
;p yf = new(nph,float)
 yf = new(nz,float)
 xf = 15
 yf = 10
;p xf(nph/2:) = xf(nph/2:) + (wi+8)
;p yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+10)
;p yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+10)
 xf = 15 ;+ ispan(0,nph-1,1)*(wi+4)
 yf = 10 + ispan(0,nz-1,1)*(he+8)


do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)
 colblu = NhlNewColor(wks,24./256,28./256,247./256)
 colred = NhlNewColor(wks,247./256,39./256,53./256)

ip = 0

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
  "mpMinLatF"      : 15. ;0.
  "mpMaxLatF"      : 80. ;50.
 end create
 set_size_mp(map,xf(ip),yf(k),wi,he)

 wind = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : var8m(0,im,1:2,k,::4,::8)
  "vfXArray" : lon(::8)
  "vfYArray" : lat(::4)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind
  "vcRefMagnitudeF" : 30.
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : l_anno(k)
  "vcGlyphStyle"    : "FillArrow"
  "vcFillArrowHeadYF" : 0.18
  "vcFillArrowHeadXF" : 0.5
;  "vcMinMagnitudeF" : 0.25
 end create
 set_size(vec1,xf(ip),yf(k),wi,he)

 contour1 = cnshadeplot(wks,lon,lat,var8m(0,im,0,k,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,var8m(0,im,3,k,:,:),cnlevt,0,"D",0,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf(k),wi,he)

; overlay(map,contour1)
 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)

 draw(map)

end do

 frame(wks)
 delete(wks)

end do

end

