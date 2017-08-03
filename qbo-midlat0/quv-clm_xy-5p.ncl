load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
load "./const_glob.ncl"

begin

;--------------------------------------- parameters
 casecode = "ra"
; casecode = "hg"

 if (casecode .eq. "ra") then
   load "./param-ra-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ;  "CGW","CTL"
 expname = param0_expname  ; "tc","tk"

 varname = (/"Q",param0_varname_z/)
 p_v = (/70,100,200,500,850/)
; p_v = (/70,100,200/)
 if (casecode .eq. "hg") then
 if ( param0_lowlev ) then
   varname(3) = param0_varname_ll
;   delete(p_v)
;   p_v = param0_z_v_ll
 end if
 end if

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = (/"tem3d_qg_nonst-p86_xyp",param0_fname_z/)
 end if
 
 vnamef = "quv-clm"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_xy-5p_"+nma+"m_hg_"+y1+"-"+y2
   if ( param0_lowlev ) then
     foutname = vnamef+"_xy-z"+(p_v/1.e3)+"_"+nma+"m_hg_"+y1+"-"+y2
   end if
 else
   foutname = vnamef+"_xy-5p_"+nma+"m_"+cexp+"_"+y1+"-"+y2
 end if

;--------------------------------------- read data
 nvi = dimsizes(varname)
 nv = 4
 nt = (y2-y1+1)*12
 nz = dimsizes(p_v)
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 iv_q = 0
 iv_z = 1

 iv_ug = 2
 iv_vg = 3
 
 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname1 = "$DATD/tem3d/ERA-I/2000/01/"+expname+"."+fname_v(0)+".2000.01.nc"
   fname0 = "/data11/data-arch/ERA-I/monthly/2000/01/era-int."+fname_v(1)+".anal.pl.200001.mm.nc"
 end if
 if (casecode .eq. "hg") then
   fname0 = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+y1+"01.nc"
 end if
 ff0 = addfile(fname0,"r")
 axisname = ff0->$varname(1)$!2
 lon = ff0->$axisname$
 axisname = ff0->$varname(1)$!1
 lat = ff0->$axisname$({10:85})
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 ff1 = addfile(fname1,"r")
 axisname = ff1->$varname(0)$!2
 lon2 = ff1->$(ff1->$varname(0)$!2)$
 axisname = ff1->$varname(0)$!1
 lat2 = ff1->$axisname$({10:85})

 if ( any(lat .ne. lat2) .or. any(lon .ne. lon2) ) then
   print("Grids are different.")
   exit
 end if

 var = new((/nv,nt,nz,ny,nx/),float,1.e20)
 l_firstmissing = False

 fname = new(nt,string)
 do iv=0, nvi-1

   n = 0
   do yy=y1, y2
   do m=1, 12
     mm = sprinti("%2.2i",m)
     yyyymm = yy+sprinti("%2.2i",m)
     if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
       fname(n) = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm+"/"+  \
                  expname_v+"."+fname_v(iv)+".anal.pl."+yyyymm+".mm.nc"
       if (iv .eq. iv_q) then
         fname(n) = "$DATD/tem3d/"+cexp_v+"/"+yy+"/"+mm+"/"+  \
                    expname_v+".tem3d_qg_nonst-p86_xyp."+yy+"."+mm+".nc"
       end if
     end if
     if (casecode .eq. "hg") then
       fname(n) = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
     end if
     if ( yyyymm .eq. "197901" .and. .not. isfilepresent(fname(n)) ) then
       l_firstmissing = True
     end if
     n = n + 1
   end do
   end do

   if ( l_firstmissing ) then
     fname(0) = fname(1)
   end if

   ff = addfiles(fname,"r")
   if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
     ListSetType(ff, "join")
   end if
   if (casecode .eq. "hg") then
     ListSetType(ff, "cat")
   end if
   var(iv,:,:,:,:) = (/ ff[:]->$varname(iv)$(:,{p_v},{10:85},:) /)
   delete(ff)

 end do

 fcor = 2.*ome_earth*sin(lat*deg2rad)
 itwodx = 1./( (lon(2) - lon(0))*deg2rad*a_earth*cos(lat*deg2rad) )
 itwody = 1./( (lat(2) - lat(0))*deg2rad*a_earth )
 itwodx_m = conform_dims((/nt,nz,ny,nx-2/),itwodx,2)
 fcor_m = conform_dims((/nt,nz,ny,nx/),fcor,2)

 ; u_g
 var(iv_ug,:,:,1:ny-2,:) = -(var(iv_z,:,:,2:,:) - var(iv_z,:,:,:ny-3,:))*itwody
 var(iv_ug,:,:,:,:) = var(iv_ug,:,:,:,:)/fcor_m
 var(iv_ug,:,:,0   ,:) = var(iv_ug,:,:,1   ,:)
 var(iv_ug,:,:,ny-1,:) = var(iv_ug,:,:,ny-2,:)
 ; v_g
 var(iv_vg,:,:,:,1:nx-2) = (var(iv_z,:,:,:,2:) - var(iv_z,:,:,:,:nx-3))*itwodx_m
 var(iv_vg,:,:,:,0   ) = (var(iv_z,:,:,:,1) - var(iv_z,:,:,:,nx-1))*  \
                         itwodx_m(:,:,:,0)
 var(iv_vg,:,:,:,nx-1) = (var(iv_z,:,:,:,0) - var(iv_z,:,:,:,nx-2))*  \
                         itwodx_m(:,:,:,0)
 var(iv_vg,:,:,:,:) = var(iv_vg,:,:,:,:)/fcor_m
 

 if ( l_firstmissing ) then
   var(:,0,:,:,:) = var@_FillValue
 end if

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmpv = var
   var(:,0,:,:,:) = tmpv(:,nt-1,:,:,:)
   var(:,1:,:,:,:) = tmpv(:,:nt-2,:,:,:)
   delete(tmpv)
 end if

 varc = new((/nv,12,nz,ny,nx/),float,1.e20)
 if (casecode .eq. "hg") then
   ; use only 50 yrs
   do im=0, 11
     varc(:,im,:,:,:) = dim_avg_n(var(:,im:50*12-1:12,:,:,:), 1)
   end do
   if ( nma .eq. 3 .and. nt .gt. 50*12 ) then
     ; exclude the first (latest) Dec.
     varc(:,0,:,:,:) = dim_avg_n(var(:,12:50*12:12,:,:,:), 1)
   end if 
 else
   do im=0, 11
     varc(:,im,:,:,:) = dim_avg_n(var(:,im::12,:,:,:), 1)
   end do
 end if

 var8m = new((/nv,12/nma,nz,ny,nx/),float)
 do im=0, 12/nma-1
   imon = im*nma
   var8m(:,im,:,:,:) = dim_avg_n(varc(:,imon:imon+nma-1,:,:,:),1)
 end do

 fcoef = ezfftf(var8m)
 fcoef(:,:,:,:,:,15:) = 0.
 var8m = ezfftb(fcoef, fcoef@xbar)
 var8m = wgt_runave_n_Wrap(var8m,(/1.,2.,1./),0,3)

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*1.0
; cnlev = fspan(-10,10,21)*0.2
 cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlevq = ispan(0,20,1)
 dcnq = where( p_v .lt. 500., 0.25e-4, 0.25e-4 )

;--------------------------------------- contour
 cnlab = new(nz,string)
 cnlab = ""
 cnlab(nz-1) = "[mm s:S:-1:N:]"
 l_anno = new(nz,logical,False)
 l_anno(nz-1) = True

 mstr = cexp

 latrng_pl = 65.
 scale_pl = 0.3
 he = latrng_pl*scale_pl
 wi = 360.*scale_pl/sqrt(2.)
 yf = new(nz,float)
 xf = 15 + ispan(0,1,1)*(wi+4)
 yf = 10 + ispan(0,nz-1,1)*(he+8)


do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)
 colgrn = NhlNewColor(wks,0,1,0)
 colblu = NhlNewColor(wks,24./256,28./256,247./256)
 colred = NhlNewColor(wks,247./256,39./256,53./256)

do ip=0, 0 ;1

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
  "mpGeophysicalLineColor" : colgrn
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 0. ;85.
  "mpMaxLonF"      : 360. ;165.
  "mpMinLatF"      : 15. ;0.
  "mpMaxLatF"      : 80. ;50.
 end create
 set_size_mp(map,xf(ip),yf(k),wi,he)

 if (ip .eq. 0) then
   wind = create "f1" vectorFieldClass defaultapp
    "vfDataArray" : var8m((/iv_ug,iv_vg/),im,k,::5,::10)
    "vfXArray" : lon(::10)
    "vfYArray" : lat(::5)
    "vfMissingUValueV" : var8m@_FillValue
    "vfMissingVValueV" : var8m@_FillValue
   end create
 end if
 if (ip .eq. 1) then
   wind = create "f1" vectorFieldClass defaultapp
    "vfDataArray" : var8m((/iv_ua,iv_va/),im,k,::5,::10)
    "vfXArray" : lon(::10)
    "vfYArray" : lat(::5)
    "vfMissingUValueV" : var8m@_FillValue
    "vfMissingVValueV" : var8m@_FillValue
   end create
 end if
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind
  "vcRefMagnitudeF" : 30.
;  "vcMinMagnitudeF" : 1.5
  "vcRefLengthF"    : 0.04
  "vcRefAnnoOn"     : l_anno(k)
  "vcRefAnnoString2On" : False
  "vcGlyphStyle"    : "FillArrow"
  "vcFillArrowWidthF" : 0.05
  "vcFillArrowHeadYF" : 0.18
  "vcFillArrowHeadXF" : 0.5
 end create
 if (ip .eq. 1) then
  setvalues vec1
   "vcRefMagnitudeF" : 3.
;   "vcMinMagnitudeF" : 0.15
  end setvalues
 end if
 set_size(vec1,xf(ip),yf(k),wi,he)

; qm = conform_dims( (/ny,nx/), dim_avg( var8m(iv_q,im,k,:,:), 0 )
; contour1 = cnshadeplot(wks,lon,lat,var8m(iv_q,im,k,:,:) - qm,cnlevq,cnfil,cnlab(k))
; contour1 = cnshadeplot(wks,lon,lat,var8m(iv_q,im,k,:,:),cnlevq*dcnq,cnfil,cnlab(k))
; cnaddline(contour1,wks,lon,lat,qm,cnlev,0,"D",0,"T",1.0,"C",1)
 contour1 = cnlineplot(wks,lon,lat,var8m(iv_q,im,k,:,:),cnlevq*dcnq(k),0)
 cnline_att(contour1,"D",0,"T",1.2,"C",colred)
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

end do  ; k

end do  ; ip

 frame(wks)
 delete(wks)

end do

end

